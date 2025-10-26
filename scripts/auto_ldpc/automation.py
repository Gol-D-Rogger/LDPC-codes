"""
High-level automation driver for LDPC matrix exploration.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

CURRENT_FILE = Path(__file__).resolve()
WORKSPACE_PARENT = CURRENT_FILE.parents[2]
if str(WORKSPACE_PARENT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_PARENT))

from scripts import submit_ldpc_bsub as legacy_submit

from scripts.auto_ldpc.commands import CommandRunner
from scripts.auto_ldpc.config import ConfigBundle, load_bundle
from scripts.auto_ldpc.jobs import SimulationJob, build_matrix_job_name, format_inner_command, submit_jobs
from scripts.auto_ldpc.log_parser import scan_logs
from scripts.auto_ldpc.state import StateStore
from scripts.auto_ldpc.xlsx_export import export_records
from scripts.auto_ldpc.plotting import plot_top_results


class AutomationRunner:
    def __init__(self, bundle: ConfigBundle, dry_run: bool = False) -> None:
        self.bundle = bundle
        self.cfg = bundle.raw
        self.paths = bundle.paths
        system_cfg = self.cfg.get("system", {})
        self.runner = CommandRunner(system_cfg.get("command_prefix"), dry_run=dry_run)
        database_cfg = self.cfg.get("database", {})
        db_path = self.bundle.expand_path(database_cfg.get("path"))
        if db_path is None:
            base = self.paths.output_root or self.paths.workspace_root
            db_path = (base / "state/state.db").resolve()
        db_path.parent.mkdir(parents=True, exist_ok=True)
        self.state = StateStore(db_path)
        self.dry_run = dry_run

        self.anchor_snr = float(self.cfg["simulation"]["anchor_snr"])
        self.top_n = int(self.cfg["selection"]["top_n"])
        self.errinj_threshold = float(self.cfg["selection"]["err_inj_threshold"])
        self.snr_list = [float(snr) for snr in self.cfg["simulation"]["snr_list"]]
        self.deep_snr_list = [float(snr) for snr in self.cfg["simulation"]["deep_snr_list"]]
        self.sleep_seconds = int(self.cfg["limits"]["sleep_seconds"])
        self.max_matrices = int(self.cfg["limits"]["max_matrices"])
        stop_flag_value = self.cfg["limits"]["stop_flag_file"]
        self.stop_flag_path = self.bundle.expand_path(stop_flag_value, base=self.paths.output_root or self.paths.workspace_root)

        self.errinj_cfg = self.cfg["err_inj"]
        self.sim_cfg = self.cfg["simulation"]
        self.gen_cfg = self.cfg["generation"]
        self.plot_cfg = self.cfg.get("plot", {})

        log_parser_cfg = self.cfg.get("log_parser", {})
        self.awgn_log_cfg = log_parser_cfg.get("awgn", {})
        self.errinj_log_cfg = log_parser_cfg.get("err_inj", {})
        self.awgn_stat_keys = self.awgn_log_cfg.get("stat_keys", ["LDPC FER"])
        self.awgn_fer_key = self.awgn_log_cfg.get("fer_key", "LDPC FER")
        self.errinj_stat_keys = self.errinj_log_cfg.get("stat_keys", ["LDPC FER"])
        self.errinj_fer_key = self.errinj_log_cfg.get("fer_key", "LDPC FER")

        # Determine matrix storage root via legacy helper
        default_args = self._build_generation_args(count=0)
        default_paths = legacy_submit.prepare_paths(default_args)
        self.sim_matrix_root = default_paths.sim_matrix_root

    def close(self) -> None:
        self.state.close()

    # --- Log parsing -----------------------------------------------------

    def refresh_awgn_logs(self) -> None:
        awgn_root = self.paths.sim_out_root / "awgn"
        if not awgn_root.exists():
            print(f"  [警告] AWGN 日志目录不存在: {awgn_root}")
            return
            
        results = scan_logs(
            log_root=awgn_root,
            model="AWGN",
            stat_keys=self.awgn_stat_keys,
            fer_key=self.awgn_fer_key,
        )
        
        parsed_count = 0
        anchor_count = 0
        snr_values_seen = set()
        for result in results:
            if result.fer_value is None:
                continue
            snr_values_seen.add(result.snr)
            is_anchor = math.isclose(result.snr, self.anchor_snr, rel_tol=1e-6, abs_tol=1e-6)
            if not is_anchor and abs(result.snr - self.anchor_snr) < 0.1:
                print(f"  [调试] matrix{result.matrix_id} SNR={result.snr} 不匹配 anchor={self.anchor_snr}, 差值={abs(result.snr - self.anchor_snr)}, 文件={result.log_path.name}")
            self.state.upsert_matrix(result.matrix_id)
            self.state.record_awgn_run(
                matrix_id=result.matrix_id,
                fer=result.fer_value,
                metrics=result.metrics,
                snr=result.snr,
                log_path=result.log_path,
                is_anchor=is_anchor,
            )
            parsed_count += 1
            if is_anchor:
                anchor_count += 1
        
        print(f"  解析到 {parsed_count} 条 AWGN 结果 (其中 {anchor_count} 条为 anchor SNR={self.anchor_snr})")
        if snr_values_seen:
            print(f"  检测到的 SNR 值: {sorted(snr_values_seen)}")

    def refresh_errinj_logs(self) -> None:
        err_root = self.paths.sim_out_root / "err_inj"
        if not err_root.exists():
            print(f"  [警告] ERR_INJ 日志目录不存在: {err_root}")
            return
            
        results = scan_logs(
            log_root=err_root,
            model="ERR_INJ",
            stat_keys=self.errinj_stat_keys,
            fer_key=self.errinj_fer_key,
            snr_override=0.0,
        )
        
        parsed_count = 0
        passed_count = 0
        for result in results:
            if result.fer_value is None:
                continue
            passes = result.fer_value <= self.errinj_threshold
            self.state.upsert_matrix(result.matrix_id)
            self.state.update_err_inj_result(
                matrix_id=result.matrix_id,
                fer=result.fer_value,
                metrics=result.metrics,
                log_path=result.log_path,
                passes_threshold=passes,
            )
            parsed_count += 1
            if passes:
                passed_count += 1
        
        print(f"  解析到 {parsed_count} 条 ERR_INJ 结果 (其中 {passed_count} 条通过阈值 ≤{self.errinj_threshold})")

    # --- Job submission helpers -----------------------------------------

    def _awgn_log_path(self, matrix_id: str, snr: float) -> Path:
        dir_path = self.paths.sim_out_root / "awgn" / f"matrix{matrix_id}"
        dir_path.mkdir(parents=True, exist_ok=True)
        return dir_path / f"snr{snr}.log"

    def _errinj_log_path(self, matrix_id: str) -> Path:
        dir_path = self.paths.sim_out_root / "err_inj" / f"matrix{matrix_id}"
        dir_path.mkdir(parents=True, exist_ok=True)
        return dir_path / "err_inj.log"

    def submit_awgn_runs(
        self,
        matrix_ids: Sequence[str],
        snr_list: Sequence[float],
        matrix_dir: Optional[Path] = None,
        job_prefix: Optional[str] = None,
    ) -> None:
        if matrix_dir is None:
            matrix_dir = self.sim_matrix_root
        jobs: List[SimulationJob] = []
        queue = self.sim_cfg.get("queue", "regr_q")
        prefix = job_prefix or self.sim_cfg.get("job_prefix", "ldpc_auto")
        for matrix_id in matrix_ids:
            for snr in snr_list:
                log_path = self._awgn_log_path(matrix_id, snr)
                err_path = log_path.with_suffix(".err")
                command = format_inner_command(
                    exec_path=self.paths.exec_path,
                    config_path=self.paths.config_path,
                    model="AWGN",
                    snr=snr,
                    matrix_id=matrix_id,
                    matrix_dir=matrix_dir,
                )
                job_name = build_matrix_job_name(prefix, matrix_id, "AWGN", snr)
                jobs.append(
                    SimulationJob(
                        matrix_id=matrix_id,
                        snr=snr,
                        model="AWGN",
                        queue=queue,
                        job_name=job_name,
                        log_path=log_path,
                        err_path=err_path,
                        command=command,
                    )
                )
                if not self.dry_run:
                    self.state.increment_total_runs(matrix_id)
        submit_jobs(self.runner, jobs)

    def submit_errinj_runs(self, matrix_ids: Sequence[str]) -> None:
        jobs: List[SimulationJob] = []
        queue = self.errinj_cfg.get("queue", self.sim_cfg.get("queue", "regr_q"))
        prefix = self.errinj_cfg.get("job_prefix", "ldpc_errinj")
        template = self.errinj_cfg.get("command")
        config_path = Path(self.errinj_cfg.get("config", "config/another.cnfg"))
        if not config_path.is_absolute():
            config_path = (self.paths.workspace_root / config_path).resolve()
        exec_dir = self.paths.exec_path.parent
        for matrix_id in matrix_ids:
            log_path = self._errinj_log_path(matrix_id)
            err_path = log_path.with_suffix(".err")
            command_inner = template.format(
                config=str(config_path),
                matrix_id=matrix_id,
                matrix_dir=str(self.sim_matrix_root),
            )
            command = f'cd "{exec_dir}" && export LDPC_MATRIX_DIR="{self.sim_matrix_root}"; {command_inner}'
            job_name = build_matrix_job_name(prefix, matrix_id, "ERRINJ")
            jobs.append(
                SimulationJob(
                    matrix_id=matrix_id,
                    snr=0.0,
                    model="ERR_INJ",
                    queue=queue,
                    job_name=job_name,
                    log_path=log_path,
                    err_path=err_path,
                    command=command,
                )
            )
            if not self.dry_run:
                self.state.mark_err_inj_submitted(matrix_id, job_name)
            if not self.dry_run:
                self.state.increment_total_runs(matrix_id)
        submit_jobs(self.runner, jobs)

    def submit_deep_runs(self, matrix_ids: Sequence[str], matrix_dir: Optional[Path] = None) -> None:
        if not self.deep_snr_list:
            return
        for matrix_id in matrix_ids:
            self.submit_awgn_runs([matrix_id], self.deep_snr_list, matrix_dir, job_prefix="ldpc_deep")
            self.state.mark_deep_runs_submitted(matrix_id)

    # --- Matrix generation and scheduling --------------------------------

    def _build_generation_args(self, count: int) -> argparse.Namespace:
        legacy_out_dir = (self.paths.output_root or self.paths.sim_out_root) / "legacy_out"
        ns = argparse.Namespace(
            exec=str(self.paths.exec_path),
            config=str(self.paths.config_path),
            gen_job_name=None,
            matrix_dir=str(self.paths.matrix_dir) if self.paths.matrix_dir else None,
            out_dir=str(legacy_out_dir),
            queue=self.sim_cfg.get("queue", "regr_q"),
            job_prefix=self.sim_cfg.get("job_prefix", "ldpc_auto"),
            cwd=None,
            snr=None,
            snr_seq=None,
            build=False,
            gen_bin=str(self.paths.gen_bin),
            gen_n=self.gen_cfg.get("gen_n"),
            gen_k=self.gen_cfg.get("gen_k"),
            gen_count=count,
            gen_queue=self.gen_cfg.get("queue"),
            gen_mode=self.gen_cfg.get("mode", "lsf"),
            gen_out_base=self._resolve_gen_out_base(),
            rename_enable=not self.gen_cfg.get("disable_rename", False),
            no_rename=False,
            sim_matrix_dir=self._resolve_anchor_matrix_dir(),
            dry_run=self.dry_run,
        )
        return ns

    def generate_and_submit_batch(self, count: int) -> None:
        if count <= 0:
            return
        args = self._build_generation_args(count)
        paths = legacy_submit.prepare_paths(args)
        result = legacy_submit.generate_matrices(args, paths, args.dry_run)
        matrix_dir = result.matrix_dir
        matrix_ids = result.new_ids or []
        if not matrix_ids:
            return
        if result.sim_root:
            self.sim_matrix_root = result.sim_root
        self.state.ensure_batch_matrices(matrix_ids)
        self.submit_awgn_runs(matrix_ids, self.snr_list, result.sim_root)

    def _resolve_gen_out_base(self) -> Optional[str]:
        value = self.gen_cfg.get("gen_out_base")
        path = self.bundle.expand_path(value, base=self.paths.output_root or self.paths.workspace_root)
        if path is None:
            base = self.paths.output_root or self.paths.matrix_dir or self.paths.workspace_root
            path = base
        return str(path)

    def _resolve_anchor_matrix_dir(self) -> Optional[str]:
        value = self.gen_cfg.get("anchor_matrix_dir")
        path = self.bundle.expand_path(value, base=self.paths.matrix_dir or self.paths.output_root or self.paths.workspace_root)
        if path is None:
            return str(self.paths.matrix_dir) if self.paths.matrix_dir else None
        return str(path)

    # --- XLSX export -----------------------------------------------------

    def export_xlsx_if_enabled(self) -> Optional[Path]:
        xlsx_cfg = self.cfg["xlsx"]
        if not xlsx_cfg.get("enable", True):
            return None
        filename = xlsx_cfg.get("filename", "simulation_results.xlsx")
        output_path = self.paths.perf_dir / filename
        records = self.state.fetch_runs_for_xlsx()
        if not records:
            return None
        export_records(records, output_path)
        return output_path

    def plot_top_if_enabled(self, matrix_ids: Sequence[str]) -> None:
        if not self.plot_cfg.get("enable", True):
            return
        if not matrix_ids:
            return
        default_output = (self.paths.output_root or self.paths.perf_dir) / "best_plot.png"
        output_value = self.plot_cfg.get("output")
        output_path = self.bundle.expand_path(output_value, base=self.paths.perf_dir) if output_value else None
        if output_path is None:
            output_path = default_output
        baseline_cfg_raw = self.plot_cfg.get("baseline", {})
        baseline_cfg = dict(baseline_cfg_raw) if isinstance(baseline_cfg_raw, dict) else {}
        if "path" in baseline_cfg:
            baseline_path = self.bundle.expand_path(baseline_cfg.get("path"), base=self.paths.perf_dir)
            baseline_cfg["path"] = str(baseline_path) if baseline_path else None
        result = plot_top_results(
            state=self.state,
            matrix_ids=matrix_ids,
            output_path=output_path,
            workspace_root=self.paths.workspace_root,
            baseline_cfg=baseline_cfg,
        )
        if result:
            print(f"[INFO] 已生成最佳矩阵 SNR-FER 曲线: {result}")

    # --- Control flow ----------------------------------------------------

    def should_stop_generation(self) -> bool:
        """检查是否应该停止生成新矩阵"""
        if self.max_matrices and self.state.count_matrices() >= self.max_matrices:
            return True
        return False

    def should_exit(self) -> bool:
        """检查是否应该退出整个循环"""
        if self.stop_flag_path and self.stop_flag_path.exists():
            print(f"检测到停止标志文件 {self.stop_flag_path} ，自动化即将结束。")
            return True
        return False

    def run_iteration(self) -> None:
        print(f"[INFO] 开始新的迭代循环...")
        
        print(f"[INFO] 刷新 AWGN 日志...")
        self.refresh_awgn_logs()
        
        print(f"[INFO] 刷新 ERR_INJ 日志...")
        self.refresh_errinj_logs()

        print("[INFO] 获取 Top-N 矩阵...")
        top_records = self.state.fetch_top_matrices(self.top_n)
        print(f"[INFO] 获取到 Top-{self.top_n} 矩阵: {len(top_records)} 个")
        
        # 显示 Top-N 矩阵详情
        for i, rec in enumerate(top_records, 1):
            print(f"  {i}. matrix_id={rec.matrix_id}, anchor_fer={rec.anchor_fer}, "
                  f"passes_threshold={rec.passes_threshold}, deep_submitted={rec.deep_runs_submitted}")
        
        pending_errinj = [rec.matrix_id for rec in top_records if rec.err_inj_status in ("pending", "failed")]
        if pending_errinj:
            print(f"[INFO] 提交 ERR_INJ 任务: {len(pending_errinj)} 个矩阵 - {pending_errinj}")
            self.submit_errinj_runs(pending_errinj)
        else:
            print(f"[INFO] 无需提交 ERR_INJ 任务")

        ready_for_deep = [rec.matrix_id for rec in top_records if rec.passes_threshold and not rec.deep_runs_submitted]
        if ready_for_deep:
            print(f"[INFO] 提交深挖任务: {len(ready_for_deep)} 个矩阵 - {ready_for_deep}")
            print(f"[DEBUG] 深挖条件详情:")
            for rec in top_records:
                if rec.matrix_id in ready_for_deep:
                    print(f"  - matrix_id={rec.matrix_id}: passes_threshold={rec.passes_threshold}, "
                          f"deep_runs_submitted={rec.deep_runs_submitted}, err_inj_fer={rec.err_inj_fer}")
            self.submit_deep_runs(ready_for_deep)
        else:
            print(f"[INFO] 无需提交深挖任务")
            print(f"[DEBUG] 深挖跳过原因:")
            for rec in top_records:
                reason = []
                if not rec.passes_threshold:
                    reason.append(f"未通过ERR_INJ阈值(fer={rec.err_inj_fer})")
                if rec.deep_runs_submitted:
                    reason.append("已提交深挖")
                if reason:
                    print(f"  - matrix_id={rec.matrix_id}: {', '.join(reason)}")

        batch_size = int(self.gen_cfg.get("batch_size") or 0)
        current_count = self.state.count_matrices()
        print(f"[INFO] 当前矩阵总数: {current_count}, batch_size: {batch_size}, max_matrices: {self.max_matrices}")
        
        if not self.should_stop_generation() and batch_size > 0:
            print(f"[INFO] 生成并提交新批次: {batch_size} 个矩阵")
            self.generate_and_submit_batch(batch_size)
        elif self.should_stop_generation() and batch_size > 0:
            print(f"[INFO] 已达到配置的最大矩阵数量，停止生成新矩阵，但继续处理现有任务...")
        elif batch_size == 0:
            print(f"[INFO] batch_size=0，不生成新矩阵，仅处理现有任务")

        matrix_ids_for_plot = [rec.matrix_id for rec in top_records if rec.anchor_fer is not None]
        print(f"[INFO] 导出 XLSX 和生成图表...")
        self.export_xlsx_if_enabled()
        self.plot_top_if_enabled(matrix_ids_for_plot)
        print(f"[INFO] 本轮迭代完成\n")

    def run(self, loop: bool = True) -> None:
        print(f"[INFO] ========== 自动化脚本启动 ==========")
        print(f"[INFO] 循环模式: {'是' if loop else '否'}")
        print(f"[INFO] 工作空间: {self.paths.workspace_root}")
        print(f"[INFO] 输出目录: {self.paths.output_root}")
        print(f"[INFO] 状态数据库: {self.state.db_path}")
        print(f"[INFO] 睡眠间隔: {self.sleep_seconds} 秒")
        print(f"[INFO] =====================================\n")
        
        try:
            iteration_count = 0
            while True:
                iteration_count += 1
                print(f"\n{'='*60}")
                print(f"第 {iteration_count} 轮迭代")
                print(f"{'='*60}")
                
                self.run_iteration()
                
                if not loop:
                    print(f"[INFO] 单次迭代模式，退出")
                    break
                    
                if self.should_exit():
                    break
                
                print(f"[INFO] 等待 {self.sleep_seconds} 秒后开始下一轮...")
                time.sleep(self.sleep_seconds)
        finally:
            print(f"\n[INFO] 自动化脚本结束，关闭数据库连接")
            self.close()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LDPC 自动化探索脚本")
    parser.add_argument("--config", required=True, help="配置文件路径 (JSON 或 YAML)")
    parser.add_argument("--dry-run", action="store_true", help="仅打印命令，不实际提交。")
    parser.add_argument("--once", action="store_true", help="只执行一次迭代。")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    config_path = Path(args.config).resolve()
    bundle = load_bundle(config_path)
    runner = AutomationRunner(bundle, dry_run=args.dry_run)
    runner.run(loop=not args.once)
    return 0


if __name__ == "__main__":
    sys.exit(main())
