"""
Local automation driver for LDPC matrix exploration.

此版本与 scripts/auto_ldpc/automation.py 保持同样的流程，只是矩阵生成、
AWGN 仿真、ERR_INJ 及深挖全部在本地直接执行，便于无需 LSF/bsub 的环境
进行调试。其他功能（状态数据库、日志解析、绘图等）与原脚本一致。
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import sys
import time
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

CURRENT_FILE = Path(__file__).resolve()
WORKSPACE_PARENT = CURRENT_FILE.parents[2]
if str(WORKSPACE_PARENT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_PARENT))

from scripts.auto_ldpc.commands import CommandRunner
from scripts.auto_ldpc.config import ConfigBundle, load_bundle
from scripts.auto_ldpc.log_parser import scan_logs
from scripts.auto_ldpc.plotting import plot_top_results
from scripts.auto_ldpc.state import StateStore
from scripts.auto_ldpc.xlsx_export import export_records


class LocalAutomationRunner:
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

        self.sim_matrix_root = self._resolve_sim_matrix_root()

    # ------------------------------------------------------------------ #
    # 基础工具                                                           #
    # ------------------------------------------------------------------ #

    def _ensure_dir(self, path: Path) -> None:
        if not self.dry_run:
            path.mkdir(parents=True, exist_ok=True)

    def _matrix_dimensions(self) -> Tuple[int, int, int]:
        n = int(self.gen_cfg.get("gen_n"))
        k = int(self.gen_cfg.get("gen_k"))
        mk = n - k
        return n, k, mk

    def _gen_out_base(self) -> Path:
        value = self.gen_cfg.get("gen_out_base")
        base = self.bundle.expand_path(value, base=self.paths.output_root or self.paths.workspace_root)
        if base is None:
            base = (self.paths.output_root or self.paths.workspace_root)
        return Path(base)

    def _matrix_dir(self) -> Path:
        n, _, mk = self._matrix_dimensions()
        return self._gen_out_base() / f"{mk}x{n}" / "matrix"

    def _mask_dir(self) -> Path:
        n, _, mk = self._matrix_dimensions()
        return self._gen_out_base() / f"{mk}x{n}" / "mask_matrix"

    def _cycle_dir(self) -> Path:
        n, _, mk = self._matrix_dimensions()
        return self._gen_out_base() / f"{mk}x{n}" / "cycle_record"

    def _resolve_sim_matrix_root(self) -> Path:
        value = self.gen_cfg.get("anchor_matrix_dir")
        path = self.bundle.expand_path(value, base=self.paths.matrix_dir or self._gen_out_base())
        if path is None:
            return self._gen_out_base()
        return Path(path)

    def _list_matrix_files(self, directory: Path) -> List[Path]:
        if not directory.exists():
            return []
        return sorted(directory.glob("*_QC_H_*.txt"))

    def _extract_new_ids(self, matrix_dir: Path, before: List[Path], after: List[Path]) -> List[str]:
        before_set = {p.resolve() for p in before}
        new_files = [p.resolve() for p in after if p.resolve() not in before_set]
        new_ids: List[str] = []
        pattern_new = re.compile(r"_([0-9]+)\.txt$")
        pattern_old = re.compile(r"_([0-9]+)_([0-9]+)\.txt$")
        for file_path in new_files:
            match = pattern_new.search(file_path.name)
            if match:
                new_ids.append(match.group(1))
        if not new_ids:
            for file_path in new_files:
                match = pattern_old.search(file_path.name)
                if match:
                    new_ids.append(f"{match.group(1)}_{match.group(2)}")
        return new_ids

    # ------------------------------------------------------------------ #
    # 矩阵生成                                                           #
    # ------------------------------------------------------------------ #

    def generate_matrices(self, count: int) -> List[str]:
        print(f"[generate_matrices] 请求生成 {count} 个矩阵")
        if count <= 0:
            print("[generate_matrices] count <= 0, 返回空列表")
            return []

        n, k, mk = self._matrix_dimensions()
        gen_out_base = self._gen_out_base()
        matrix_dir = self._matrix_dir()
        print(f"[generate_matrices] 矩阵目录: {matrix_dir}")
        self._ensure_dir(matrix_dir)
        self._ensure_dir(self._mask_dir())
        self._ensure_dir(self._cycle_dir())

        files_before = self._list_matrix_files(matrix_dir)
        print(f"[generate_matrices] 生成前文件数: {len(files_before)}")

        env = os.environ.copy()
        env["GENLDPC_OUT_DIR"] = str(gen_out_base)

        for idx in range(1, count + 1):
            cmd = [
                str(self.paths.gen_bin),
                str(n),
                str(k),
                str(idx),
                str(gen_out_base),
            ]
            print(f"[generate_matrices] 运行生成命令 {idx}/{count}: {' '.join(cmd)}")
            self.runner.run(
                cmd,
                capture_output=False,
                cwd=str(self.paths.gen_bin.parent),
                env=env,
            )

        if not self.gen_cfg.get("disable_rename", False):
            rename_cmd = [
                "bash",
                str(self.paths.rename_script),
                str(mk),
                str(n),
                "--apply",
                "--base-dir",
                str(gen_out_base),
            ]
            print(f"[generate_matrices] 运行重命名: {' '.join(rename_cmd)}")
            self.runner.run(rename_cmd, cwd=str(self.paths.rename_script.parent))

        files_after = self._list_matrix_files(matrix_dir)
        print(f"[generate_matrices] 生成后文件数: {len(files_after)}")
        
        new_ids = self._extract_new_ids(matrix_dir, files_before, files_after)
        print(f"[generate_matrices] 提取的新矩阵 IDs: {new_ids}")
        
        if new_ids and not self.dry_run:
            self.state.ensure_batch_matrices(new_ids)
            print(f"[generate_matrices] 已将 {len(new_ids)} 个矩阵添加到状态数据库")
        
        return new_ids

    # ------------------------------------------------------------------ #
    # 仿真执行                                                           #
    # ------------------------------------------------------------------ #

    def _awgn_log_paths(self, matrix_id: str, snr: float) -> Tuple[Path, Path]:
        log_dir = self.paths.sim_out_root / "awgn" / f"matrix{matrix_id}"
        self._ensure_dir(log_dir)
        # snr_str = f"{snr}".replace(".", "_")
        return log_dir / f"snr{snr}.log", log_dir / f"snr{snr}.err"

    def _errinj_log_paths(self, matrix_id: str) -> Tuple[Path, Path]:
        log_dir = self.paths.sim_out_root / "err_inj" / f"matrix{matrix_id}"
        self._ensure_dir(log_dir)
        return log_dir / "err_inj.log", log_dir / "err_inj.err"

    def _run_and_capture(
        self,
        cmd: Sequence[str],
        cwd: Optional[str],
        env: Optional[dict],
        log_path: Path,
        err_path: Path,
    ) -> None:
        result = self.runner.run(cmd, capture_output=True, cwd=cwd, env=env, check=True)
        if self.dry_run:
            return
        log_path.write_text(result.stdout or "", encoding="utf-8")
        err_path.write_text(result.stderr or "", encoding="utf-8")

    def _run_awgn_for_snrs(self, matrix_ids: Sequence[str], snrs: Sequence[float]) -> None:
        print(f"[_run_awgn_for_snrs] 矩阵: {matrix_ids}, SNR: {snrs}")
        if not matrix_ids or not snrs:
            print("[_run_awgn_for_snrs] 矩阵或SNR列表为空，跳过")
            return
        
        env = os.environ.copy()
        env["LDPC_MATRIX_DIR"] = str(self.sim_matrix_root)
        exec_dir = str(self.paths.exec_path.parent)
        print(f"[_run_awgn_for_snrs] 执行目录: {exec_dir}")
        print(f"[_run_awgn_for_snrs] LDPC_MATRIX_DIR: {self.sim_matrix_root}")

        for matrix_id in matrix_ids:
            for snr in snrs:
                log_path, err_path = self._awgn_log_paths(matrix_id, snr)
                cmd = [
                    str(self.paths.exec_path),
                    "LDPC",
                    str(self.paths.config_path),
                    "AWGN",
                    f"{snr}",
                    str(matrix_id),
                    str(self.sim_matrix_root),
                ]
                print(f"[_run_awgn_for_snrs] 运行: matrix_id={matrix_id}, SNR={snr}")
                print(f"[_run_awgn_for_snrs] 命令: {' '.join(cmd)}")
                print(f"[_run_awgn_for_snrs] 日志: {log_path}")
                
                self._run_and_capture(cmd, exec_dir, env, log_path, err_path)
                print(f"[_run_awgn_for_snrs] 完成: matrix_id={matrix_id}, SNR={snr}")

    def run_awgn(self, matrix_ids: Sequence[str]) -> None:
        self._run_awgn_for_snrs(matrix_ids, self.snr_list)

    def run_deep_awgn(self, matrix_ids: Sequence[str]) -> None:
        self._run_awgn_for_snrs(matrix_ids, self.deep_snr_list)

    def run_errinj(self, matrix_ids: Sequence[str]) -> None:
        if not matrix_ids:
            return
        env = os.environ.copy()
        env["LDPC_MATRIX_DIR"] = str(self.sim_matrix_root)
        exec_dir = str(self.paths.exec_path.parent)
        template = self.errinj_cfg.get("command")
        if not template:
            raise RuntimeError("err_inj.command 未配置。")
        
        # 获取 ERR_INJ 专用的配置文件路径
        errinj_config_value = self.errinj_cfg.get("config")
        if errinj_config_value:
            errinj_config_path = self.bundle.expand_path(errinj_config_value, base=self.paths.workspace_root)
            if errinj_config_path is None:
                # 如果扩展失败，尝试相对于 workspace_root 解析
                errinj_config_path = self.paths.workspace_root / errinj_config_value
        else:
            # 如果没有指定，使用默认的 config_path
            errinj_config_path = self.paths.config_path
        
        print(f"[run_errinj] ERR_INJ 配置文件: {errinj_config_path}")

        for matrix_id in matrix_ids:
            log_path, err_path = self._errinj_log_paths(matrix_id)
            formatted = template.format(
                config=str(errinj_config_path),
                matrix_id=str(matrix_id),
                matrix_dir=str(self.sim_matrix_root),
            )
            cmd = shlex.split(formatted)
            print(f"[run_errinj] 运行命令: {' '.join(cmd)}")
            self._run_and_capture(cmd, exec_dir, env, log_path, err_path)

    # ------------------------------------------------------------------ #
    # 日志刷新 / 绘图                                                    #
    # ------------------------------------------------------------------ #

    def refresh_awgn_logs(self) -> None:
        awgn_root = self.paths.sim_out_root / "awgn"
        results = scan_logs(
            log_root=awgn_root,
            model="AWGN",
            stat_keys=self.awgn_stat_keys,
            fer_key=self.awgn_fer_key,
        )
        for result in results:
            if result.fer_value is None:
                continue
            self.state.upsert_matrix(result.matrix_id)
            self.state.record_awgn_run(
                matrix_id=result.matrix_id,
                fer=result.fer_value,
                metrics=result.metrics,
                snr=result.snr,
                log_path=result.log_path,
                is_anchor=math.isclose(result.snr, self.anchor_snr, rel_tol=1e-6, abs_tol=1e-6),
            )

    def refresh_errinj_logs(self) -> None:
        err_root = self.paths.sim_out_root / "err_inj"
        results = scan_logs(
            log_root=err_root,
            model="ERR_INJ",
            stat_keys=self.errinj_stat_keys,
            fer_key=self.errinj_fer_key,
            snr_override=0.0,
        )
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

    def export_xlsx_if_enabled(self) -> None:
        xlsx_cfg = self.cfg["xlsx"]
        if not xlsx_cfg.get("enable", True):
            return
        filename = xlsx_cfg.get("filename", "simulation_results.xlsx")
        output_path = self.paths.perf_dir / filename
        records = self.state.fetch_runs_for_xlsx()
        if not records:
            return
        export_records(records, output_path)

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

    # ------------------------------------------------------------------ #
    # 主循环                                                             #
    # ------------------------------------------------------------------ #

    def should_stop(self) -> bool:
        if self.stop_flag_path and self.stop_flag_path.exists():
            print(f"检测到停止标志文件 {self.stop_flag_path} ，自动化即将结束。")
            return True
        if self.max_matrices and self.state.count_matrices() >= self.max_matrices:
            print("已达到配置的最大矩阵数量，停止提交新任务。")
            return True
        return False

    def run_iteration(self) -> None:
        print("\n" + "="*60)
        print("[DEBUG] 开始新的迭代")
        print("="*60)
        
        print("[DEBUG] 刷新 AWGN 日志...")
        self.refresh_awgn_logs()
        
        print("[DEBUG] 刷新 ERR_INJ 日志...")
        self.refresh_errinj_logs()

        print("[DEBUG] 获取 Top-N 矩阵...")
        top_records = self.state.fetch_top_matrices(self.top_n)
        print(f"[DEBUG] Top-{self.top_n} 矩阵数量: {len(top_records)}")
        
        pending_errinj = [rec.matrix_id for rec in top_records if rec.err_inj_status in ("pending", "failed")]
        if pending_errinj:
            print(f"[DEBUG] 运行 ERR_INJ: {pending_errinj}")
            self.run_errinj(pending_errinj)
        else:
            print("[DEBUG] 无需运行 ERR_INJ")

        ready_for_deep = [rec.matrix_id for rec in top_records if rec.passes_threshold and not rec.deep_runs_submitted]
        if ready_for_deep:
            print(f"[DEBUG] 运行深挖 AWGN: {ready_for_deep}")
            self.run_deep_awgn(ready_for_deep)
            for matrix_id in ready_for_deep:
                self.state.mark_deep_runs_submitted(matrix_id)
        else:
            print("[DEBUG] 无需运行深挖")

        batch_size = int(self.gen_cfg.get("batch_size") or 0)
        current_count = self.state.count_matrices()
        print(f"[DEBUG] 当前矩阵数: {current_count}, batch_size: {batch_size}, max_matrices: {self.max_matrices}")
        print(f"[DEBUG] should_stop(): {self.should_stop()}")
        
        if not self.should_stop() and batch_size > 0:
            print(f"[DEBUG] 生成 {batch_size} 个新矩阵...")
            new_ids = self.generate_matrices(batch_size)
            print(f"[DEBUG] 生成的矩阵 IDs: {new_ids}")
            if new_ids:
                print(f"[DEBUG] 运行 AWGN 仿真，矩阵: {new_ids}, SNR列表: {self.snr_list}")
                self.run_awgn(new_ids)
                print("[DEBUG] AWGN 仿真完成")
            else:
                print("[DEBUG] 警告：未生成新矩阵，跳过 AWGN")
        else:
            print(f"[DEBUG] 跳过矩阵生成 (should_stop={self.should_stop()}, batch_size={batch_size})")

        matrix_ids_for_plot = [rec.matrix_id for rec in top_records if rec.anchor_fer is not None]
        print(f"[DEBUG] 可绘图矩阵数: {len(matrix_ids_for_plot)}")
        self.export_xlsx_if_enabled()
        self.plot_top_if_enabled(matrix_ids_for_plot)
        print("[DEBUG] 迭代完成\n")

    def run(self, loop: bool = True) -> None:
        try:
            while True:
                self.run_iteration()
                if not loop or self.should_stop():
                    break
                time.sleep(self.sleep_seconds)
        finally:
            self.state.close()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LDPC 自动化（本地模式）")
    parser.add_argument("--config", required=True, help="配置文件路径 (JSON/YAML)")
    parser.add_argument("--dry-run", action="store_true", help="仅打印命令，不执行。")
    parser.add_argument("--once", action="store_true", help="仅运行一轮循环。")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    config_path = Path(args.config).resolve()
    bundle = load_bundle(config_path)
    runner = LocalAutomationRunner(bundle, dry_run=args.dry_run)
    runner.run(loop=not args.once)
    return 0


if __name__ == "__main__":
    sys.exit(main())

