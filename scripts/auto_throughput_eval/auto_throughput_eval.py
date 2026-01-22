#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

try:
    import tomllib  # Python 3.11+
except Exception:  # pragma: no cover
    tomllib = None  # type: ignore


def _ensure_repo_on_syspath() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    p = str(repo_root)
    if p not in sys.path:
        sys.path.insert(0, p)


_ensure_repo_on_syspath()

from auto_fer_eval.executors import LocalExecutor, LsfExecutor
from auto_fer_eval.log_parser import parse_log_text
from auto_fer_eval.manifest import RunRecord, find_run, load_manifest, save_manifest, upsert_run


@dataclass(frozen=True)
class CaseConfig:
    name: str
    workdir: str
    exe: str
    sim_mode: str
    config: str
    ch_model: str
    cmd_extra_args: list[str]
    out_dir: str
    log_prefix: str
    config_tag: str
    job_name_prefix: str

    snr_low: float
    snr_high: float
    snr_step: float

    max_sim_num: int
    max_err_num: int


@dataclass(frozen=True)
class RunnerConfig:
    executor: str = "local"  # local | lsf
    max_in_flight: int = 4
    poll_sec: float = 2.0
    fail_fast: bool = True


@dataclass(frozen=True)
class LsfConfig:
    queue: str = ""
    log_base_dir: str = ""
    bsub_extra: list[str] = field(default_factory=list)
    bjobs_extra: list[str] = field(default_factory=list)
    bkill_extra: list[str] = field(default_factory=list)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(prog="auto_throughput_eval", description="Fixed-grid throughput (aver_iter) runner.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run throughput grid for all cases")
    run.add_argument("--config", required=True, help="TOML config path")
    run.add_argument("--dry-run", action="store_true", help="Only print the plan; do not execute commands")

    exp = sub.add_parser("export", help="Parse existing logs/manifest and export csv")
    exp.add_argument("--config", required=True, help="TOML config path")

    ns = ap.parse_args(argv)
    if ns.cmd == "run":
        cases, runner, lsf = load_config(Path(ns.config))
        for c in cases:
            if runner.executor == "local":
                run_case_local(c, runner=runner, dry_run=bool(ns.dry_run))
            else:
                run_case_lsf(c, runner=runner, lsf=lsf, dry_run=bool(ns.dry_run))
            if not bool(ns.dry_run):
                export_case_csv(c, lsf=lsf)
        return 0
    if ns.cmd == "export":
        cases, _runner, lsf = load_config(Path(ns.config))
        for c in cases:
            export_case_csv(c, lsf=lsf)
        return 0

    raise SystemExit(f"Unknown cmd: {ns.cmd}")


def load_config(path: Path) -> tuple[list[CaseConfig], RunnerConfig, LsfConfig]:
    if tomllib is None:
        raise RuntimeError("tomllib not available in this Python; use Python>=3.11")
    data = tomllib.loads(path.read_text(encoding="utf-8"))

    defaults = data.get("defaults", {})
    runner_raw = data.get("runner", {})
    lsf_raw = data.get("lsf", {})
    case_list = data.get("cases", [])
    if not isinstance(case_list, list) or not case_list:
        raise ValueError("Config must contain [[cases]] list.")

    cases: list[CaseConfig] = []
    for raw in case_list:
        merged = dict(defaults)
        merged.update(raw or {})
        cases.append(_parse_case(merged))

    runner = _parse_runner(runner_raw)
    lsf = _parse_lsf(lsf_raw)
    return cases, runner, lsf


def _parse_runner(raw: Any) -> RunnerConfig:
    if not raw:
        return RunnerConfig()
    if not isinstance(raw, dict):
        raise ValueError("Config [runner] must be a table.")
    executor = str(raw.get("executor", "local")).strip().lower()
    if executor not in {"local", "lsf"}:
        raise ValueError("runner.executor must be 'local' or 'lsf'.")
    return RunnerConfig(
        executor=executor,
        max_in_flight=int(raw.get("max_in_flight", 4)),
        poll_sec=float(raw.get("poll_sec", 2.0)),
        fail_fast=bool(raw.get("fail_fast", True)),
    )


def _parse_lsf(raw: Any) -> LsfConfig:
    if not raw:
        return LsfConfig(queue="", log_base_dir="", bsub_extra=[], bjobs_extra=[], bkill_extra=[])
    if not isinstance(raw, dict):
        raise ValueError("Config [lsf] must be a table.")
    queue = str(raw.get("queue", "")).strip()
    log_base_dir = str(raw.get("log_base_dir", "")).strip()
    bsub_extra = [str(x) for x in raw.get("bsub_extra", [])]
    bjobs_extra = [str(x) for x in raw.get("bjobs_extra", [])]
    bkill_extra = [str(x) for x in raw.get("bkill_extra", [])]
    return LsfConfig(queue=queue, log_base_dir=log_base_dir, bsub_extra=bsub_extra, bjobs_extra=bjobs_extra, bkill_extra=bkill_extra)


def _parse_case(d: dict[str, Any]) -> CaseConfig:
    name = str(d["name"]).strip()
    if not name:
        raise ValueError("case.name must be non-empty.")

    workdir = str(d.get("workdir", "."))
    exe = str(d["exe"])
    sim_mode = str(d.get("sim_mode", "LDPC"))
    config = str(d["config"])
    ch_model = str(d["ch_model"])
    cmd_extra_args_raw = d.get("cmd_extra_args", [])
    if cmd_extra_args_raw is None:
        cmd_extra_args_raw = []
    if isinstance(cmd_extra_args_raw, str) or not isinstance(cmd_extra_args_raw, list):
        raise ValueError('cmd_extra_args must be a TOML array, e.g. cmd_extra_args=["0","./matrix"].')
    cmd_extra_args = [str(x) for x in cmd_extra_args_raw]

    out_dir = str(d.get("out_dir", "./perf_auto_throughput"))
    log_prefix = str(d.get("log_prefix", "tp"))
    config_tag = str(d.get("config_tag", "")).strip()
    if not config_tag:
        raise ValueError("config_tag must be provided (you define it in TOML).")
    job_name_prefix = str(d.get("job_name_prefix", "Throughput Data")).strip()

    snr_low = float(d["snr_low"])
    snr_high = float(d["snr_high"])
    snr_step = float(d["snr_step"])
    if not math.isfinite(snr_low) or not math.isfinite(snr_high) or not math.isfinite(snr_step) or snr_step <= 0:
        raise ValueError("snr_low/snr_high/snr_step must be finite, with snr_step>0.")
    if snr_high < snr_low - 1e-12:
        raise ValueError("Require snr_high >= snr_low.")

    max_sim_num = int(d.get("max_sim_num", 0))
    max_err_num = int(d.get("max_err_num", 0))
    if max_sim_num <= 0:
        raise ValueError("max_sim_num must be > 0.")
    if max_err_num < 0:
        raise ValueError("max_err_num must be >= 0.")

    return CaseConfig(
        name=name,
        workdir=workdir,
        exe=exe,
        sim_mode=sim_mode,
        config=config,
        ch_model=ch_model,
        cmd_extra_args=cmd_extra_args,
        out_dir=out_dir,
        log_prefix=log_prefix,
        config_tag=config_tag,
        job_name_prefix=job_name_prefix,
        snr_low=snr_low,
        snr_high=snr_high,
        snr_step=snr_step,
        max_sim_num=max_sim_num,
        max_err_num=max_err_num,
    )


def run_case_local(c: CaseConfig, *, runner: RunnerConfig, dry_run: bool) -> None:
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    snrs = list(build_snr_grid(c.snr_low, c.snr_high, c.snr_step))
    todo: list[tuple[float, Path, list[str]]] = []
    cached = 0
    for snr in snrs:
        log_path = log_path_for(c, snr)
        cmd = build_cmd(c, snr)

        found = find_run(manifest, "snr", snr)
        if found and found.get("raw_ber") is not None and found.get("retry_dec_avg_iter") is not None:
            cached += 1
            continue

        rec_done = parse_done_log(snr=snr, log_path=log_path, cmd=cmd)
        if rec_done is not None:
            upsert_run(manifest, rec_done)
            cached += 1
            continue
        todo.append((snr, log_path, cmd))

    save_manifest(manifest_path, manifest)

    print(f"\n[auto_throughput_eval] case={c.name} grid={len(snrs)} cached={cached} todo={len(todo)} out={case_dir}")
    if dry_run:
        if todo:
            print(
                f"[dry-run] would patch config in-place: {c.config} "
                f"(max_sim_num={c.max_sim_num}, max_err_num={c.max_err_num})"
            )
        for snr, log_path, cmd in todo:
            print(f"[dry-run] snr={snr} log={log_path}")
            print(f"          cmd: {' '.join(cmd)}")
        return

    cfg_abs = resolve_cfg_path(c)
    orig_cfg_text: Optional[str] = None
    if todo:
        orig_cfg_text = cfg_abs.read_text(encoding="utf-8")
        patch_config_max_sim_num(cfg_abs, cfg_abs, c.max_sim_num, max_error_num=c.max_err_num)

    ex = LocalExecutor()
    max_in_flight = max(1, int(runner.max_in_flight))
    poll_sec = max(0.05, float(runner.poll_sec))

    in_flight: dict[float, dict[str, Any]] = {}

    def cancel_all(reason: str) -> None:
        for xv, info in list(in_flight.items()):
            try:
                ex.cancel(info["job"])
            except Exception:
                pass
            print(f"[auto_throughput_eval] cancel snr={xv} reason={reason}")
            in_flight.pop(xv, None)

    try:
        while todo or in_flight:
            while todo and len(in_flight) < max_in_flight:
                snr, log_path, cmd = todo.pop(0)
                job = ex.submit(cmd, cwd=c.workdir, log_path=log_path)
                in_flight[snr] = {"job": job, "log_path": log_path, "cmd": cmd}
                print(f"[auto_throughput_eval] submit snr={snr} -> {log_path.name}")

            if not in_flight:
                break

            time.sleep(poll_sec)
            for snr, info in list(in_flight.items()):
                st = ex.poll(info["job"])
                if not st.done:
                    continue

                in_flight.pop(snr, None)
                if not st.ok:
                    msg = f"job EXIT: snr={snr} log={info['log_path']}"
                    print(f"[auto_throughput_eval] {msg}")
                    if runner.fail_fast:
                        cancel_all(reason=msg)
                        save_manifest(manifest_path, manifest)
                        raise SystemExit(msg)
                    continue

                rec = parse_done_log(snr=snr, log_path=Path(info["log_path"]), cmd=list(info["cmd"]))
                if rec is None:
                    msg = f"parse failed: snr={snr} log={info['log_path']}"
                    print(f"[auto_throughput_eval] {msg}")
                    if runner.fail_fast:
                        cancel_all(reason=msg)
                        save_manifest(manifest_path, manifest)
                        raise SystemExit(msg)
                    continue

                upsert_run(manifest, rec)
                save_manifest(manifest_path, manifest)
                print(f"[auto_throughput_eval] done snr={snr} RBER={rec.raw_ber} aver_iter={rec.retry_dec_avg_iter}")
    finally:
        if orig_cfg_text is not None:
            cfg_abs.write_text(orig_cfg_text, encoding="utf-8")


def run_case_lsf(c: CaseConfig, *, runner: RunnerConfig, lsf: LsfConfig, dry_run: bool) -> None:
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    cfg_abs = resolve_cfg_path(c)
    matrix_size = parse_matrix_size(cfg_abs)
    snrs = list(build_snr_grid(c.snr_low, c.snr_high, c.snr_step))

    todo: list[tuple[float, Path, list[str], str]] = []
    cached = 0
    for snr in snrs:
        log_path = log_path_for_submission(c, snr, lsf=lsf)
        cmd = build_cmd(c, snr)
        job_name = build_job_name(c, matrix_size=matrix_size, snr=snr)

        found = find_run(manifest, "snr", snr)
        if found and found.get("raw_ber") is not None and found.get("retry_dec_avg_iter") is not None:
            cached += 1
            continue

        rec_done = parse_done_log(snr=snr, log_path=log_path, cmd=cmd)
        if rec_done is not None:
            upsert_run(manifest, rec_done)
            cached += 1
            continue

        todo.append((snr, log_path, cmd, job_name))

    save_manifest(manifest_path, manifest)

    print(f"\n[auto_throughput_eval] case={c.name} grid={len(snrs)} cached={cached} todo={len(todo)} out={case_dir}")
    if dry_run:
        if todo:
            print(
                f"[dry-run] would patch config in-place: {c.config} "
                f"(max_sim_num={c.max_sim_num}, max_err_num={c.max_err_num})"
            )
        cwd_abs = str(Path(c.workdir).resolve())
        for snr, log_path, cmd, job_name in todo:
            bsub_cmd = build_bsub_command_line(
                cmd,
                cwd_abs=cwd_abs,
                log_path=str(log_path),
                job_name=job_name,
                queue=lsf.queue,
                bsub_extra=lsf.bsub_extra,
            )
            print(f"[dry-run] snr={snr}")
            print(f"          {bsub_cmd}")
        return

    orig_cfg_text: Optional[str] = None
    if todo:
        orig_cfg_text = cfg_abs.read_text(encoding="utf-8")
        patch_config_max_sim_num(cfg_abs, cfg_abs, c.max_sim_num, max_error_num=c.max_err_num)

    ex = LsfExecutor(
        queue=lsf.queue,
        bsub_extra=lsf.bsub_extra,
        bjobs_extra=lsf.bjobs_extra,
        bkill_extra=lsf.bkill_extra,
        log_base_dir="",  # log path mapping is handled by log_path_for_submission()
    )
    max_in_flight = max(1, int(runner.max_in_flight))
    poll_sec = max(0.2, float(runner.poll_sec))

    in_flight: dict[float, dict[str, Any]] = {}

    def cancel_all(reason: str) -> None:
        for xv, info in list(in_flight.items()):
            try:
                ex.cancel(info["job"])
            except Exception:
                pass
            print(f"[auto_throughput_eval] cancel snr={xv} reason={reason}")
            in_flight.pop(xv, None)
    try:
        while todo or in_flight:
            while todo and len(in_flight) < max_in_flight:
                snr, log_path, cmd, job_name = todo.pop(0)
                job = ex.submit(cmd, cwd=c.workdir, log_path=Path(log_path), job_name=job_name, queue=lsf.queue)
                in_flight[snr] = {"job": job, "log_path": log_path, "cmd": cmd, "job_name": job_name}
                jid = getattr(job, "job_id", "?")
                print(f"[auto_throughput_eval] bsub snr={snr} job_id={jid} job_name={job_name}")

            if not in_flight:
                break

            time.sleep(poll_sec)
            for snr, info in list(in_flight.items()):
                st = ex.poll(info["job"])
                if not st.done:
                    continue

                in_flight.pop(snr, None)
                if not st.ok:
                    msg = f"job EXIT: snr={snr} log={info['log_path']} job_name={info.get('job_name')}"
                    print(f"[auto_throughput_eval] {msg}")
                    if runner.fail_fast:
                        cancel_all(reason=msg)
                        save_manifest(manifest_path, manifest)
                        raise SystemExit(msg)
                    continue

                rec = parse_done_log(snr=snr, log_path=Path(info["log_path"]), cmd=list(info["cmd"]))
                if rec is None:
                    msg = f"parse failed: snr={snr} log={info['log_path']} job_name={info.get('job_name')}"
                    print(f"[auto_throughput_eval] {msg}")
                    if runner.fail_fast:
                        cancel_all(reason=msg)
                        save_manifest(manifest_path, manifest)
                        raise SystemExit(msg)
                    continue

                upsert_run(manifest, rec)
                save_manifest(manifest_path, manifest)
                print(f"[auto_throughput_eval] done snr={snr} RBER={rec.raw_ber} aver_iter={rec.retry_dec_avg_iter}")
    finally:
        if orig_cfg_text is not None:
            cfg_abs.write_text(orig_cfg_text, encoding="utf-8")


def patch_config_max_sim_num(src: Path, dst: Path, max_sim_num: int, max_error_num: Optional[int] = None) -> None:
    txt = src.read_text(encoding="utf-8")
    lines = txt.splitlines()

    idx_sim = None
    for i, line in enumerate(lines):
        if "maximum simulation number" in line:
            idx_sim = i
            break
    if idx_sim is None:
        idx_sim = 8

    line = lines[idx_sim]
    if "//" in line:
        left, right = line.split("//", 1)
        comment = "//" + right
    else:
        left, comment = line, ""
    left = left.rstrip("\n")

    import re

    new_left = re.sub(r"^\s*\d+", str(int(max_sim_num)), left)
    lines[idx_sim] = (new_left.rstrip() + (" " + comment if comment else "")).rstrip()

    if max_error_num is not None:
        idx_err = None
        for i, line in enumerate(lines):
            if "maximum error number" in line:
                idx_err = i
                break
        if idx_err is None:
            idx_err = 9

        line = lines[idx_err]
        if "//" in line:
            left, right = line.split("//", 1)
            comment = "//" + right
        else:
            left, comment = line, ""
        left = left.rstrip("\n")

        new_left = re.sub(r"^\s*\d+", str(int(max_error_num)), left)
        lines[idx_err] = (new_left.rstrip() + (" " + comment if comment else "")).rstrip()

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8")


def resolve_cfg_path(c: CaseConfig) -> Path:
    p = Path(c.config)
    if p.is_absolute():
        return p
    return (Path(c.workdir) / p).resolve()


def build_cmd(c: CaseConfig, snr: float) -> list[str]:
    # Keep the config path as-is from TOML (typically "config/xxx.cnfg"),
    # so that the bsub command does not depend on out_dir/log_base_dir.
    cmd: list[str] = [c.exe, c.sim_mode, c.config, c.ch_model]
    if c.ch_model.upper() != "CLEAN":
        cmd.append(str(snr))
    cmd.extend(c.cmd_extra_args)
    return cmd


def parse_done_log(*, snr: float, log_path: Path, cmd: list[str]) -> Optional[RunRecord]:
    if not log_path.exists():
        return None
    txt = log_path.read_text(encoding="utf-8", errors="replace")
    m = parse_log_text(txt)
    if m.raw_ber is None or m.retry_dec_avg_iter is None:
        return None
    return RunRecord(
        axis_type="snr",
        axis_value=float(snr),
        log_path=str(log_path),
        cmd=cmd,
        raw_ber=m.raw_ber,
        theo_rber=m.theo_rber,
        ldpc_fer=m.ldpc_fer,
        fail_cw=m.fail_cw,
        total_packets=m.total_packets,
        retry_dec_avg_iter=m.retry_dec_avg_iter,
    )


def log_path_for(c: CaseConfig, snr: float) -> Path:
    case_dir = Path(c.out_dir) / c.name
    return (case_dir / f"{c.log_prefix}_{format_snr(snr)}.log").resolve()


def format_snr(x: float) -> str:
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return f"snr{s}"


def build_snr_grid(low: float, high: float, step: float) -> list[float]:
    origin = float(low)
    x = quantize_axis(origin, step, origin=origin)
    out: list[float] = []
    while x <= high + 1e-12:
        out.append(float(x))
        x = quantize_axis(x + step, step, origin=origin)
        if len(out) > 1000000:
            raise RuntimeError("SNR grid too large; check snr_step.")
    return out


def step_decimals(step: float) -> int:
    if not math.isfinite(step) or step <= 0:
        return 12
    s = f"{step:.12g}"
    if "e" in s or "E" in s:
        return 12
    if "." not in s:
        return 0
    return len(s.split(".")[1].rstrip("0"))


def quantize_axis(x: float, step: float, *, origin: float = 0.0) -> float:
    if step <= 0:
        return x
    q = origin + round((x - origin) / step) * step
    q = round(q, step_decimals(step))
    if abs(q) < 1e-12:
        q = 0.0
    return float(q)


def log_path_for_submission(c: CaseConfig, snr: float, *, lsf: LsfConfig) -> Path:
    """
    Determine the log path that LSF should write to.

    If lsf.log_base_dir is set, redirect logs to:
      <log_base_dir>/<case_name>/<log_file>
    Otherwise use:
      <out_dir>/<case_name>/<log_file>
    """
    log_path = log_path_for(c, snr)
    if not lsf.log_base_dir:
        return log_path

    case_dir = (Path(c.out_dir) / c.name).resolve()
    try:
        rel = log_path.resolve().relative_to(case_dir)
    except Exception:
        rel = Path(log_path.name)
    return (Path(lsf.log_base_dir) / c.name / rel).resolve()


def parse_matrix_size(cfg_path: Path) -> str:
    """
    Parse base-matrix size from a *.cnfg file.
    Prefer matching marker comments; fallback to line 4/5 (1-based).
    """
    txt = cfg_path.read_text(encoding="utf-8", errors="replace")
    lines = txt.splitlines()

    def _first_int(s: str) -> Optional[int]:
        import re

        m = re.search(r"^\s*([0-9]+)\b", s)
        if not m:
            return None
        try:
            return int(m.group(1))
        except Exception:
            return None

    row = None
    col = None
    for ln in lines:
        if row is None and "row number of base matrix" in ln:
            row = _first_int(ln)
        if col is None and "column number of base matrix" in ln:
            col = _first_int(ln)
        if row is not None and col is not None:
            break

    if row is None and len(lines) >= 4:
        row = _first_int(lines[3])
    if col is None and len(lines) >= 5:
        col = _first_int(lines[4])

    if row is None or col is None:
        raise ValueError(f"Cannot parse matrix size from cnfg: {cfg_path}")
    return f"{row}x{col}"


def build_job_name(c: CaseConfig, *, matrix_size: str, snr: float) -> str:
    snr_s = format_snr_value(snr)
    prefix = c.job_name_prefix.strip() or "Throughput Data"
    return f"{prefix} {matrix_size} {c.config_tag} snr-{snr_s}"


def format_snr_value(x: float) -> str:
    return f"{x:.3f}".rstrip("0").rstrip(".")


def build_bsub_command_line(
    cmd: list[str],
    *,
    cwd_abs: str,
    log_path: str,
    job_name: str,
    queue: str,
    bsub_extra: list[str],
) -> str:
    parts: list[str] = ["bsub"]
    if queue:
        parts.extend(["-q", quote_shell(queue)])
    if job_name:
        parts.extend(["-J", quote_shell(job_name)])
    parts.extend(["-o", quote_shell(log_path), "-cwd", quote_shell(cwd_abs)])
    parts.extend([quote_shell(x) for x in bsub_extra])
    # The simulation command is shown as a single field (one pair of double quotes)
    # for readability and to match the user's preferred bsub CLI style.
    cmdline = " ".join([str(x) for x in cmd])
    parts.append(quote_shell(cmdline))
    return " ".join(parts)


def quote_shell(s: str) -> str:
    # Display helper for dry-run printing (not used for subprocess execution).
    s2 = str(s)
    # Use double quotes for all fields (for consistent bsub command display).
    # Escape chars that are special in double quotes or interactive shells.
    s2 = (
        s2.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("$", "\\$")
        .replace("`", "\\`")
        .replace("!", "\\!")
    )
    return f'"{s2}"'


def export_case_csv(c: CaseConfig, *, lsf: LsfConfig) -> Path:
    """
    Export snr/RBER/aver_iter to a CSV under the case directory.
    The data source is manifest.json (filled from parsing logs if needed).
    """
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    snrs = list(build_snr_grid(c.snr_low, c.snr_high, c.snr_step))
    rows: list[list[object]] = []
    missing: list[float] = []

    for snr in snrs:
        rec = find_run(manifest, "snr", snr)
        raw_ber = rec.get("raw_ber") if rec else None
        avg_it = rec.get("retry_dec_avg_iter") if rec else None

        if raw_ber is None or avg_it is None:
            log_p = None
            if rec and rec.get("log_path"):
                log_p = Path(str(rec["log_path"]))
            if log_p is None:
                log_p = log_path_for_submission(c, snr, lsf=lsf)
            cmd = list(rec.get("cmd")) if rec and isinstance(rec.get("cmd"), list) else build_cmd(c, snr)
            rec_done = parse_done_log(snr=snr, log_path=log_p, cmd=cmd)
            if rec_done is not None:
                upsert_run(manifest, rec_done)
                raw_ber = rec_done.raw_ber
                avg_it = rec_done.retry_dec_avg_iter

        if raw_ber is None or avg_it is None:
            missing.append(float(snr))
        rows.append([float(snr), None if raw_ber is None else float(raw_ber), None if avg_it is None else float(avg_it)])

    save_manifest(manifest_path, manifest)

    csv_path = (case_dir / "throughput.csv").resolve()
    write_csv(csv_path, headers=["snr", "RBER", "aver_iter"], rows=rows)
    print(f"[auto_throughput_eval] export csv: {csv_path}")

    if missing:
        miss_s = ", ".join([format_snr_value(x) for x in missing[:20]])
        more = "" if len(missing) <= 20 else f" ...(+{len(missing)-20})"
        raise SystemExit(f"Missing metrics for SNR points: {miss_s}{more}. See {manifest_path}")

    return csv_path


def write_csv(path: Path, *, headers: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(list(headers))
        for r in rows:
            out: list[str] = []
            for v in r:
                if v is None:
                    out.append("")
                elif isinstance(v, (int, float, bool)):
                    out.append(format_number(v))
                else:
                    out.append(str(v))
            w.writerow(out)


def format_number(v: object) -> str:
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # Keep enough precision for typical BER values while avoiding long tails.
        return f"{v:.15g}"
    return str(v)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
