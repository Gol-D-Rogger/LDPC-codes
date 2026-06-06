#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
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
from auto_fer_eval.job_db import JobDB
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

    axis_type: str  # "snr" or "k"
    x_low: float
    x_high: float
    x_step: float

    max_sim_num: int
    max_err_num: int


@dataclass(frozen=True)
class RunnerConfig:
    executor: str = "local"  # local | lsf
    max_in_flight: int = 4
    poll_sec: float = 2.0
    fail_fast: bool = False
    # Log flush grace window (seconds) after LSF/local reports DONE/EXIT.
    # Shared filesystems (NFS) can delay creation/flush of bsub -o logs.
    timeout_log_grace_sec: float = 10.0
    max_retries: int = 3


@dataclass(frozen=True)
class LsfConfig:
    queue: str = ""
    use_cwd: bool = True
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
    # Support both [runner] and [throughput] for backward compatibility
    runner_raw = data.get("runner", data.get("throughput", {}))
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

    fail_fast_raw = raw.get("fail_fast", False)
    # TOML should provide a proper boolean (true/false). Be defensive:
    # - allow 0/1
    # - allow strings like "false"/"true" (but recommend using boolean)
    fail_fast: bool
    if isinstance(fail_fast_raw, bool):
        fail_fast = fail_fast_raw
    elif isinstance(fail_fast_raw, (int, float)) and fail_fast_raw in (0, 1):
        fail_fast = bool(fail_fast_raw)
    elif isinstance(fail_fast_raw, str):
        s = fail_fast_raw.strip().lower()
        if s in {"false", "0", "no", "n"}:
            fail_fast = False
        elif s in {"true", "1", "yes", "y"}:
            fail_fast = True
        else:
            raise ValueError("runner.fail_fast must be a boolean true/false (do not quote).")
    else:
        raise ValueError("runner.fail_fast must be a boolean true/false (do not quote).")

    return RunnerConfig(
        executor=executor,
        max_in_flight=int(raw.get("max_in_flight", 4)),
        poll_sec=float(raw.get("poll_sec", 2.0)),
        fail_fast=fail_fast,
        timeout_log_grace_sec=float(raw.get("timeout_log_grace_sec", 10.0)),
        max_retries=int(raw.get("max_retries", 3)),
    )


def _parse_lsf(raw: Any) -> LsfConfig:
    if not raw:
        return LsfConfig(queue="", use_cwd=True, log_base_dir="", bsub_extra=[], bjobs_extra=[], bkill_extra=[])
    if not isinstance(raw, dict):
        raise ValueError("Config [lsf] must be a table.")
    queue = str(raw.get("queue", "")).strip()
    use_cwd_raw = raw.get("use_cwd", True)
    if isinstance(use_cwd_raw, bool):
        use_cwd = use_cwd_raw
    elif isinstance(use_cwd_raw, (int, float)) and use_cwd_raw in (0, 1):
        use_cwd = bool(use_cwd_raw)
    elif isinstance(use_cwd_raw, str):
        s = use_cwd_raw.strip().lower()
        if s in {"true", "1", "yes", "y"}:
            use_cwd = True
        elif s in {"false", "0", "no", "n"}:
            use_cwd = False
        else:
            raise ValueError("lsf.use_cwd must be a boolean true/false (or 1/0).")
    else:
        raise ValueError("lsf.use_cwd must be a boolean true/false (or 1/0).")
    log_base_dir = str(raw.get("log_base_dir", "")).strip()
    bsub_extra = [str(x) for x in raw.get("bsub_extra", [])]
    bjobs_extra = [str(x) for x in raw.get("bjobs_extra", [])]
    bkill_extra = [str(x) for x in raw.get("bkill_extra", [])]
    return LsfConfig(
        queue=queue,
        use_cwd=use_cwd,
        log_base_dir=log_base_dir,
        bsub_extra=bsub_extra,
        bjobs_extra=bjobs_extra,
        bkill_extra=bkill_extra,
    )


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

    axis_type = str(d.get("axis_type", "")).strip().lower()
    if not axis_type:
        axis_type = "k" if str(ch_model).upper() == "ERR_INJ" else "snr"
    if axis_type not in {"snr", "k"}:
        raise ValueError(f"Unsupported axis_type: {axis_type}")

    # Support both snr_min/snr_max/step and k_min/k_max/step formats (like auto_fer_eval)
    if axis_type == "k":
        x_low = float(d.get("k_min", d.get("k_low", d.get("snr_min", d.get("snr_low")))))
        x_high = float(d.get("k_max", d.get("k_high", d.get("snr_max", d.get("snr_high")))))
        x_step = float(d.get("step", d.get("k_step", d.get("snr_step"))))
    else:
        x_low = float(d.get("snr_min", d.get("snr_low")))
        x_high = float(d.get("snr_max", d.get("snr_high")))
        x_step = float(d.get("step", d.get("snr_step")))

    if not math.isfinite(x_low) or not math.isfinite(x_high) or not math.isfinite(x_step) or x_step <= 0:
        raise ValueError("axis low/high/step must be finite, with step>0.")
    if x_high < x_low - 1e-12:
        raise ValueError("Require x_high >= x_low.")

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
        axis_type=axis_type,
        x_low=x_low,
        x_high=x_high,
        x_step=x_step,
        max_sim_num=max_sim_num,
        max_err_num=max_err_num,
    )


def run_case_local(c: CaseConfig, *, runner: RunnerConfig, dry_run: bool) -> None:
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    cfg_abs = resolve_cfg_path(c)
    runtime_cfg = (case_dir / "config_runtime.cnfg").resolve()

    xs = list(build_axis_grid(c.x_low, c.x_high, c.x_step, axis_type=c.axis_type))
    todo: list[tuple[float, Path, list[str]]] = []
    cached = 0
    for xv in xs:
        log_path = log_path_for(c, xv)
        cmd = build_cmd(c, xv, config_arg=str(runtime_cfg))

        rec_done = _resolve_complete_record(
            axis_type=c.axis_type,
            axis_value=xv,
            manifest_rec=find_run(manifest, c.axis_type, xv),
            fallback_log_path=log_path,
            fallback_cmd=cmd,
        )
        if rec_done is not None:
            upsert_run(manifest, rec_done)
            cached += 1
            continue
        todo.append((xv, log_path, cmd))

    save_manifest(manifest_path, manifest)

    print(f"\n[auto_throughput_eval] case={c.name} axis={c.axis_type} grid={len(xs)} cached={cached} todo={len(todo)} out={case_dir}")
    if dry_run:
        if todo:
            print(
                f"[dry-run] would generate runtime config: {runtime_cfg} (from {cfg_abs}) "
                f"(max_sim_num={c.max_sim_num}, max_err_num={c.max_err_num})"
            )
        for xv, log_path, cmd in todo:
            print(f"[dry-run] axis={c.axis_type} x={xv} log={log_path}")
            print(f"          cmd: {' '.join(cmd)}")
        return

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
            print(f"[auto_throughput_eval] cancel axis={c.axis_type} x={xv} reason={reason}")
            in_flight.pop(xv, None)

    try:
        if todo:
            patch_config_max_sim_num(cfg_abs, runtime_cfg, c.max_sim_num, max_error_num=c.max_err_num)

        submit_attempts: dict[float, int] = {}
        while todo or in_flight:
            while todo and len(in_flight) < max_in_flight:
                xv, log_path, cmd = todo.pop(0)
                if log_path.exists():
                    try:
                        log_path.unlink()
                    except Exception:
                        pass
                job = ex.submit(cmd, cwd=c.workdir, log_path=log_path)
                in_flight[xv] = {"job": job, "log_path": log_path, "cmd": cmd}
                submit_attempts[xv] = int(submit_attempts.get(xv, 0)) + 1
                print(
                    f"[auto_throughput_eval] submit axis={c.axis_type} x={xv} "
                    f"attempt={submit_attempts[xv]} -> {log_path.name}"
                )

            if not in_flight:
                break

            time.sleep(poll_sec)
            for xv, info in list(in_flight.items()):
                st = ex.poll(info["job"])
                if not st.done:
                    continue

                in_flight.pop(xv, None)

                rec = parse_done_log_with_grace(
                    axis_type=c.axis_type,
                    axis_value=xv,
                    log_path=Path(info["log_path"]),
                    cmd=list(info["cmd"]),
                    grace_sec=float(runner.timeout_log_grace_sec),
                )

                retry_reason: Optional[str] = None
                if not st.ok:
                    retry_reason = f"job EXIT: axis={c.axis_type} x={xv} log={info['log_path']}"
                elif rec is None:
                    retry_reason = f"log incomplete: axis={c.axis_type} x={xv} log={info['log_path']}"

                if retry_reason is not None:
                    print(f"[auto_throughput_eval] {retry_reason}")
                    if runner.fail_fast:
                        cancel_all(reason=retry_reason)
                        save_manifest(manifest_path, manifest)
                        raise SystemExit(retry_reason)
                    next_attempt = int(submit_attempts.get(xv, 0)) + 1
                    if next_attempt > runner.max_retries:
                        print(
                            f"[auto_throughput_eval] SKIP axis={c.axis_type} x={xv}: "
                            f"exceeded max_retries={runner.max_retries}"
                        )
                    else:
                        todo.append((float(xv), Path(info["log_path"]).resolve(), list(info["cmd"])))
                        print(
                            f"[auto_throughput_eval] retry axis={c.axis_type} x={xv} "
                            f"next_attempt={next_attempt}"
                        )
                    continue

                upsert_run(manifest, rec)
                save_manifest(manifest_path, manifest)
                print(f"[auto_throughput_eval] done axis={c.axis_type} x={xv} RBER={rec.raw_ber} aver_iter={rec.retry_dec_avg_iter}")
    finally:
        pass


def run_case_lsf(c: CaseConfig, *, runner: RunnerConfig, lsf: LsfConfig, dry_run: bool) -> None:
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)
    job_db_path = (case_dir / "jobs.json").resolve()
    jobs_by_point = _load_job_records_by_point(job_db_path, axis_type=c.axis_type)

    cfg_abs = resolve_cfg_path(c)
    runtime_cfg = (case_dir / "config_runtime.cnfg").resolve()
    matrix_size = parse_matrix_size(cfg_abs)
    xs = list(build_axis_grid(c.x_low, c.x_high, c.x_step, axis_type=c.axis_type))

    if dry_run:
        todo: list[tuple[float, Path, list[str], str]] = []
        cached = 0
        active_recorded = 0
        for xv in xs:
            log_path = log_path_for_submission(c, xv, lsf=lsf)
            cmd = build_cmd(c, xv, config_arg=str(runtime_cfg))
            job_name = build_job_name(c, matrix_size=matrix_size, axis_value=xv)
            point_key = _throughput_point_key(c.axis_type, xv)

            rec_done = _resolve_complete_record(
                axis_type=c.axis_type,
                axis_value=xv,
                manifest_rec=find_run(manifest, c.axis_type, xv),
                fallback_log_path=log_path,
                fallback_cmd=cmd,
            )
            if rec_done is not None:
                upsert_run(manifest, rec_done)
                cached += 1
                continue

            if any(_job_record_is_active(rec) for rec in jobs_by_point.get(point_key, [])):
                active_recorded += 1
                continue
            todo.append((xv, log_path, cmd, job_name))

        save_manifest(manifest_path, manifest)
        print(
            f"\n[auto_throughput_eval] case={c.name} axis={c.axis_type} grid={len(xs)} "
            f"cached={cached} active_recorded={active_recorded} todo={len(todo)} out={case_dir}"
        )
        if todo:
            print(
                f"[dry-run] would generate runtime config: {runtime_cfg} (from {cfg_abs}) "
                f"(max_sim_num={c.max_sim_num}, max_err_num={c.max_err_num})"
            )
        if active_recorded:
            print(f"[dry-run] recorded active jobs from jobs.json: {active_recorded}")
        cwd_abs = str(Path(c.workdir).resolve())
        for xv, log_path, cmd, job_name in todo:
            bsub_cmd = build_bsub_command_line(
                cmd,
                cwd_abs=cwd_abs,
                log_path=str(log_path),
                job_name=job_name,
                queue=lsf.queue,
                bsub_extra=lsf.bsub_extra,
                include_cwd=lsf.use_cwd,
            )
            print(f"[dry-run] axis={c.axis_type} x={xv}")
            print(f"          {bsub_cmd}")
        return

    ex = LsfExecutor(
        queue=lsf.queue,
        bsub_extra=lsf.bsub_extra,
        bjobs_extra=lsf.bjobs_extra,
        bkill_extra=lsf.bkill_extra,
        include_cwd=lsf.use_cwd,
        log_base_dir="",  # log path mapping is handled by log_path_for_submission()
    )
    job_db = JobDB(job_db_path)
    max_in_flight = max(1, int(runner.max_in_flight))
    poll_sec = max(0.2, float(runner.poll_sec))

    in_flight: dict[float, dict[str, Any]] = {}
    todo: list[tuple[float, Path, list[str], str, str]] = []
    cached = 0
    adopted = 0

    for xv in xs:
        log_path = log_path_for_submission(c, xv, lsf=lsf)
        cmd = build_cmd(c, xv, config_arg=str(runtime_cfg))
        job_name = build_job_name(c, matrix_size=matrix_size, axis_value=xv)
        point_key = _throughput_point_key(c.axis_type, xv)

        rec_done = _resolve_complete_record(
            axis_type=c.axis_type,
            axis_value=xv,
            manifest_rec=find_run(manifest, c.axis_type, xv),
            fallback_log_path=log_path,
            fallback_cmd=cmd,
        )
        if rec_done is not None:
            upsert_run(manifest, rec_done)
            cached += 1
            continue

        reused = False
        for job_rec in jobs_by_point.get(point_key, []):
            job_key = str(job_rec.get("key") or "")
            job_log_path = Path(job_rec.get("log_path") or log_path).resolve()
            job_cmd_raw = job_rec.get("cmd")
            job_cmd = [str(x) for x in job_cmd_raw] if isinstance(job_cmd_raw, list) else list(cmd)

            rec_done = parse_done_log(
                axis_type=c.axis_type,
                axis_value=xv,
                log_path=job_log_path,
                cmd=job_cmd,
                debug=False,
            )
            if rec_done is not None:
                upsert_run(manifest, rec_done)
                cached += 1
                now = time.time()
                if job_key:
                    job_db.set_state(job_key, "DONE", t=now)
                    job_db.upsert({"key": job_key, "finished_at": now})
                reused = True
                break

            if not _job_record_is_active(job_rec):
                continue

            job = _resume_executor_job(ex, cmd=job_cmd, cwd=c.workdir, log_path=job_log_path, rec=job_rec)
            if job is None:
                continue

            now = time.time()
            st = ex.poll(job)
            if job_key:
                job_db.set_state(job_key, st.state, t=now)

            if st.done:
                rec_done = parse_done_log_with_grace(
                    axis_type=c.axis_type,
                    axis_value=xv,
                    log_path=job_log_path,
                    cmd=job_cmd,
                    grace_sec=max(float(runner.timeout_log_grace_sec), 30.0),
                )
                if rec_done is not None:
                    upsert_run(manifest, rec_done)
                    cached += 1
                    if job_key:
                        job_db.set_state(job_key, "DONE", t=now)
                        job_db.upsert({"key": job_key, "finished_at": now})
                    reused = True
                    break
                if job_key:
                    job_db.upsert({"key": job_key, "finished_at": now})
                continue

            in_flight[xv] = {
                "job": job,
                "job_key": job_key,
                "log_path": job_log_path,
                "cmd": job_cmd,
                "job_name": str(job_rec.get("job_name") or job_name),
                "adopted": True,
            }
            adopted += 1
            reused = True
            jid = job_rec.get("job_id")
            print(f"[auto_throughput_eval] adopt axis={c.axis_type} x={xv} job_id={jid} state={st.state}")
            break

        if reused:
            continue

        todo.append((xv, log_path, cmd, job_name, point_key))

    save_manifest(manifest_path, manifest)
    job_db.flush()
    print(
        f"\n[auto_throughput_eval] case={c.name} axis={c.axis_type} grid={len(xs)} "
        f"cached={cached} adopted={adopted} todo={len(todo)} out={case_dir}"
    )

    if todo:
        patch_config_max_sim_num(cfg_abs, runtime_cfg, c.max_sim_num, max_error_num=c.max_err_num)

    def cancel_all(reason: str) -> None:
        for xv, info in list(in_flight.items()):
            try:
                ex.cancel(info["job"])
            except Exception:
                pass
            job_key = str(info.get("job_key") or "")
            if job_key:
                job_db.mark_cancel(job_key, reason=reason)
            print(f"[auto_throughput_eval] cancel axis={c.axis_type} x={xv} reason={reason}")
            in_flight.pop(xv, None)
        job_db.flush()
    try:
        submit_attempts: dict[float, int] = {}
        while todo or in_flight:
            while todo and len(in_flight) < max_in_flight:
                xv, log_path, cmd, job_name, point_key = todo.pop(0)
                if Path(log_path).exists():
                    try:
                        Path(log_path).unlink()
                    except Exception:
                        pass
                job = ex.submit(cmd, cwd=c.workdir, log_path=Path(log_path), job_name=job_name, queue=lsf.queue)
                backend, jid = _job_identity(job)
                job_key = f"{backend}:{jid}"
                t_submit = time.time()
                job_db.upsert(
                    {
                        "key": job_key,
                        "backend": backend,
                        "job_id": jid,
                        "job_name": job_name,
                        "point_key": point_key,
                        "stage": "throughput",
                        "axis_type": c.axis_type,
                        "axis_value": float(xv),
                        "log_path": str(Path(log_path).resolve()),
                        "cmd": list(cmd),
                        "submitted_at": t_submit,
                        "required_max_sim_num": int(c.max_sim_num),
                        "required_max_err_num": int(c.max_err_num),
                    }
                )
                job_db.set_state(job_key, "SUBMITTED", t=t_submit)
                job_db.flush()
                in_flight[xv] = {
                    "job": job,
                    "job_key": job_key,
                    "log_path": Path(log_path).resolve(),
                    "cmd": list(cmd),
                    "job_name": job_name,
                    "point_key": point_key,
                    "adopted": False,
                }
                jid = getattr(job, "job_id", "?")
                submit_attempts[xv] = int(submit_attempts.get(xv, 0)) + 1
                print(
                    f"[auto_throughput_eval] bsub axis={c.axis_type} x={xv} "
                    f"attempt={submit_attempts[xv]} job_id={jid} job_name={job_name}"
                )

            if not in_flight:
                break

            time.sleep(poll_sec)
            for xv, info in list(in_flight.items()):
                st = ex.poll(info["job"])
                job_key = str(info.get("job_key") or "")
                now = time.time()
                if job_key:
                    job_db.set_state(job_key, st.state, t=now)
                if not st.done:
                    continue

                in_flight.pop(xv, None)
                rec = parse_done_log_with_grace(
                    axis_type=c.axis_type,
                    axis_value=xv,
                    log_path=Path(info["log_path"]),
                    cmd=list(info["cmd"]),
                    grace_sec=max(float(runner.timeout_log_grace_sec), 30.0),
                )
                if rec is not None:
                    upsert_run(manifest, rec)
                    save_manifest(manifest_path, manifest)
                    if job_key:
                        job_db.set_state(job_key, "DONE", t=now)
                        job_db.upsert({"key": job_key, "finished_at": now})
                        job_db.flush()
                    print(f"[auto_throughput_eval] done axis={c.axis_type} x={xv} RBER={rec.raw_ber} aver_iter={rec.retry_dec_avg_iter}")
                    continue

                if job_key:
                    job_db.upsert({"key": job_key, "finished_at": now})
                    job_db.flush()

                retry_reason: Optional[str] = None
                if not st.ok:
                    retry_reason = (
                        f"job EXIT: axis={c.axis_type} x={xv} "
                        f"log={info['log_path']} job_name={info.get('job_name')}"
                    )
                else:
                    retry_reason = (
                        f"log incomplete: axis={c.axis_type} x={xv} "
                        f"log={info['log_path']} job_name={info.get('job_name')}"
                    )

                print(f"[auto_throughput_eval] {retry_reason}")
                if runner.fail_fast:
                    cancel_all(reason=retry_reason)
                    save_manifest(manifest_path, manifest)
                    raise SystemExit(retry_reason)

                next_attempt = int(submit_attempts.get(xv, 0)) + 1
                if next_attempt > runner.max_retries:
                    print(
                        f"[auto_throughput_eval] SKIP axis={c.axis_type} x={xv}: "
                        f"exceeded max_retries={runner.max_retries}"
                    )
                else:
                    point_key = str(info.get("point_key") or _throughput_point_key(c.axis_type, float(xv)))
                    todo.append(
                        (
                            float(xv),
                            Path(info["log_path"]).resolve(),
                            list(info["cmd"]),
                            str(info.get("job_name") or build_job_name(c, matrix_size=matrix_size, axis_value=float(xv))),
                            point_key,
                        )
                    )
                    print(f"[auto_throughput_eval] retry axis={c.axis_type} x={xv} next_attempt={next_attempt}")
                continue

    finally:
        job_db.flush()


def _throughput_point_key(axis_type: str, axis_value: float) -> str:
    return f"throughput:{axis_type}:{format_axis_tag(axis_value, axis_type=axis_type)}"


def _job_sort_key(rec: dict[str, Any]) -> float:
    t = rec.get("submitted_at")
    if t is None:
        t = rec.get("last_update")
    try:
        return float(t) if t is not None else 0.0
    except (TypeError, ValueError):
        return 0.0


def _load_job_records_by_point(job_db_path: Path, *, axis_type: str) -> dict[str, list[dict[str, Any]]]:
    if not job_db_path.exists():
        return {}
    try:
        data = json.loads(job_db_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    jobs = data.get("jobs", [])
    if not isinstance(jobs, list):
        return {}

    out: dict[str, list[dict[str, Any]]] = {}
    for rec in jobs:
        if not isinstance(rec, dict):
            continue
        if str(rec.get("axis_type") or "") != str(axis_type):
            continue
        point_key = str(rec.get("point_key") or "").strip()
        if not point_key:
            xv = rec.get("axis_value")
            try:
                point_key = _throughput_point_key(axis_type, float(xv))
            except (TypeError, ValueError):
                continue
        rec2 = dict(rec)
        rec2["point_key"] = point_key
        out.setdefault(point_key, []).append(rec2)

    for point_key in list(out.keys()):
        out[point_key].sort(key=_job_sort_key, reverse=True)
    return out


def _job_record_is_active(rec: dict[str, Any]) -> bool:
    state = str(rec.get("last_state") or "").upper()
    return state in {"SUBMITTED", "PEND", "RUN", "UNKNOWN", "PSUSP", "UNKWN"}


def _job_identity(job: Any) -> tuple[str, str]:
    if hasattr(job, "job_id"):
        return ("lsf", str(getattr(job, "job_id")))
    if hasattr(job, "proc") and hasattr(getattr(job, "proc"), "pid"):
        return ("local", str(getattr(getattr(job, "proc"), "pid")))
    return ("unknown", str(id(job)))


def _resume_executor_job(executor: Any, *, cmd: list[str], cwd: str, log_path: Path, rec: dict[str, Any]) -> Optional[Any]:
    resume_job = getattr(executor, "resume_job", None)
    if not callable(resume_job):
        return None
    return resume_job(
        cmd=list(cmd),
        cwd=str(cwd),
        log_path=log_path.resolve(),
        backend=str(rec.get("backend") or ""),
        job_id=str(rec.get("job_id") or ""),
    )


def _resolve_complete_record(
    *,
    axis_type: str,
    axis_value: float,
    manifest_rec: Optional[dict[str, Any]],
    fallback_log_path: Path,
    fallback_cmd: list[str],
) -> Optional[RunRecord]:
    candidates: list[tuple[Path, list[str]]] = []
    if manifest_rec is not None:
        rec_log_raw = manifest_rec.get("log_path")
        rec_log_path = Path(str(rec_log_raw)).resolve() if rec_log_raw else fallback_log_path.resolve()
        rec_cmd_raw = manifest_rec.get("cmd")
        rec_cmd = [str(x) for x in rec_cmd_raw] if isinstance(rec_cmd_raw, list) else list(fallback_cmd)
        candidates.append((rec_log_path, rec_cmd))

    fallback = (fallback_log_path.resolve(), list(fallback_cmd))
    seen: set[tuple[str, tuple[str, ...]]] = set()
    for log_path, cmd in candidates + [fallback]:
        key = (str(log_path), tuple(str(x) for x in cmd))
        if key in seen:
            continue
        seen.add(key)
        rec_done = parse_done_log(
            axis_type=axis_type,
            axis_value=axis_value,
            log_path=log_path,
            cmd=list(cmd),
            debug=False,
        )
        if rec_done is not None:
            return rec_done
    return None


def _text_has_final_statistics(text: str) -> bool:
    return re.search(r"^\s*\[STATISTICS\]\s+LDPC\s+FER\s*:\s*[0-9eE+\-\.]+\s*$", text, flags=re.M) is not None


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


def build_cmd(c: CaseConfig, axis_value: float, *, config_arg: Optional[str] = None) -> list[str]:
    cfg = str(config_arg) if config_arg is not None else c.config
    cmd: list[str] = [c.exe, c.sim_mode, cfg, c.ch_model]
    if c.ch_model.upper() != "CLEAN":
        if c.axis_type == "k":
            cmd.append(str(int(round(axis_value))))
        else:
            cmd.append(str(axis_value))
    cmd.extend(c.cmd_extra_args)
    return cmd


def parse_done_log(
    *,
    axis_type: str,
    axis_value: float,
    log_path: Path,
    cmd: list[str],
    debug: bool = True,
) -> Optional[RunRecord]:
    if not log_path.exists():
        if debug:
            print(f"[auto_throughput_eval] DEBUG: Log file does not exist: {log_path}")
        return None

    txt = log_path.read_text(encoding="utf-8", errors="replace")
    if not _text_has_final_statistics(txt):
        if debug:
            print(f"[auto_throughput_eval] DEBUG: Log missing final [STATISTICS]: {log_path}")
            print(f"[auto_throughput_eval] DEBUG: Log file size: {log_path.stat().st_size} bytes")
            print(f"[auto_throughput_eval] DEBUG: Last 500 chars of log:")
            print(txt[-500:] if len(txt) > 500 else txt)
        return None

    m = parse_log_text(txt)

    if m.raw_ber is None or m.retry_dec_avg_iter is None:
        if debug:
            print(
                f"[auto_throughput_eval] DEBUG: Parse incomplete - raw_ber={m.raw_ber}, retry_dec_avg_iter={m.retry_dec_avg_iter}"
            )
            print(f"[auto_throughput_eval] DEBUG: Log file size: {log_path.stat().st_size} bytes")
            print(f"[auto_throughput_eval] DEBUG: Last 500 chars of log:")
            print(txt[-500:] if len(txt) > 500 else txt)
        return None

    return RunRecord(
        axis_type=str(axis_type),
        axis_value=float(axis_value),
        log_path=str(log_path),
        cmd=cmd,
        raw_ber=m.raw_ber,
        theo_rber=m.theo_rber,
        ldpc_fer=m.ldpc_fer,
        fail_cw=m.fail_cw,
        total_packets=m.total_packets,
        retry_dec_avg_iter=m.retry_dec_avg_iter,
    )


def parse_done_log_with_grace(
    *,
    axis_type: str,
    axis_value: float,
    log_path: Path,
    cmd: list[str],
    grace_sec: float,
) -> Optional[RunRecord]:
    """
    Best-effort parse with a bounded grace window.

    Rationale: on shared filesystems (NFS), bjobs may report DONE/EXIT before the
    bsub -o log is created/flushed. We retry parsing for a short grace window.
    """
    t_end = time.time() + max(0.0, float(grace_sec))
    sleep_s = 0.2
    while True:
        rec = parse_done_log(axis_type=axis_type, axis_value=axis_value, log_path=log_path, cmd=cmd, debug=False)
        if rec is not None:
            return rec
        if time.time() >= t_end:
            break
        time.sleep(sleep_s)
        sleep_s = min(2.0, sleep_s * 1.5)
    # Final attempt with debug enabled.
    return parse_done_log(axis_type=axis_type, axis_value=axis_value, log_path=log_path, cmd=cmd, debug=True)


def log_path_for(c: CaseConfig, axis_value: float) -> Path:
    case_dir = Path(c.out_dir) / c.name
    return (case_dir / f"{c.log_prefix}_{format_axis_tag(axis_value, axis_type=c.axis_type)}.log").resolve()


def format_axis_tag(x: float, *, axis_type: str) -> str:
    if axis_type == "k":
        return f"k{int(round(x))}"
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return f"snr{s}"


def build_axis_grid(low: float, high: float, step: float, *, axis_type: str) -> list[float]:
    origin = float(low)
    x = quantize_axis(origin, step, axis_type=axis_type, origin=origin)
    out: list[float] = []
    while x <= high + 1e-12:
        out.append(float(x))
        x = quantize_axis(x + step, step, axis_type=axis_type, origin=origin)
        if len(out) > 1000000:
            raise RuntimeError("Axis grid too large; check step.")
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


def quantize_axis(x: float, step: float, *, axis_type: str, origin: float = 0.0) -> float:
    if step <= 0:
        return x
    if axis_type == "k":
        step_i = max(1, int(round(step)))
        q = origin + round((x - origin) / float(step_i)) * float(step_i)
        return float(int(round(q)))
    q = origin + round((x - origin) / step) * step
    q = round(q, step_decimals(step))
    if abs(q) < 1e-12:
        q = 0.0
    return float(q)


def log_path_for_submission(c: CaseConfig, axis_value: float, *, lsf: LsfConfig) -> Path:
    """
    Determine the log path that LSF should write to.

    If lsf.log_base_dir is set, redirect logs to:
      <log_base_dir>/<case_name>/<log_file>
    Otherwise use:
      <out_dir>/<case_name>/<log_file>
    """
    log_path = log_path_for(c, axis_value)
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


def build_job_name(c: CaseConfig, *, matrix_size: str, axis_value: float) -> str:
    x_s = format_axis_value(axis_value, axis_type=c.axis_type)
    # Simplified job name format: case_name + axis value
    return f"{c.name} {c.axis_type}-{x_s}"


def format_axis_value(x: float, *, axis_type: str) -> str:
    if axis_type == "k":
        return str(int(round(x)))
    return f"{x:.3f}".rstrip("0").rstrip(".")


def build_bsub_command_line(
    cmd: list[str],
    *,
    cwd_abs: str,
    log_path: str,
    job_name: str,
    queue: str,
    bsub_extra: list[str],
    include_cwd: bool = True,
) -> str:
    parts: list[str] = ["bsub"]
    if queue:
        parts.extend(["-q", quote_shell(queue)])
    if job_name:
        parts.extend(["-J", quote_shell(job_name)])
    parts.extend(["-o", quote_shell(log_path)])
    if include_cwd:
        parts.extend(["-cwd", quote_shell(cwd_abs)])
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
    Export axis/RBER/aver_iter to a CSV under the case directory.
    The data source is manifest.json (filled from parsing logs if needed).
    """
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    xs = list(build_axis_grid(c.x_low, c.x_high, c.x_step, axis_type=c.axis_type))
    rows: list[list[object]] = []
    missing: list[float] = []

    for xv in xs:
        rec = find_run(manifest, c.axis_type, xv)
        log_p = log_path_for_submission(c, xv, lsf=lsf)
        cmd = build_cmd(c, xv)
        rec_done = _resolve_complete_record(
            axis_type=c.axis_type,
            axis_value=xv,
            manifest_rec=rec,
            fallback_log_path=log_p,
            fallback_cmd=cmd,
        )
        raw_ber = rec_done.raw_ber if rec_done is not None else None
        avg_it = rec_done.retry_dec_avg_iter if rec_done is not None else None
        if rec_done is not None:
            upsert_run(manifest, rec_done)

        if raw_ber is None or avg_it is None:
            missing.append(float(xv))
        rows.append([float(xv), None if raw_ber is None else float(raw_ber), None if avg_it is None else float(avg_it)])

    save_manifest(manifest_path, manifest)

    csv_path = (case_dir / "throughput.csv").resolve()
    write_csv(csv_path, headers=[c.axis_type, "RBER", "aver_iter"], rows=rows)
    print(f"[auto_throughput_eval] export csv: {csv_path}")

    if missing:
        miss_s = ", ".join([format_axis_value(x, axis_type=c.axis_type) for x in missing[:20]])
        more = "" if len(missing) <= 20 else f" ...(+{len(missing)-20})"
        raise SystemExit(f"Missing metrics for {c.axis_type} points: {miss_s}{more}. See {manifest_path}")

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
