#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Optional

try:
    import tomllib  # Python 3.11+
except Exception:  # pragma: no cover
    tomllib = None  # type: ignore

# Keep imports runnable via:
#   python3 auto_fer_eval/auto_fer_eval.py ...
# without requiring package installation.
try:
    from log_parser import parse_log_text  # type: ignore
    from manifest import RunRecord, find_run, load_manifest, save_manifest, upsert_run  # type: ignore
    from planner import IntervalPlan, StepPolicy  # type: ignore
    from executors import DryRunExecutor, LocalExecutor, LsfExecutor  # type: ignore
    from job_db import JobDB  # type: ignore
    from xlsx_min import write_xlsx  # type: ignore
except ImportError:  # pragma: no cover
    from .log_parser import parse_log_text  # type: ignore
    from .manifest import RunRecord, find_run, load_manifest, save_manifest, upsert_run  # type: ignore
    from .planner import IntervalPlan, StepPolicy  # type: ignore
    from .executors import DryRunExecutor, LocalExecutor, LsfExecutor  # type: ignore
    from .job_db import JobDB  # type: ignore
    from .xlsx_min import write_xlsx  # type: ignore


def _abs_time_str() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def _status(msg: str) -> None:
    print(f"{msg} [{_abs_time_str()}]")


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
    x_start: float
    x_stop: float
    direction: int
    pilot_step: float
    low_margin: float
    snr_span_cap: float

    fer_hi: float
    fer_lo: float
    fer_pilot_stop: float
    fer_pilot_too_low: float
    max_points: int
    slope_max_decades: float
    step_policy: StepPolicy
    fit_enable: bool
    fit_fer_hi: float
    fit_fer_lo: float
    fit_target_fer: float
    fit_stop_margin: float
    fit_min_points: int
    fit_max_extend: float
    gate_enable: bool
    gate_zero_k: int


@dataclass(frozen=True)
class LoopConfig:
    name: str
    max_sim_num: Optional[int] = None
    fer_pilot_stop: Optional[float] = None
    snr_span_cap: Optional[float] = None
    fer_lo: Optional[float] = None
    fit_target_fer: Optional[float] = None
    fit_enable: Optional[bool] = None


@dataclass(frozen=True)
class AdaptiveConfig:
    enable: bool = False
    executor: str = "local"  # local | lsf

    # Stage budgets (patching "maximum simulation number" in *.cnfg)
    pilot_max_sim_num: int = 100
    main_max_sim_num: int = 1000

    # Main scan (None = use step_default from CaseConfig)
    main_step: Optional[float] = None
    # Trigger FER (None = use fit_fer_lo from CaseConfig)
    trigger_fer: Optional[float] = None

    # Scheduling
    max_in_flight: int = 20
    poll_sec: float = 30.0
    kill_margin: float = 0.0
    # Progress export (xlsx): during deep scan, export every N seconds (0 disables).
    export_xlsx_sec: float = 3600.0

    # Circuit breaker / fail-fast:
    # - startup_ok_required: only allow full concurrency after this many *successful*
    #   points are parsed (prevents mass submission when args/matrix path is wrong).
    # - fail_fast: abort on the first invalid log / failed job and cancel in-flight jobs.
    startup_ok_required: int = 1
    fail_fast: bool = True

    # Reliability guard for adaptive decisions (trigger / fit):
    # if FAIL CW is parsed and is < this value, treat the point as "too uncertain"
    # for decisions (but still record it in manifest). This mainly matters for
    # partial logs from timeouts / early termination.
    min_fail_cw_for_decision: int = 2

    # Timeout (seconds). If >0, kill jobs that run longer than this (RUN-time only for LSF).
    max_job_runtime_sec: float = 0.0
    # After bkill/terminate, wait and retry log parsing for this many seconds (best effort).
    timeout_log_grace_sec: float = 10.0

    # Optional initial seed SNR (otherwise uses case.snr_min)
    snr_init: Optional[float] = None

    # Backfill step (None = use step_mid from CaseConfig). Set to 0 to disable backfill.
    backfill_step: Optional[float] = None
    # Maximum FER decade jump before triggering backfill (e.g., 2.0 means 100x jump).
    backfill_decade_threshold: float = 2.0


@dataclass(frozen=True)
class LsfConfig:
    queue: str = ""
    queue_slow: str = ""  # For pilot/main/backfill (fallback to queue)
    queue_fast: str = ""  # For deep scan (fallback to queue)
    log_base_dir: str = ""  # LSF log output base directory (preserves case structure)
    bsub_extra: list[str] = field(default_factory=list)
    bjobs_extra: list[str] = field(default_factory=list)
    bkill_extra: list[str] = field(default_factory=list)


@dataclass
class PlanLimiter:
    max_new: int
    new_count: int = 0

    def on_new_cmd(self) -> None:
        if self.max_new <= 0:
            return
        self.new_count += 1
        if self.new_count >= self.max_new:
            raise PlanLimitReached


class PlanLimitReached(Exception):
    pass


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(prog="auto_fer_eval", description="Auto pilot+refine runner for FER interval.")
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--config", required=True, help="TOML config path")
    common.add_argument("--dry-run", action="store_true", help="Do not execute commands; only print the plan")
    common.add_argument(
        "--max-new",
        type=int,
        default=0,
        help="In --dry-run mode, stop after printing this many new commands (0 = no limit).",
    )
    common.add_argument(
        "--loop",
        default="",
        help="Only run a single loop by name (e.g., L0/L1/L2). Default: run all loops from config.",
    )

    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("run", parents=[common], help="Run pilot+refine for all cases")
    sub.add_parser("pilot", parents=[common], help="Pilot only (locate interval)")
    sub.add_parser("adaptive", parents=[common], help="Adaptive scheduling (pilot->scan->kill->deep)")

    ns = ap.parse_args(argv)
    limiter = PlanLimiter(max_new=int(ns.max_new)) if ns.dry_run and int(ns.max_new) > 0 else None

    cfg = load_config(Path(ns.config))
    cases = cfg["cases"]
    loops_all = cfg["loops"]
    adaptive_cfg: AdaptiveConfig = cfg.get("adaptive", AdaptiveConfig())
    lsf_cfg: LsfConfig = cfg.get("lsf", LsfConfig())

    if ns.cmd == "adaptive":
        if not adaptive_cfg.enable:
            raise SystemExit("Config missing [adaptive] enable=true (or set enable=true).")
        _run_adaptive_all_cases(
            cases,
            adaptive=adaptive_cfg,
            lsf=lsf_cfg,
            dry_run=bool(ns.dry_run),
            limiter=limiter,
        )
        return 0
    loops = loops_all
    if str(ns.loop).strip():
        want = str(ns.loop).strip()
        loops = [l for l in loops_all if l.name == want]
        if not loops:
            names = ", ".join([l.name for l in loops_all if l.name] or ["<default>"])
            raise SystemExit(f"Unknown loop '{want}'. Available: {names}")
    name_to_loop_idx = {l.name: i for i, l in enumerate(loops_all)}

    for base_case in cases:
        # Reuse the pilot interval from the first loop as the seed plan for later loops.
        ref_interval: Optional[IntervalPlan] = None
        if loops_all:
            ref_loop_name = loops_all[0].name
            ref_dir = Path(base_case.out_dir) / base_case.name
            ref_path = (ref_dir / ref_loop_name / "interval.json") if ref_loop_name else (ref_dir / "interval.json")
            if ref_path.exists():
                ref_interval = load_interval(ref_path)

        for loop in loops:
            prev_case_dir: Optional[Path] = None
            idx = name_to_loop_idx.get(loop.name)
            if idx is not None and idx > 0:
                prev_loop = loops_all[idx - 1]
                prev_case_dir = (Path(base_case.out_dir) / base_case.name / prev_loop.name).resolve()

            c = prepare_case_for_loop(base_case, loop)
            print()
            _status(
                f"[auto_fer_eval] case={c.name} model={c.ch_model} axis={c.axis_type} "
                f"start={c.x_start} stop={c.x_stop} dir={c.direction}"
            )
            try:
                if ns.cmd == "pilot":
                    plan = _pilot_only(
                        c,
                        dry_run=ns.dry_run,
                        limiter=limiter,
                        seed_plan=ref_interval,
                        prev_case_dir=prev_case_dir,
                    )
                else:
                    plan = _run_case(
                        c,
                        dry_run=ns.dry_run,
                        limiter=limiter,
                        seed_plan=ref_interval,
                        prev_case_dir=prev_case_dir,
                    )
                # Update ref_interval after the first loop finishes (or when no ref file existed).
                if ref_interval is None:
                    ref_interval = plan
            except PlanLimitReached:
                _status("[auto_fer_eval] reached --max-new limit; stop planning for now.")
                return 0
    return 0


def load_config(path: Path) -> dict[str, Any]:
    if tomllib is None:
        raise RuntimeError("tomllib not available in this Python; use Python>=3.11")
    data = tomllib.loads(path.read_text(encoding="utf-8"))

    defaults = data.get("defaults", {})
    loops_raw = data.get("loops", [])
    case_list = data.get("cases", [])
    adaptive_raw = data.get("adaptive", {})
    lsf_raw = data.get("lsf", {})
    if not isinstance(case_list, list) or not case_list:
        raise ValueError("Config must contain [[cases]] list.")

    cases: list[CaseConfig] = []
    for raw in case_list:
        merged = dict(defaults)
        merged.update(raw or {})
        cases.append(_parse_case(merged))

    loops = _parse_loops(loops_raw)
    adaptive = _parse_adaptive(adaptive_raw)
    lsf = _parse_lsf(lsf_raw)
    return {"cases": cases, "loops": loops, "adaptive": adaptive, "lsf": lsf}


def _parse_adaptive(raw: Any) -> AdaptiveConfig:
    if not raw:
        return AdaptiveConfig(enable=False)
    if not isinstance(raw, dict):
        raise ValueError("Config [adaptive] must be a table.")

    enable = bool(raw.get("enable", True))
    executor = str(raw.get("executor", "local")).strip().lower()
    if executor not in {"local", "lsf"}:
        raise ValueError("adaptive.executor must be 'local' or 'lsf'.")

    snr_init = raw.get("snr_init", None)
    snr_init_f = float(snr_init) if snr_init is not None else None

    main_step = raw.get("main_step", None)
    main_step_f = float(main_step) if main_step is not None else None

    trigger_fer = raw.get("trigger_fer", None)
    trigger_fer_f = float(trigger_fer) if trigger_fer is not None else None

    backfill_step = raw.get("backfill_step", None)
    backfill_step_f = float(backfill_step) if backfill_step is not None else None

    export_xlsx_sec = raw.get("export_xlsx_sec", None)
    export_xlsx_sec_f = float(export_xlsx_sec) if export_xlsx_sec is not None else 3600.0

    return AdaptiveConfig(
        enable=enable,
        executor=executor,
        pilot_max_sim_num=int(raw.get("pilot_max_sim_num", 100)),
        main_max_sim_num=int(raw.get("main_max_sim_num", 1000)),
        main_step=main_step_f,
        trigger_fer=trigger_fer_f,
        max_in_flight=int(raw.get("max_in_flight", 20)),
        poll_sec=float(raw.get("poll_sec", 30.0)),
        kill_margin=float(raw.get("kill_margin", 0.0)),
        export_xlsx_sec=export_xlsx_sec_f,
        startup_ok_required=int(raw.get("startup_ok_required", 1)),
        fail_fast=bool(raw.get("fail_fast", True)),
        min_fail_cw_for_decision=int(raw.get("min_fail_cw_for_decision", 2)),
        max_job_runtime_sec=float(raw.get("max_job_runtime_sec", 0.0)),
        timeout_log_grace_sec=float(raw.get("timeout_log_grace_sec", 10.0)),
        snr_init=snr_init_f,
        backfill_step=backfill_step_f,
        backfill_decade_threshold=float(raw.get("backfill_decade_threshold", 2.0)),
    )


def _parse_lsf(raw: Any) -> LsfConfig:
    if not raw:
        return LsfConfig()
    if not isinstance(raw, dict):
        raise ValueError("Config [lsf] must be a table.")
    queue = str(raw.get("queue", "")).strip()
    queue_slow = str(raw.get("queue_slow", "")).strip()
    queue_fast = str(raw.get("queue_fast", "")).strip()
    log_base_dir = str(raw.get("log_base_dir", "")).strip()
    bsub_extra = [str(x) for x in raw.get("bsub_extra", [])]
    bjobs_extra = [str(x) for x in raw.get("bjobs_extra", [])]
    bkill_extra = [str(x) for x in raw.get("bkill_extra", [])]
    return LsfConfig(
        queue=queue,
        queue_slow=queue_slow,
        queue_fast=queue_fast,
        log_base_dir=log_base_dir,
        bsub_extra=bsub_extra,
        bjobs_extra=bjobs_extra,
        bkill_extra=bkill_extra,
    )


def _parse_loops(raw: Any) -> list[LoopConfig]:
    if not raw:
        # Backward-compatible default: no loop directory suffix.
        return [LoopConfig(name="", max_sim_num=None)]
    if not isinstance(raw, list):
        raise ValueError("Config loops must be a list ([[loops]]).")
    out: list[LoopConfig] = []
    for item in raw:
        if not isinstance(item, dict):
            raise ValueError("Each loop must be a table ([[loops]]).")
        name = str(item.get("name") or item.get("id") or "").strip()
        if not name:
            raise ValueError("Loop must have a non-empty name (use name=\"L0\" etc).")
        max_sim_num = item.get("max_sim_num")
        max_sim_num_i = int(max_sim_num) if max_sim_num is not None else None
        fer_pilot_stop = item.get("fer_pilot_stop")
        fer_pilot_stop_f = float(fer_pilot_stop) if fer_pilot_stop is not None else None
        snr_span_cap = item.get("snr_span_cap")
        snr_span_cap_f = float(snr_span_cap) if snr_span_cap is not None else None
        fer_lo = item.get("fer_lo")
        fer_lo_f = float(fer_lo) if fer_lo is not None else None
        fit_target_fer = item.get("fit_target_fer")
        fit_target_fer_f = float(fit_target_fer) if fit_target_fer is not None else None
        fit_enable = item.get("fit_enable")
        fit_enable_b = bool(fit_enable) if fit_enable is not None else None
        out.append(
            LoopConfig(
                name=name,
                max_sim_num=max_sim_num_i,
                fer_pilot_stop=fer_pilot_stop_f,
                snr_span_cap=snr_span_cap_f,
                fer_lo=fer_lo_f,
                fit_target_fer=fit_target_fer_f,
                fit_enable=fit_enable_b,
            )
        )
    return out


def _parse_case(d: dict[str, Any]) -> CaseConfig:
    name = str(d["name"])
    workdir = str(d.get("workdir", "."))
    exe = str(d["exe"])
    sim_mode = str(d.get("sim_mode", "LDPC"))
    config = str(d["config"])
    ch_model = str(d["ch_model"])
    cmd_extra_args_raw = d.get("cmd_extra_args", [])
    if cmd_extra_args_raw is None:
        cmd_extra_args_raw = []
    if isinstance(cmd_extra_args_raw, str):
        raise ValueError(
            "cmd_extra_args must be a TOML array, e.g. "
            'cmd_extra_args = ["0", "./matrix"] (matrix_id, matrix_dir).'
        )
    if not isinstance(cmd_extra_args_raw, list):
        raise ValueError("cmd_extra_args must be a TOML array (list of strings).")
    cmd_extra_args = [str(x) for x in cmd_extra_args_raw]
    out_dir = str(d.get("out_dir", "./perf_auto"))
    log_prefix = str(d.get("log_prefix", "run"))

    # Axis type: snr for AWGN-like; k for ERR_INJ.
    axis_type = str(d.get("axis_type", "")).strip().lower()
    if not axis_type:
        axis_type = "k" if ch_model.upper() == "ERR_INJ" else "snr"
    if axis_type not in {"snr", "k"}:
        raise ValueError(f"Unsupported axis_type: {axis_type}")

    if axis_type == "snr":
        x_start = float(d["snr_min"])
        x_stop = float(d["snr_max"])
        direction = +1
    else:
        k_min = float(d["k_min"])
        k_max = float(d["k_max"])
        x_start = k_max
        x_stop = k_min
        direction = -1

    pilot_step = float(d.get("pilot_step", 0.5))
    low_margin = float(d.get("low_margin", 0.1))
    snr_span_cap = float(d.get("snr_span_cap", 0.5))
    fer_hi = float(d.get("fer_hi", 0.9))
    fer_lo = float(d.get("fer_lo", 1e-6))
    fer_pilot_stop = float(d.get("fer_pilot_stop", 1e-3))
    fer_pilot_too_low = float(d.get("fer_pilot_too_low", 0.1))
    max_points = int(d.get("max_points", 300))
    slope_max_decades = float(d.get("slope_max_decades", 0.5))

    step_policy = StepPolicy(
        step_default=float(d.get("step_default", 0.1)),
        step_mid=float(d.get("step_mid", 0.05)),
        step_low=float(d.get("step_low", 0.025)),
        fer_mid=float(d.get("fer_mid", 1e-5)),
        fer_low=float(d.get("fer_low", 1e-6)),
    )

    fit_enable = bool(d.get("fit_enable", False))
    fit_fer_hi = float(d.get("fit_fer_hi", 1e-2))
    fit_fer_lo = float(d.get("fit_fer_lo", 1e-4))
    fit_target_fer = float(d.get("fit_target_fer", fer_lo))
    fit_stop_margin = float(d.get("fit_stop_margin", 0.05))
    fit_min_points = int(d.get("fit_min_points", 3))
    fit_max_extend = float(d.get("fit_max_extend", 1.0))
    gate_enable = bool(d.get("gate_enable", True))
    gate_zero_k = int(d.get("gate_zero_k", 2))

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
        x_start=x_start,
        x_stop=x_stop,
        direction=direction,
        pilot_step=pilot_step,
        low_margin=low_margin,
        snr_span_cap=snr_span_cap,
        fer_hi=fer_hi,
        fer_lo=fer_lo,
        fer_pilot_stop=fer_pilot_stop,
        fer_pilot_too_low=fer_pilot_too_low,
        max_points=max_points,
        slope_max_decades=slope_max_decades,
        step_policy=step_policy,
        fit_enable=fit_enable,
        fit_fer_hi=fit_fer_hi,
        fit_fer_lo=fit_fer_lo,
        fit_target_fer=fit_target_fer,
        fit_stop_margin=fit_stop_margin,
        fit_min_points=fit_min_points,
        fit_max_extend=fit_max_extend,
        gate_enable=gate_enable,
        gate_zero_k=gate_zero_k,
    )


def prepare_case_for_loop(base: CaseConfig, loop: LoopConfig) -> CaseConfig:
    # Create a sub-case name to avoid collisions across loops.
    name = base.name if not loop.name else f"{base.name}/{loop.name}"
    out_dir = base.out_dir

    fer_pilot_stop = loop.fer_pilot_stop if loop.fer_pilot_stop is not None else base.fer_pilot_stop
    snr_span_cap = loop.snr_span_cap if loop.snr_span_cap is not None else base.snr_span_cap
    fer_lo = loop.fer_lo if loop.fer_lo is not None else base.fer_lo
    fit_enable = loop.fit_enable if loop.fit_enable is not None else base.fit_enable
    fit_target_fer = base.fit_target_fer
    if loop.fit_target_fer is not None:
        fit_target_fer = loop.fit_target_fer
    elif loop.fer_lo is not None:
        # If a loop overrides fer_lo but not fit_target_fer, default the fit target to the loop's fer_lo.
        fit_target_fer = fer_lo

    cfg_path = base.config
    if loop.max_sim_num is not None:
        case_dir = Path(out_dir) / name
        case_dir.mkdir(parents=True, exist_ok=True)

        base_cfg = Path(cfg_path)
        if not base_cfg.is_absolute():
            base_cfg = (Path(base.workdir) / base_cfg).resolve()

        loop_cfg = (case_dir / f"config_{loop.name}.cnfg").resolve()
        patch_config_max_sim_num(base_cfg, loop_cfg, loop.max_sim_num)
        cfg_path = str(loop_cfg)

    return CaseConfig(
        name=name,
        workdir=base.workdir,
        exe=base.exe,
        sim_mode=base.sim_mode,
        config=cfg_path,
        ch_model=base.ch_model,
        cmd_extra_args=list(base.cmd_extra_args),
        out_dir=out_dir,
        log_prefix=base.log_prefix,
        axis_type=base.axis_type,
        x_start=base.x_start,
        x_stop=base.x_stop,
        direction=base.direction,
        pilot_step=base.pilot_step,
        low_margin=base.low_margin,
        snr_span_cap=snr_span_cap,
        fer_hi=base.fer_hi,
        fer_lo=fer_lo,
        fer_pilot_stop=fer_pilot_stop,
        fer_pilot_too_low=base.fer_pilot_too_low,
        max_points=base.max_points,
        slope_max_decades=base.slope_max_decades,
        step_policy=base.step_policy,
        fit_enable=fit_enable,
        fit_fer_hi=base.fit_fer_hi,
        fit_fer_lo=base.fit_fer_lo,
        fit_target_fer=fit_target_fer,
        fit_stop_margin=base.fit_stop_margin,
        fit_min_points=base.fit_min_points,
        fit_max_extend=base.fit_max_extend,
        gate_enable=base.gate_enable,
        gate_zero_k=base.gate_zero_k,
    )


def patch_config_max_sim_num(src: Path, dst: Path, max_sim_num: int, max_error_num: Optional[int] = None) -> None:
    """
    Patch the config file's "maximum simulation number" and optionally "maximum error number" lines.
    If the marker comments are not found, fallback to line 9 and 10 (1-based).

    Args:
        src: Source config file
        dst: Destination config file
        max_sim_num: Maximum simulation number (line 9)
        max_error_num: Maximum error number (line 10). If None, keep original value.
    """
    txt = src.read_text(encoding="utf-8")
    lines = txt.splitlines()

    # Find and patch maximum simulation number (line 9)
    idx_sim = None
    for i, line in enumerate(lines):
        if "maximum simulation number" in line:
            idx_sim = i
            break
    if idx_sim is None:
        idx_sim = 8  # line 9 (0-based)

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

    # Find and patch maximum error number (line 10) if specified
    if max_error_num is not None:
        idx_err = None
        for i, line in enumerate(lines):
            if "maximum error number" in line:
                idx_err = i
                break
        if idx_err is None:
            idx_err = 9  # line 10 (0-based)

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


def _pilot_only(
    c: CaseConfig,
    dry_run: bool,
    limiter: Optional[PlanLimiter],
    seed_plan: Optional[IntervalPlan],
    prev_case_dir: Optional[Path],
) -> IntervalPlan:
    case_dir = Path(c.out_dir) / c.name
    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    try:
        plan0 = seed_plan if seed_plan is not None else pilot_find_interval(c, manifest, dry_run=dry_run, limiter=limiter)
        plan = gate_plan_start_from_prev_zeros(c, plan0, prev_case_dir)
        _status(f"[auto_fer_eval] pilot interval: {c.axis_type} start={plan.start} stop={plan.stop} dir={plan.direction}")
        save_interval(case_dir / "interval.json", plan)
    finally:
        save_manifest(manifest_path, manifest)
    return plan


def _run_case(
    c: CaseConfig,
    dry_run: bool,
    limiter: Optional[PlanLimiter],
    seed_plan: Optional[IntervalPlan],
    prev_case_dir: Optional[Path],
) -> IntervalPlan:
    case_dir = Path(c.out_dir) / c.name
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    try:
        base_plan0 = seed_plan if seed_plan is not None else pilot_find_interval(c, manifest, dry_run=dry_run, limiter=limiter)
        base_plan = gate_plan_start_from_prev_zeros(c, base_plan0, prev_case_dir)
        stop_x = base_plan.stop
        plan = base_plan
        _status(f"[auto_fer_eval] interval: {c.axis_type} start={plan.start} stop={stop_x} dir={plan.direction}")
        save_interval(case_dir / "interval.json", IntervalPlan(start=plan.start, stop=stop_x, axis_type=plan.axis_type, direction=plan.direction))

        # Refine inside interval.
        x = plan.start
        origin = plan.start
        n_points = 0
        last_fer_eff: Optional[float] = None
        prev_x: Optional[float] = None
        prev_fer_eff: Optional[float] = None

        while _in_range(x, stop_x, plan.direction):
            n_points += 1
            if n_points > c.max_points:
                _status(f"[auto_fer_eval] reached max_points={c.max_points}, stop.")
                break

            rec = run_point(c, manifest, x, dry_run=dry_run, limiter=limiter)
            last_fer_eff = fer_effective(rec.ldpc_fer, rec.total_packets)

            if last_fer_eff is not None and last_fer_eff <= c.fer_lo:
                _status(f"[auto_fer_eval] reached fer_lo={c.fer_lo:g}, stop.")
                break

            if c.fit_enable and c.axis_type == "snr":
                est = predict_stop_snr_from_manifest(c, manifest)
                if est is not None:
                    # Only extend, never shrink.
                    new_stop = min(c.x_stop, est)
                    if c.fit_max_extend > 0:
                        new_stop = min(new_stop, stop_x + c.fit_max_extend)
                    if new_stop > stop_x + 1e-12:
                        stop_x = new_stop
                        plan = IntervalPlan(start=plan.start, stop=stop_x, axis_type=plan.axis_type, direction=plan.direction)
                        save_interval(case_dir / "interval.json", plan)

            step = choose_step(
                c,
                cur_x=x,
                cur_fer_eff=last_fer_eff,
                prev_x=prev_x,
                prev_fer_eff=prev_fer_eff,
            )

            prev_x = x
            prev_fer_eff = last_fer_eff
            x = _quantize_axis(x + plan.direction * step, step, axis_type=c.axis_type, origin=origin)

        # Ensure the high endpoint is present.
        end_x = stop_x
        end_x = _quantize_axis(end_x, c.step_policy.step_low, axis_type=c.axis_type, origin=end_x)
        if find_run(manifest, c.axis_type, end_x) is None:
            run_point(c, manifest, end_x, dry_run=dry_run, limiter=limiter)
    finally:
        save_manifest(manifest_path, manifest)
    return IntervalPlan(start=plan.start, stop=stop_x, axis_type=plan.axis_type, direction=plan.direction)


def _run_adaptive_all_cases(
    cases: list[CaseConfig],
    *,
    adaptive: AdaptiveConfig,
    lsf: LsfConfig,
    dry_run: bool,
    limiter: Optional[PlanLimiter],
) -> None:
    if dry_run:
        # Use DryRunExecutor to show what would be submitted
        ex = DryRunExecutor(
            queue_slow=lsf.queue_slow or lsf.queue,
            queue_fast=lsf.queue_fast or lsf.queue,
            log_base_dir=lsf.log_base_dir,
        )
    elif adaptive.executor == "lsf":
        ex = LsfExecutor(
            queue=lsf.queue,
            queue_slow=lsf.queue_slow,
            queue_fast=lsf.queue_fast,
            log_base_dir=lsf.log_base_dir,
            bsub_extra=lsf.bsub_extra,
            bjobs_extra=lsf.bjobs_extra,
            bkill_extra=lsf.bkill_extra,
        )
    else:
        ex = LocalExecutor()

    for base in cases:
        if base.axis_type != "snr":
            raise RuntimeError("adaptive mode currently supports axis_type=snr only.")
        _run_adaptive_case(base, adaptive=adaptive, executor=ex, lsf=lsf, dry_run=dry_run, limiter=limiter)


def _resolve_cfg_path(c: CaseConfig) -> Path:
    p = Path(c.config)
    if not p.is_absolute():
        p = (Path(c.workdir) / p).resolve()
    return p


def _run_adaptive_case(
    base: CaseConfig,
    *,
    adaptive: AdaptiveConfig,
    executor: Any,
    lsf: LsfConfig,
    dry_run: bool,
    limiter: Optional[PlanLimiter],
) -> None:
    """
    A single-case adaptive workflow aligned with your characterization goal:

      1) Pilot (small max_sim) bidirectionally to locate the start near FER≈fer_hi.
      2) Main scan with fixed step (e.g., 0.1 dB) and bounded concurrency.
      3) When any completed point reaches FER<=trigger_fer (e.g., 1e-4), kill all in-flight jobs with higher SNR.
      4) Optionally, use log-log fit (RAW_BER vs FER) to predict a deeper stop SNR, then continue with finer steps.
    """

    case_root = (Path(base.out_dir) / base.name / "adaptive").resolve()
    pilot_dir = case_root / "pilot"
    main_dir = case_root / "main"
    pilot_dir.mkdir(parents=True, exist_ok=True)
    main_dir.mkdir(parents=True, exist_ok=True)

    # Patch configs for pilot/main.
    base_cfg = _resolve_cfg_path(base)
    pilot_cfg = case_root / "config_pilot.cnfg"
    main_cfg = case_root / "config_main.cnfg"
    # Pilot: fixed packet count (max_error_num=0)
    patch_config_max_sim_num(base_cfg, pilot_cfg, adaptive.pilot_max_sim_num, max_error_num=0)
    # Main/Backfill/Deep: run until 10 failures (max_error_num=10)
    patch_config_max_sim_num(base_cfg, main_cfg, adaptive.main_max_sim_num, max_error_num=10)

    pilot_case = replace(
        base,
        name=str(Path(base.name) / "adaptive" / "pilot"),
        config=str(pilot_cfg),
        out_dir=str(base.out_dir),
        log_prefix=f"{base.log_prefix}_pilot",
    )
    main_case = replace(
        base,
        name=str(Path(base.name) / "adaptive" / "main"),
        config=str(main_cfg),
        out_dir=str(base.out_dir),
        log_prefix=base.log_prefix,
    )

    print()
    _status(
        f"[auto_fer_eval][adaptive] case={base.name} model={base.ch_model} "
        f"snr_range=[{base.x_start},{base.x_stop}] pilot_max={adaptive.pilot_max_sim_num} main_max={adaptive.main_max_sim_num}"
    )

    # ---- Stage 1: pilot (bidirectional)
    pilot_manifest_path = pilot_dir / "manifest.json"
    pilot_manifest = load_manifest(pilot_manifest_path)
    try:
        snr_init = adaptive.snr_init if adaptive.snr_init is not None else base.x_start
        start_snr = pilot_find_start_bidirectional(
            pilot_case,
            pilot_manifest,
            snr_init=snr_init,
            dry_run=dry_run,
            limiter=limiter,
        )
    finally:
        save_manifest(pilot_manifest_path, pilot_manifest)

    if dry_run:
        _status(f"[auto_fer_eval][adaptive][dry-run] pilot start_snr≈{start_snr}")
        # Continue to show what main scan would submit

    # ---- Resolve adaptive parameters with defaults from CaseConfig
    main_step = adaptive.main_step if adaptive.main_step is not None else main_case.step_policy.step_default
    trigger_fer = adaptive.trigger_fer if adaptive.trigger_fer is not None else main_case.fit_fer_lo
    backfill_step = adaptive.backfill_step if adaptive.backfill_step is not None else main_case.step_policy.step_mid

    # ---- Stage 2: main scan
    main_manifest_path = main_dir / "manifest.json"
    manifest = load_manifest(main_manifest_path)
    try:
        snr_trigger = adaptive_main_scan(
            main_case,
            manifest,
            executor=executor,
            start_snr=start_snr,
            step=main_step,
            trigger_fer=trigger_fer,
            max_in_flight=adaptive.max_in_flight,
            poll_sec=adaptive.poll_sec,
            kill_margin=adaptive.kill_margin,
            startup_ok_required=adaptive.startup_ok_required,
            fail_fast=adaptive.fail_fast,
            min_fail_cw_for_decision=adaptive.min_fail_cw_for_decision,
            max_job_runtime_sec=adaptive.max_job_runtime_sec,
            timeout_log_grace_sec=adaptive.timeout_log_grace_sec,
            queue=lsf.queue_slow,  # main scan uses slow queue
        )

        # ---- Stage 2.5: cleanup invalid logs from killed jobs
        # Points beyond snr_trigger with no fail_cw or fail_cw < threshold are invalid.
        # Delete their logs and remove from manifest so deep scan can re-run them.
        _cleanup_invalid_logs_after_main_scan(
            main_case,
            manifest,
            snr_trigger=snr_trigger,
            min_fail_cw=adaptive.min_fail_cw_for_decision,
        )

        # ---- Stage 3: backfill gaps where FER jumps > threshold
        # (Run before deep scan so backfill points can participate in fit)
        if backfill_step > 0 and adaptive.backfill_decade_threshold > 0:
            backfill_count = adaptive_backfill_gaps(
                main_case,
                manifest,
                executor=executor,
                backfill_step=backfill_step,
                decade_threshold=adaptive.backfill_decade_threshold,
                poll_sec=adaptive.poll_sec,
                fail_fast=adaptive.fail_fast,
                timeout_log_grace_sec=adaptive.timeout_log_grace_sec,
                max_in_flight=adaptive.max_in_flight,
                queue=lsf.queue_slow,  # backfill uses slow queue
            )
            if backfill_count > 0:
                _status(f"[auto_fer_eval][adaptive] backfill done: {backfill_count} points added")

        # ---- Stage 4: deep scan (optional, sequential)
        # Deep scan can start from:
        #   a) snr_trigger (if trigger was hit), or
        #   b) the lowest FER point in manifest (if no trigger but we have data)
        deep_start = snr_trigger
        if deep_start is None and main_case.fit_enable:
            # Find the point with lowest FER as starting point for deep scan
            deep_start = _find_lowest_fer_snr(manifest)
            if deep_start is not None:
                _status(f"[auto_fer_eval][adaptive] no trigger hit, using lowest FER point snr={deep_start} for deep scan")

        if deep_start is not None and main_case.fit_enable:
            est_stop = predict_stop_snr_from_manifest(
                main_case,
                manifest,
                min_fail_cw_for_decision=adaptive.min_fail_cw_for_decision,
            )
            if est_stop is not None:
                deep_stop = min(main_case.x_stop, est_stop)
                adaptive_deep_scan(
                    main_case,
                    manifest,
                    executor=executor,
                    start_snr=deep_start,
                    stop_snr=deep_stop,
                    poll_sec=adaptive.poll_sec,
                    fail_fast=adaptive.fail_fast,
                    min_fail_cw_for_decision=adaptive.min_fail_cw_for_decision,
                    # Per your workflow: only coarse(main scan) points are watchdog-timed.
                    # Deep (fine) points must be allowed to run to completion.
                    max_job_runtime_sec=0.0,
                    timeout_log_grace_sec=adaptive.timeout_log_grace_sec,
                    max_in_flight=adaptive.max_in_flight,
                    queue=lsf.queue_fast,  # deep scan uses fast queue
                    export_xlsx_sec=adaptive.export_xlsx_sec,
                )
            else:
                # Print reason for skipping deep scan
                fit_pts = _count_fit_points(main_case, manifest, adaptive.min_fail_cw_for_decision)
                _status(
                    f"[auto_fer_eval][adaptive] skip deep scan: fit failed "
                    f"(only {fit_pts} points in [{main_case.fit_fer_lo:g},{main_case.fit_fer_hi:g}], need {main_case.fit_min_points})"
                )

        # ---- Stage 4.5: finalize fit-window points (re-run incomplete logs to get [STATISTICS])
        # Per your requirement: temporarily allow partial logs for fit/trigger decisions,
        # but eventually make fit-window points complete. Only enforce for snr <= snr_trigger.
        if not dry_run:
            finalized = adaptive_finalize_fit_logs(
                main_case,
                manifest,
                executor=executor,
                snr_trigger=snr_trigger,
                poll_sec=adaptive.poll_sec,
                fail_fast=adaptive.fail_fast,
                min_fail_cw_for_decision=adaptive.min_fail_cw_for_decision,
                timeout_log_grace_sec=adaptive.timeout_log_grace_sec,
                max_in_flight=adaptive.max_in_flight,
                queue=lsf.queue_slow,  # finalize uses slow queue (same as main/backfill)
            )
            if finalized > 0:
                _status(f"[auto_fer_eval][adaptive] finalize done: {finalized} point(s) re-run to completion")
    finally:
        save_manifest(main_manifest_path, manifest)

    # ---- Print FER summary
    _print_fer_summary(manifest, case_name=base.name)


def pilot_find_start_bidirectional(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    snr_init: float,
    dry_run: bool,
    limiter: Optional[PlanLimiter],
) -> float:
    """
    Locate a starting SNR near the "FER≈fer_hi" boundary, but allow the search to start
    from an arbitrary seed SNR (snr_init), and move in either direction.

    The returned value is a conservative start (shifted by low_margin toward worse channel).
    """
    if c.axis_type != "snr":
        raise RuntimeError("pilot_find_start_bidirectional requires axis_type=snr.")

    # In dry-run mode, skip actual pilot search and use snr_init as start
    if dry_run:
        start = _quantize_axis(snr_init, c.step_policy.step_default, axis_type="snr", origin=snr_init)
        _status(f"[auto_fer_eval][adaptive][dry-run] pilot skipped, using snr_init={start} as start_snr")
        return float(start)

    origin = float(snr_init)
    x = _quantize_axis(origin, c.pilot_step, axis_type="snr", origin=origin)
    rec = run_point(c, manifest, x, dry_run=dry_run, limiter=limiter)
    if rec.ldpc_fer is None and not dry_run:
        raise RuntimeError(
            "pilot failed: cannot parse 'LDPC FER' from the first log. "
            f"Check your executable args / matrix_dir. log={rec.log_path}"
        )
    fer_eff = fer_effective(rec.ldpc_fer, rec.total_packets)
    if fer_eff is None:
        return max(c.x_start, min(c.x_stop, x - c.low_margin))

    found_low_at: Optional[float] = None

    if fer_eff <= c.fer_hi:
        # Too good -> go to worse SNR (decrease) until FER>fer_hi.
        x_good = x
        while True:
            x_bad = _quantize_axis(x_good - c.pilot_step, c.pilot_step, axis_type="snr", origin=origin)
            if x_bad < c.x_start - 1e-12:
                # Hit boundary; accept the best we can.
                found_low_at = x_good
                break
            rec_bad = run_point(c, manifest, x_bad, dry_run=dry_run, limiter=limiter)
            if rec_bad.ldpc_fer is None and not dry_run:
                raise RuntimeError(
                    "pilot failed: cannot parse 'LDPC FER' from log. "
                    f"Check your executable args / matrix_dir. log={rec_bad.log_path}"
                )
            fer_bad = fer_effective(rec_bad.ldpc_fer, rec_bad.total_packets)
            if fer_bad is not None and fer_bad > c.fer_hi:
                refine_steps = [c.step_policy.step_default, c.step_policy.step_mid, c.step_policy.step_low]
                found_low_at = refine_pilot_start(
                    c,
                    manifest,
                    x_bad=x_bad,
                    x_good=x_good,
                    fer_hi=c.fer_hi,
                    fer_too_low=c.fer_pilot_too_low,
                    steps=refine_steps,
                    dry_run=dry_run,
                    limiter=limiter,
                )
                break
            x_good = x_bad
    else:
        # Too bad -> go to better SNR (increase) until FER<=fer_hi.
        x_bad = x
        while True:
            x_good = _quantize_axis(x_bad + c.pilot_step, c.pilot_step, axis_type="snr", origin=origin)
            if x_good > c.x_stop + 1e-12:
                raise RuntimeError(f"pilot bidir failed: never reached FER<={c.fer_hi} within range.")
            rec_good = run_point(c, manifest, x_good, dry_run=dry_run, limiter=limiter)
            if rec_good.ldpc_fer is None and not dry_run:
                raise RuntimeError(
                    "pilot failed: cannot parse 'LDPC FER' from log. "
                    f"Check your executable args / matrix_dir. log={rec_good.log_path}"
                )
            fer_good = fer_effective(rec_good.ldpc_fer, rec_good.total_packets)
            if fer_good is not None and fer_good <= c.fer_hi:
                refine_steps = [c.step_policy.step_default, c.step_policy.step_mid, c.step_policy.step_low]
                found_low_at = refine_pilot_start(
                    c,
                    manifest,
                    x_bad=x_bad,
                    x_good=x_good,
                    fer_hi=c.fer_hi,
                    fer_too_low=c.fer_pilot_too_low,
                    steps=refine_steps,
                    dry_run=dry_run,
                    limiter=limiter,
                )
                break
            x_bad = x_good

    assert found_low_at is not None
    start = found_low_at - c.low_margin
    start = max(c.x_start, min(c.x_stop, start))
    start = _quantize_axis(start, c.step_policy.step_default, axis_type="snr", origin=start)
    _status(f"[auto_fer_eval][adaptive] pilot start: snr≈{start} (found_low_at={found_low_at})")
    return float(start)


def _log_path_for(c: CaseConfig, axis_value: float) -> Path:
    log_dir = Path(c.out_dir) / c.name
    axis_str = _format_axis(axis_value, c.axis_type)
    return (log_dir / f"{c.log_prefix}_{axis_str}.log").resolve()


def _parse_done_log(c: CaseConfig, axis_value: float, log_path: Path, cmd: list[str]) -> Optional[RunRecord]:
    if not log_path.exists():
        return None
    txt = log_path.read_text(encoding="utf-8", errors="replace")
    m = parse_log_text(txt)
    if m.ldpc_fer is None:
        return None
    return RunRecord(
        axis_type=c.axis_type,
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


def adaptive_main_scan(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    executor: Any,
    start_snr: float,
    step: float,
    trigger_fer: float,
    max_in_flight: int,
    poll_sec: float,
    kill_margin: float,
    startup_ok_required: int,
    fail_fast: bool,
    min_fail_cw_for_decision: int,
    max_job_runtime_sec: float,
    timeout_log_grace_sec: float,
    queue: str = "",
) -> Optional[float]:
    if c.axis_type != "snr":
        raise RuntimeError("adaptive_main_scan requires axis_type=snr.")

    max_in_flight = max(1, int(max_in_flight))
    startup_ok_required = max(1, int(startup_ok_required))
    min_fail_cw_for_decision = max(0, int(min_fail_cw_for_decision))
    poll_sec = max(0.2, float(poll_sec))

    origin = float(start_snr)
    step_cur = float(step)  # coarse step (fixed in main scan)
    x_next = _quantize_axis(origin, step_cur, axis_type="snr", origin=origin)
    x_stop_coarse = float(c.x_stop)

    in_flight: dict[float, dict[str, Any]] = {}
    snr_trigger: Optional[float] = None
    ok_count = 0
    min_fail = min_fail_cw_for_decision

    job_db = JobDB((Path(c.out_dir) / c.name / "jobs.json").resolve())

    def _job_identity(job: Any) -> tuple[str, str]:
        # LsfJob has job_id; LocalJob has proc.pid. Duck-typed to avoid imports here.
        if hasattr(job, "job_id"):
            return ("lsf", str(getattr(job, "job_id")))
        if hasattr(job, "proc") and hasattr(getattr(job, "proc"), "pid"):
            return ("local", str(getattr(getattr(job, "proc"), "pid")))
        return ("unknown", str(id(job)))

    def _is_reliable_for_decision(fail_cw: Optional[int]) -> bool:
        if min_fail <= 0:
            return True
        if fail_cw is None:
            return True
        return int(fail_cw) >= min_fail

    def _is_reliable_for_trigger(fail_cw: Optional[int], fer_eff: Optional[float]) -> bool:
        """
        Trigger is a *stop* decision, so we allow a conservative 0-fail upper bound:
        if fer_eff already meets the trigger threshold, accept even when FAIL CW < K.
        """
        if fer_eff is not None and math.isfinite(fer_eff) and fer_eff <= trigger_fer:
            return True
        return _is_reliable_for_decision(fail_cw)

    def submit_one(xv: float) -> None:
        nonlocal in_flight, snr_trigger, ok_count
        log_path = _log_path_for(c, xv)
        cmd = build_cmd(c, xv)

        # First check manifest for existing record (higher priority than log file).
        cached = find_run(manifest, "snr", xv)
        if cached and cached.get("ldpc_fer") is not None:
            ok_count += 1
            fer_eff = fer_effective(float(cached["ldpc_fer"]), cached.get("total_packets"))
            # Print cached result
            fcw_s = str(cached.get("fail_cw")) if cached.get("fail_cw") is not None else "?"
            tp_s = str(cached.get("total_packets")) if cached.get("total_packets") is not None else "?"
            fer_s = f"{cached['ldpc_fer']:.6g}"
            _status(f"[auto_fer_eval][adaptive] cached: snr={xv} FER={fer_s} FAIL_CW={fcw_s} N={tp_s}")
            if fer_eff is not None and fer_eff <= trigger_fer:
                if _is_reliable_for_trigger(cached.get("fail_cw"), fer_eff):
                    if snr_trigger is None or xv < snr_trigger:
                        snr_trigger = float(xv)
                        _status(
                            f"[auto_fer_eval][adaptive] trigger(from cache): snr={snr_trigger} fer_eff={fer_eff:g} <= {trigger_fer:g}"
                        )
                else:
                    fcw = cached.get("fail_cw")
                    _status(
                        f"[auto_fer_eval][adaptive] skip trigger(from cache): snr={xv} FAIL_CW={fcw} < {min_fail}"
                    )
            return

        # Fallback: try to parse log file if manifest has no record.
        rec_done = _parse_done_log(c, xv, log_path, cmd)
        if rec_done is not None:
            upsert_run(manifest, rec_done)
            ok_count += 1
            fer_eff = fer_effective(rec_done.ldpc_fer, rec_done.total_packets)
            # Print cached result
            fcw_s = str(rec_done.fail_cw) if rec_done.fail_cw is not None else "?"
            tp_s = str(rec_done.total_packets) if rec_done.total_packets is not None else "?"
            fer_s = f"{rec_done.ldpc_fer:.6g}" if rec_done.ldpc_fer is not None else "?"
            _status(f"[auto_fer_eval][adaptive] cached(log): snr={xv} FER={fer_s} FAIL_CW={fcw_s} N={tp_s}")
            if fer_eff is not None and fer_eff <= trigger_fer:
                if _is_reliable_for_trigger(rec_done.fail_cw, fer_eff):
                    if snr_trigger is None or xv < snr_trigger:
                        snr_trigger = float(xv)
                        _status(
                            f"[auto_fer_eval][adaptive] trigger(from log): snr={snr_trigger} fer_eff={fer_eff:g} <= {trigger_fer:g}"
                        )
                else:
                    fcw = rec_done.fail_cw
                    _status(
                        f"[auto_fer_eval][adaptive] skip trigger(from log): snr={xv} FAIL_CW={fcw} < {min_fail}"
                    )
            return

        job_name = f"{Path(c.name).name}_{c.log_prefix}_{_format_axis(xv, 'snr')}"
        job = executor.submit(cmd, cwd=c.workdir, log_path=log_path, job_name=job_name, queue=queue)
        backend, jid = _job_identity(job)
        job_key = f"{backend}:{jid}"
        t_submit = time.time()
        job_db.upsert(
            {
                "key": job_key,
                "backend": backend,
                "job_id": jid,
                "job_name": job_name,
                "stage": "main",
                "axis_type": "snr",
                "axis_value": float(xv),
                "log_path": str(log_path),
                "cmd": cmd,
                "submitted_at": t_submit,
            }
        )
        job_db.set_state(job_key, "SUBMITTED", t=t_submit)
        job_db.flush()
        in_flight[xv] = {
            "job": job,
            "job_key": job_key,
            "submitted_at": time.time(),
            "running_since": None,
            "log_path": log_path,
            "cmd": cmd,
        }
        _status(f"[auto_fer_eval][adaptive] submit snr={xv} -> {log_path.name}")

    def cancel_all_in_flight(*, reason: str) -> None:
        for xv, info in list(in_flight.items()):
            job = info["job"]
            job_key = str(info.get("job_key") or "")
            try:
                executor.cancel(job)
            except Exception:
                pass
            if job_key:
                job_db.mark_cancel(job_key, reason=reason)
            in_flight.pop(xv, None)
        job_db.flush()

    def current_capacity() -> int:
        # Circuit breaker warmup: only allow full concurrency after we have parsed
        # at least N successful points. This prevents mass submission when the
        # very first job is misconfigured (e.g., wrong matrix_dir).
        if ok_count < startup_ok_required:
            return 1
        return max_in_flight

    # Submit/poll loop.
    while True:
        # Fill the pipeline.
        while len(in_flight) < current_capacity() and _in_range(x_next, x_stop_coarse, +1):
            if snr_trigger is not None and x_next > snr_trigger + kill_margin + 1e-12:
                break
            submit_one(x_next)
            x_next = _quantize_axis(x_next + step_cur, step_cur, axis_type="snr", origin=origin)

        # Check if no more work.
        if not in_flight and (snr_trigger is not None or not _in_range(x_next, x_stop_coarse, +1)):
            break

        # Poll all in-flight jobs.
        now = time.time()
        done_any = False
        for xv, info in list(in_flight.items()):
            # The in-flight set may be modified mid-loop (e.g., trigger kill or steep backoff).
            if xv not in in_flight:
                continue
            job = info["job"]
            job_key = str(info.get("job_key") or "")
            st = executor.poll(job)
            if job_key:
                job_db.set_state(job_key, st.state, t=now)

            # Track RUN start for timeout accounting (LSF: do not count PEND time).
            if st.state == "RUN" and info.get("running_since") is None:
                info["running_since"] = now
                if job_key:
                    job_db.upsert({"key": job_key, "running_since": now})
                    job_db.flush()

            # Timeout -> kill and parse best-effort partial log, then reduce step.
            if (
                not st.done
                and max_job_runtime_sec > 0
                and st.state == "RUN"
                and info.get("running_since") is not None
                and (now - float(info["running_since"])) > max_job_runtime_sec
            ):
                log_path = Path(info["log_path"])
                cmd = list(info["cmd"])
                _status(
                    f"[auto_fer_eval][adaptive] TIMEOUT: snr={xv} runtime>{max_job_runtime_sec:.0f}s -> kill and parse {log_path.name}"
                )
                try:
                    executor.cancel(job)
                except Exception:
                    pass
                if job_key:
                    job_db.mark_cancel(job_key, reason="timeout")

                # Cancel higher-SNR jobs: if this point is too expensive, higher SNR is unlikely cheaper.
                # Also try to parse partial results and delete invalid logs.
                for xv2, info2 in list(in_flight.items()):
                    if xv2 > xv + 1e-12:
                        try:
                            executor.cancel(info2["job"])
                        except Exception:
                            pass
                        job_key2 = str(info2.get("job_key") or "")
                        if job_key2:
                            job_db.mark_cancel(job_key2, reason="timeout_higher_snr")

                        # Try to parse partial results from killed job
                        log_path2 = Path(info2.get("log_path") or _log_path_for(c, xv2))
                        cmd2 = list(info2.get("cmd") or build_cmd(c, xv2))
                        partial2 = _parse_killed_job_log(log_path2, timeout_grace_sec=2.0)
                        if partial2 is not None and partial2.get("ldpc_fer") is not None:
                            rec2 = RunRecord(
                                axis_type=c.axis_type,
                                axis_value=float(xv2),
                                log_path=str(log_path2),
                                cmd=cmd2,
                                raw_ber=partial2.get("raw_ber"),
                                ldpc_fer=partial2.get("ldpc_fer"),
                                fail_cw=partial2.get("fail_cw"),
                                total_packets=partial2.get("total_packets"),
                            )
                            upsert_run(manifest, rec2)
                            fcw_s2 = str(partial2.get("fail_cw")) if partial2.get("fail_cw") is not None else "?"
                            _status(
                                f"[auto_fer_eval][adaptive] kill(timeout_higher) snr={xv2} "
                                f"partial: FER={partial2.get('ldpc_fer'):.2g} FAIL_CW={fcw_s2}"
                            )
                        else:
                            # No valid partial result - delete the invalid log file
                            if log_path2.exists():
                                try:
                                    log_path2.unlink()
                                    _status(f"[auto_fer_eval][adaptive] kill(timeout_higher) snr={xv2} -> deleted invalid log")
                                except Exception:
                                    _status(f"[auto_fer_eval][adaptive] kill(timeout_higher) snr={xv2} -> log invalid")
                            else:
                                _status(f"[auto_fer_eval][adaptive] kill(timeout_higher) snr={xv2}")

                        in_flight.pop(xv2, None)
                job_db.flush()

                # Best-effort: wait for log to flush and parse the last complete [SIM] block.
                # Use longer grace period for LSF (network filesystem latency).
                grace_sec = max(timeout_log_grace_sec, 30.0)  # At least 30s for LSF
                rec_done = None
                t_end = time.time() + grace_sec
                while time.time() < t_end:
                    rec_done = _parse_done_log(c, xv, log_path, cmd)
                    if rec_done is not None:
                        break
                    time.sleep(1.0)  # Longer sleep for LSF

                if rec_done is None:
                    # Log not parsable yet - don't fail, just warn and continue
                    # The log may be parsable on next run (LSF filesystem latency)
                    _status(
                        f"[auto_fer_eval][adaptive] WARN: timeout log not parsable yet (missing LDPC FER): {log_path}"
                    )
                    _status(
                        f"[auto_fer_eval][adaptive] WARN: snr={xv} will be skipped this run, may be parsable on restart"
                    )
                    # Don't raise error - just continue without this point
                    # Set trigger to this point to stop further submissions
                    if snr_trigger is None or xv < snr_trigger:
                        snr_trigger = float(xv)
                        _status(f"[auto_fer_eval][adaptive] setting snr_trigger={snr_trigger} due to timeout")
                else:
                    upsert_run(manifest, rec_done)
                    ok_count += 1
                    fer_eff = fer_effective(rec_done.ldpc_fer, rec_done.total_packets)
                    fcw = rec_done.fail_cw
                    tp = rec_done.total_packets
                    fcw_s = str(int(fcw)) if fcw is not None else "?"
                    tp_s = str(int(tp)) if tp is not None else "?"
                    _status(
                        f"[auto_fer_eval][adaptive] timeout parsed: snr={xv} FAIL_CW={fcw_s} N={tp_s} FER={rec_done.ldpc_fer:g}"
                    )
                    if fer_eff is not None and fer_eff <= trigger_fer:
                        if _is_reliable_for_trigger(rec_done.fail_cw, fer_eff):
                            if snr_trigger is None or xv < snr_trigger:
                                snr_trigger = float(xv)
                                _status(
                                    f"[auto_fer_eval][adaptive] trigger(from timeout): snr={snr_trigger} fer_eff={fer_eff:g} <= {trigger_fer:g}"
                                )
                        else:
                            _status(
                                f"[auto_fer_eval][adaptive] skip trigger(from timeout): snr={xv} FAIL_CW={rec_done.fail_cw} < {min_fail}"
                            )
                    # Coarse-only watchdog: if we time out at some SNR, do not schedule
                    # even higher SNR coarse points (they are unlikely cheaper).
                    x_stop_coarse = min(x_stop_coarse, float(xv))

                in_flight.pop(xv, None)
                done_any = True
                continue

            if not st.done:
                continue

            # DRY_RUN_DONE: skip log parsing, just mark as done
            if st.state == "DRY_RUN_DONE":
                in_flight.pop(xv, None)
                done_any = True
                continue

            # DONE/EXIT: try parse log
            cmd = list(info.get("cmd") or build_cmd(c, xv))
            log_path = Path(info.get("log_path") or _log_path_for(c, xv))
            if st.state == "EXIT" and fail_fast:
                msg = f"[auto_fer_eval][adaptive] ERROR: job exited abnormally: snr={xv} log={log_path}"
                _status(msg)
                cancel_all_in_flight(reason="job_exit")
                raise RuntimeError(msg)
            # DONE/EXIT log flush grace:
            # On some LSF setups, bjobs may report DONE before the -o log file
            # is fully flushed to the shared filesystem. Use a bounded grace
            # window (configurable) before declaring the log "not parsable".
            rec_done = None
            is_lsf = hasattr(job, "job_id")
            grace_sec = float(timeout_log_grace_sec)
            if is_lsf:
                grace_sec = max(grace_sec, 30.0)
            else:
                grace_sec = max(grace_sec, 1.0)

            t_end = time.time() + max(0.0, grace_sec)
            sleep_s = 0.2
            while True:
                rec_done = _parse_done_log(c, xv, log_path, cmd)
                if rec_done is not None:
                    break
                if time.time() >= t_end:
                    break
                time.sleep(sleep_s)
                sleep_s = min(2.0, sleep_s * 1.5)

            if rec_done is None:
                msg = f"[auto_fer_eval][adaptive] ERROR: job done but log not parsable (missing LDPC FER): {log_path}"
                _status(msg)
                if fail_fast:
                    cancel_all_in_flight(reason="log_not_parsable")
                    raise RuntimeError(msg)
                done_any = True
            else:
                upsert_run(manifest, rec_done)
                ok_count += 1
                fer_eff = fer_effective(rec_done.ldpc_fer, rec_done.total_packets)
                # Print completed job info immediately
                fcw_s = str(rec_done.fail_cw) if rec_done.fail_cw is not None else "?"
                tp_s = str(rec_done.total_packets) if rec_done.total_packets is not None else "?"
                fer_s = f"{rec_done.ldpc_fer:.6g}" if rec_done.ldpc_fer is not None else "?"
                _status(f"[auto_fer_eval][adaptive] done: snr={xv} FER={fer_s} FAIL_CW={fcw_s} N={tp_s}")
                if fer_eff is not None and fer_eff <= trigger_fer:
                    if _is_reliable_for_trigger(rec_done.fail_cw, fer_eff):
                        if snr_trigger is None or xv < snr_trigger:
                            snr_trigger = float(xv)
                            _status(
                                f"[auto_fer_eval][adaptive] trigger: snr={snr_trigger} fer_eff={fer_eff:g} <= {trigger_fer:g}"
                            )
                    else:
                        _status(
                            f"[auto_fer_eval][adaptive] skip trigger: snr={xv} FAIL_CW={rec_done.fail_cw} < {min_fail}"
                        )
                done_any = True
            if job_key:
                job_db.upsert({"key": job_key, "finished_at": now})
                job_db.flush()
            in_flight.pop(xv, None)

        # Triggered: kill higher-SNR jobs immediately (no margin by default).
        # Also try to parse partial results from killed jobs.
        if snr_trigger is not None:
            for xv, info in list(in_flight.items()):
                if xv > snr_trigger + kill_margin + 1e-12:
                    executor.cancel(info["job"])
                    job_key2 = str(info.get("job_key") or "")
                    if job_key2:
                        job_db.mark_cancel(job_key2, reason="trigger_kill")

                    # Try to parse partial results from killed job
                    log_path = Path(info.get("log_path") or _log_path_for(c, xv))
                    cmd = list(info.get("cmd") or build_cmd(c, xv))
                    partial = _parse_killed_job_log(log_path, timeout_grace_sec=2.0)
                    if partial is not None and partial.get("ldpc_fer") is not None:
                        rec = RunRecord(
                            axis_type=c.axis_type,
                            axis_value=float(xv),
                            log_path=str(log_path),
                            cmd=cmd,
                            raw_ber=partial.get("raw_ber"),
                            ldpc_fer=partial.get("ldpc_fer"),
                            fail_cw=partial.get("fail_cw"),
                            total_packets=partial.get("total_packets"),
                        )
                        upsert_run(manifest, rec)
                        fcw_s = str(partial.get("fail_cw")) if partial.get("fail_cw") is not None else "?"
                        tp_s = str(partial.get("total_packets")) if partial.get("total_packets") is not None else "?"
                        _status(
                            f"[auto_fer_eval][adaptive] kill snr={xv} (>{snr_trigger}+{kill_margin}) "
                            f"partial: FER={partial.get('ldpc_fer'):.2g} FAIL_CW={fcw_s} N={tp_s}"
                        )
                    else:
                        # No valid partial result - delete the invalid log file
                        if log_path.exists():
                            try:
                                log_path.unlink()
                                _status(f"[auto_fer_eval][adaptive] kill snr={xv} (>{snr_trigger}+{kill_margin}) -> deleted invalid log")
                            except Exception:
                                _status(f"[auto_fer_eval][adaptive] kill snr={xv} (>{snr_trigger}+{kill_margin}) -> log invalid")
                        else:
                            _status(f"[auto_fer_eval][adaptive] kill snr={xv} (>{snr_trigger}+{kill_margin})")

                    in_flight.pop(xv, None)
            job_db.flush()

        if not done_any:
            time.sleep(poll_sec)

    if snr_trigger is not None:
        _status(f"[auto_fer_eval][adaptive] main scan done. snr_trigger={snr_trigger}")
    else:
        _status("[auto_fer_eval][adaptive] main scan done. snr_trigger=<none>")
    return snr_trigger


def adaptive_deep_scan(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    executor: Any,
    start_snr: float,
    stop_snr: float,
    poll_sec: float,
    fail_fast: bool,
    min_fail_cw_for_decision: int,
    max_job_runtime_sec: float,
    timeout_log_grace_sec: float,
    max_in_flight: int = 10,
    queue: str = "",
    export_xlsx_sec: float = 0.0,
) -> None:
    """
    Parallel deep scan from start_snr to stop_snr.

    Instead of sequential step-by-step scanning, we:
    1. Generate all points from start_snr to stop_snr using step_low
    2. Submit them in parallel
    """
    if c.axis_type != "snr":
        return

    poll_sec = max(0.2, float(poll_sec))
    max_in_flight = max(1, int(max_in_flight))
    min_fail = max(0, int(min_fail_cw_for_decision))
    job_db_path = (Path(c.out_dir) / c.name / "jobs.json").resolve()
    job_db = JobDB(job_db_path)

    export_sec = float(export_xlsx_sec) if export_xlsx_sec is not None else 0.0
    xlsx_path: Optional[Path] = None
    next_export_at: Optional[float] = None
    if export_sec > 0:
        # Place xlsx next to "adaptive/" (i.e., parent of "main/") so each case has one file.
        xlsx_path = (Path(c.out_dir) / c.name).resolve().parent / f"{c.log_prefix}_progress.xlsx"
        next_export_at = time.time()

    def _job_identity(job: Any) -> tuple[str, str]:
        if hasattr(job, "job_id"):
            return ("lsf", str(getattr(job, "job_id")))
        if hasattr(job, "proc") and hasattr(getattr(job, "proc"), "pid"):
            return ("local", str(getattr(getattr(job, "proc"), "pid")))
        return ("unknown", str(id(job)))

    # Generate all deep scan points.
    deep_points: list[float] = []
    x = float(start_snr)
    origin = x
    step = c.step_policy.step_low

    while _in_range(x, stop_snr, +1):
        x = _quantize_axis(x, step, axis_type="snr", origin=origin)
        if not _in_range(x, stop_snr, +1):
            break
        # Skip if already exists with sufficient fail_cw.
        cached = find_run(manifest, "snr", x)
        if cached and cached.get("ldpc_fer") is not None:
            # Check if fail_cw is sufficient for decision
            fcw = cached.get("fail_cw")
            if fcw is not None:
                try:
                    if int(fcw) >= min_fail:
                        # Valid cached result - skip
                        x = _quantize_axis(x + step, step, axis_type="snr", origin=origin)
                        continue
                except (TypeError, ValueError):
                    pass
            # fail_cw is None or insufficient - need to re-run
            _status(f"[auto_fer_eval][adaptive] deep scan: snr={x} has insufficient fail_cw={fcw}, will re-run")
        deep_points.append(x)
        x = _quantize_axis(x + step, step, axis_type="snr", origin=origin)

    if not deep_points:
        _status(f"[auto_fer_eval][adaptive] deep scan: no new points needed (start={start_snr} stop={stop_snr})")
        if xlsx_path is not None:
            export_progress_xlsx(
                c,
                manifest,
                xlsx_path=xlsx_path,
                job_db_path=job_db_path,
                planned_snrs=[],
                pending_snrs=[],
                in_flight_snrs=[],
            )
        return

    _status(f"[auto_fer_eval][adaptive] deep scan: start={start_snr} stop={stop_snr} ({len(deep_points)} points)")

    # Submit and poll in parallel.
    in_flight: dict[float, dict[str, Any]] = {}
    pending = list(deep_points)
    completed = 0

    # Initial snapshot export once deep scan starts.
    if xlsx_path is not None and next_export_at is not None and time.time() >= next_export_at:
        export_progress_xlsx(
            c,
            manifest,
            xlsx_path=xlsx_path,
            job_db_path=job_db_path,
            planned_snrs=deep_points,
            pending_snrs=pending,
            in_flight_snrs=list(in_flight.keys()),
        )
        next_export_at = time.time() + export_sec

    while pending or in_flight:
        # Periodic export (even if logs are incomplete).
        if xlsx_path is not None and next_export_at is not None and time.time() >= next_export_at:
            export_progress_xlsx(
                c,
                manifest,
                xlsx_path=xlsx_path,
                job_db_path=job_db_path,
                planned_snrs=deep_points,
                pending_snrs=pending,
                in_flight_snrs=list(in_flight.keys()),
            )
            next_export_at = time.time() + export_sec

        # Fill pipeline.
        while pending and len(in_flight) < max_in_flight:
            xv = pending.pop(0)
            log_path = _log_path_for(c, xv)
            cmd = build_cmd(c, xv)
            job_name = f"{Path(c.name).name}_{c.log_prefix}_{_format_axis(xv, 'snr')}"
            job = executor.submit(cmd, cwd=c.workdir, log_path=log_path, job_name=job_name, queue=queue)
            backend, jid = _job_identity(job)
            job_key = f"{backend}:{jid}"
            t_submit = time.time()
            job_db.upsert(
                {
                    "key": job_key,
                    "backend": backend,
                    "job_id": jid,
                    "job_name": job_name,
                    "stage": "deep",
                    "axis_type": "snr",
                    "axis_value": float(xv),
                    "log_path": str(log_path),
                    "cmd": cmd,
                    "submitted_at": t_submit,
                }
            )
            job_db.set_state(job_key, "SUBMITTED", t=t_submit)
            job_db.flush()
            in_flight[xv] = {"job": job, "job_key": job_key, "log_path": log_path, "cmd": cmd}
            _status(f"[auto_fer_eval][adaptive] submit(deep) snr={xv} -> {log_path.name}")

        if not in_flight:
            break

        # Poll all in-flight jobs.
        done_any = False
        for xv, info in list(in_flight.items()):
            job = info["job"]
            job_key = str(info["job_key"])
            st = executor.poll(job)
            job_db.set_state(job_key, st.state, t=time.time())

            if not st.done:
                continue

            # DRY_RUN_DONE: skip log parsing, just mark as done
            if st.state == "DRY_RUN_DONE":
                in_flight.pop(xv, None)
                done_any = True
                continue

            log_path = Path(info["log_path"])
            cmd = list(info["cmd"])

            if st.state == "EXIT" and fail_fast:
                raise RuntimeError(f"[auto_fer_eval][adaptive] ERROR: deep job exited abnormally: snr={xv} log={log_path}")

            rec_done = _parse_done_log(c, xv, log_path, cmd)
            if rec_done is None:
                _status(f"[auto_fer_eval][adaptive] warn: deep log incomplete: {log_path}")
            else:
                if min_fail > 0 and rec_done.fail_cw is not None and int(rec_done.fail_cw) < min_fail:
                    _status(
                        f"[auto_fer_eval][adaptive] warn: deep log has too few FAIL CW for decision: "
                        f"snr={xv} FAIL_CW={rec_done.fail_cw} < {min_fail}"
                    )
                upsert_run(manifest, rec_done)
                completed += 1
                fcw_s = str(rec_done.fail_cw) if rec_done.fail_cw is not None else "?"
                tp_s = str(rec_done.total_packets) if rec_done.total_packets is not None else "?"
                fer_s = f"{rec_done.ldpc_fer:.6g}" if rec_done.ldpc_fer is not None else "?"
                _status(f"[auto_fer_eval][adaptive] done(deep): snr={xv} FER={fer_s} FAIL_CW={fcw_s} N={tp_s}")

                # Check if we've reached fer_lo.
                fer_eff = fer_effective(rec_done.ldpc_fer, rec_done.total_packets)
                if fer_eff is not None and fer_eff <= c.fer_lo:
                    _status(f"[auto_fer_eval][adaptive] reached fer_lo={c.fer_lo:g}, stop deep scan.")
                    # Cancel remaining jobs.
                    for xv2, info2 in list(in_flight.items()):
                        if xv2 != xv:
                            try:
                                executor.cancel(info2["job"])
                            except Exception:
                                pass
                            job_db.mark_cancel(str(info2["job_key"]), reason="fer_lo_reached")
                    job_db.flush()
                    if xlsx_path is not None:
                        export_progress_xlsx(
                            c,
                            manifest,
                            xlsx_path=xlsx_path,
                            job_db_path=job_db_path,
                            planned_snrs=deep_points,
                            pending_snrs=pending,
                            in_flight_snrs=list(in_flight.keys()),
                        )
                    return

            job_db.upsert({"key": job_key, "finished_at": time.time()})
            job_db.flush()
            in_flight.pop(xv, None)
            done_any = True

        if not done_any:
            time.sleep(poll_sec)

    # Final export after deep scan completes normally.
    if xlsx_path is not None:
        export_progress_xlsx(
            c,
            manifest,
            xlsx_path=xlsx_path,
            job_db_path=job_db_path,
            planned_snrs=deep_points,
            pending_snrs=[],
            in_flight_snrs=[],
        )


def _count_fit_points(c: CaseConfig, manifest: dict[str, Any], min_fail_cw: int) -> int:
    """Count how many points fall within the fit FER range."""
    count = 0
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        fer = r.get("ldpc_fer")
        fcw = r.get("fail_cw")
        tp = r.get("total_packets")
        if fer is None:
            continue
        try:
            fer_f = float(fer)
            fcw_i = int(fcw) if fcw is not None else None
            tp_i = int(tp) if tp is not None else None
        except (TypeError, ValueError):
            continue
        if min_fail_cw > 0 and fcw_i is not None and fcw_i < min_fail_cw:
            continue
        fer_eff = fer_effective(fer_f, tp_i)
        if fer_eff is None:
            continue
        if c.fit_fer_lo <= fer_eff <= c.fit_fer_hi:
            count += 1
    return count


def _find_lowest_fer_snr(manifest: dict[str, Any]) -> Optional[float]:
    """Find the SNR with the lowest FER in manifest."""
    best_snr: Optional[float] = None
    best_fer: Optional[float] = None
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        snr = r.get("axis_value")
        fer = r.get("ldpc_fer")
        tp = r.get("total_packets")
        if snr is None or fer is None:
            continue
        try:
            snr_f = float(snr)
            fer_f = float(fer)
            tp_i = int(tp) if tp is not None else None
        except (TypeError, ValueError):
            continue
        fer_eff = fer_effective(fer_f, tp_i)
        if fer_eff is None:
            continue
        if best_fer is None or fer_eff < best_fer:
            best_fer = fer_eff
            best_snr = snr_f
    return best_snr


def _cleanup_invalid_logs_after_main_scan(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    snr_trigger: Optional[float],
    min_fail_cw: int,
) -> int:
    """
    Clean up invalid/incomplete logs after main scan ends.

    For points beyond snr_trigger that have:
    - No fail_cw (log was killed before any complete [SIM] block)
    - fail_cw < min_fail_cw (too few failures for reliable statistics)

    We delete the log file and remove from manifest, so deep scan can re-run them.

    Also scans the log directory for orphan log files (exist but not in manifest).

    Returns the number of cleaned up points.
    """
    if snr_trigger is None:
        return 0

    cleaned = 0
    runs_to_remove: list[dict[str, Any]] = []

    # Part 1: Clean up invalid entries in manifest
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        snr = r.get("axis_value")
        if snr is None:
            continue
        try:
            snr_f = float(snr)
        except (TypeError, ValueError):
            continue

        # Only clean up points beyond trigger
        if snr_f <= snr_trigger + 1e-12:
            continue

        # Check if this entry is invalid/incomplete
        fcw = r.get("fail_cw")
        fer = r.get("ldpc_fer")

        # Invalid if: no FER, or no fail_cw, or fail_cw < threshold
        is_invalid = False
        if fer is None:
            is_invalid = True
        elif fcw is None:
            is_invalid = True
        elif min_fail_cw > 0:
            try:
                if int(fcw) < min_fail_cw:
                    is_invalid = True
            except (TypeError, ValueError):
                is_invalid = True

        if is_invalid:
            runs_to_remove.append(r)
            # Delete the log file
            log_path = r.get("log_path")
            if log_path:
                p = Path(log_path)
                if p.exists():
                    try:
                        p.unlink()
                        _status(f"[auto_fer_eval][adaptive] cleanup: deleted invalid log {p.name} (snr={snr_f})")
                        cleaned += 1
                    except Exception as e:
                        _status(f"[auto_fer_eval][adaptive] cleanup: failed to delete {p.name}: {e}")

    # Remove from manifest
    for r in runs_to_remove:
        try:
            manifest.get("runs", []).remove(r)
        except ValueError:
            pass

    # Part 2: Scan log directory for orphan files (beyond trigger, not in manifest)
    log_dir = Path(c.out_dir) / c.name
    if log_dir.exists():
        # Get all SNR values in manifest
        manifest_snrs = set()
        for r in manifest.get("runs", []):
            if r.get("axis_type") == "snr" and r.get("axis_value") is not None:
                try:
                    manifest_snrs.add(float(r["axis_value"]))
                except (TypeError, ValueError):
                    pass

        # Find orphan log files
        import re
        for log_file in log_dir.glob(f"{c.log_prefix}_snr*.log"):
            # Extract SNR from filename like "demo_snr3.85.log"
            m = re.search(r"snr([0-9.]+)\.log$", log_file.name)
            if m:
                try:
                    snr_f = float(m.group(1))
                    # Only clean up orphans beyond trigger
                    if snr_f > snr_trigger + 1e-12 and snr_f not in manifest_snrs:
                        try:
                            log_file.unlink()
                            _status(f"[auto_fer_eval][adaptive] cleanup: deleted orphan log {log_file.name} (snr={snr_f})")
                            cleaned += 1
                        except Exception as e:
                            _status(f"[auto_fer_eval][adaptive] cleanup: failed to delete orphan {log_file.name}: {e}")
                except (TypeError, ValueError):
                    pass

    if cleaned > 0:
        _status(f"[auto_fer_eval][adaptive] cleanup: removed {cleaned} invalid/orphan entries")

    return cleaned


def _log_has_complete_statistics(log_path: Path) -> bool:
    """
    A "complete" log means it reached the final [STATISTICS] section.

    We intentionally key off statistics lines (not SIM blocks), so that partial logs
    from timeout/kill can still be used for decisions temporarily, but later we can
    re-run selected points to produce a complete, comparable dataset.
    """
    if not log_path.exists():
        return False
    try:
        txt = log_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    import re

    has_fer = re.search(
        r"^\s*\[STATISTICS\]\s+LDPC\s+FER\s*:\s*[0-9eE+\-\.]+\s*$",
        txt,
        flags=re.M,
    )
    return bool(has_fer)


def adaptive_finalize_fit_logs(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    executor: Any,
    snr_trigger: Optional[float],
    poll_sec: float,
    fail_fast: bool,
    min_fail_cw_for_decision: int,
    timeout_log_grace_sec: float,
    max_in_flight: int = 10,
    queue: str = "",
) -> int:
    """
    Finalize (re-run) points in the fit window so their logs end with [STATISTICS].

    Motivation:
      - During main scan we may parse partial logs (last complete [SIM] block) from timeout/kill.
      - Those partial points can be temporarily used for trigger/fit decisions.
      - But for final reporting we want the key fit-window points to be "complete" (have [STATISTICS]).

    Policy (per your latest clarification):
      - Only consider points with SNR <= snr_trigger (if snr_trigger is None, use c.x_stop).
      - Only consider points whose FER_eff is within [c.fit_fer_lo, c.fit_fer_hi].
      - If the log is missing [STATISTICS], delete and re-run it to completion.

    Returns: number of points successfully re-run and parsed.
    """
    if c.axis_type != "snr":
        return 0

    poll_sec = max(0.2, float(poll_sec))
    max_in_flight = max(1, int(max_in_flight))
    min_fail = max(0, int(min_fail_cw_for_decision))
    snr_limit = float(snr_trigger) if snr_trigger is not None else float(c.x_stop)
    grace_sec = max(float(timeout_log_grace_sec), 30.0)

    # Collect candidate SNR points to finalize.
    cand: dict[float, dict[str, Any]] = {}
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        snr = r.get("axis_value")
        fer = r.get("ldpc_fer")
        tp = r.get("total_packets")
        if snr is None or fer is None:
            continue
        try:
            snr_f = float(snr)
            fer_f = float(fer)
            tp_i = int(tp) if tp is not None else None
        except (TypeError, ValueError):
            continue

        if snr_f > snr_limit + 1e-12:
            continue

        fer_eff = fer_effective(fer_f, tp_i)
        if fer_eff is None:
            continue
        if fer_eff < c.fit_fer_lo or fer_eff > c.fit_fer_hi:
            continue

        log_path = Path(r.get("log_path") or _log_path_for(c, snr_f)).resolve()
        if _log_has_complete_statistics(log_path):
            continue

        cand[snr_f] = {"log_path": log_path, "cmd": build_cmd(c, snr_f)}

    if not cand:
        return 0

    snrs = sorted(cand.keys())
    _status(
        f"[auto_fer_eval][adaptive] finalize: {len(snrs)} fit-window points need completion "
        f"(SNR<= {snr_limit:.3f}, FER_eff in [{c.fit_fer_lo:g},{c.fit_fer_hi:g}])"
    )

    job_db = JobDB((Path(c.out_dir) / c.name / "jobs.json").resolve())

    def _job_identity(job: Any) -> tuple[str, str]:
        backend = job.__class__.__name__.replace("Job", "").lower()
        jid = getattr(job, "job_id", None)
        if jid is None:
            jid = getattr(job, "pid", None)
        return backend, str(jid)

    pending = list(snrs)
    in_flight: dict[float, dict[str, Any]] = {}
    completed = 0

    while pending or in_flight:
        while pending and len(in_flight) < max_in_flight:
            xv = pending.pop(0)
            info = cand[xv]
            log_path = Path(info["log_path"])
            cmd = list(info["cmd"])

            if log_path.exists():
                try:
                    log_path.unlink()
                except Exception:
                    pass

            job_name = f"{Path(c.name).name}_{c.log_prefix}_{_format_axis(xv, 'snr')}_finalize"
            job = executor.submit(cmd, cwd=c.workdir, log_path=log_path, job_name=job_name, queue=queue)
            backend, jid = _job_identity(job)
            job_key = f"{backend}:{jid}"
            t_submit = time.time()
            job_db.upsert(
                {
                    "key": job_key,
                    "backend": backend,
                    "job_id": jid,
                    "job_name": job_name,
                    "stage": "finalize_fit",
                    "axis_type": "snr",
                    "axis_value": float(xv),
                    "log_path": str(log_path),
                    "cmd": cmd,
                    "submitted_at": t_submit,
                }
            )
            job_db.set_state(job_key, "SUBMITTED", t=t_submit)
            job_db.flush()
            in_flight[xv] = {"job": job, "job_key": job_key, "log_path": log_path, "cmd": cmd}
            _status(f"[auto_fer_eval][adaptive] submit(finalize) snr={xv:.3f} -> {log_path.name}")

        if not in_flight:
            break

        done_any = False
        now = time.time()
        for xv, info in list(in_flight.items()):
            job = info["job"]
            job_key = str(info["job_key"])
            st = executor.poll(job)
            job_db.set_state(job_key, st.state, t=now)

            if not st.done:
                continue

            if st.state == "DRY_RUN_DONE":
                in_flight.pop(xv, None)
                done_any = True
                continue

            log_path = Path(info["log_path"])
            cmd = list(info["cmd"])

            if st.state == "EXIT" and fail_fast:
                raise RuntimeError(f"[auto_fer_eval][adaptive] ERROR: finalize job exited abnormally: snr={xv} log={log_path}")

            rec_done = None
            t_end = time.time() + grace_sec
            while time.time() < t_end:
                rec_done = _parse_done_log(c, xv, log_path, cmd)
                if rec_done is not None and _log_has_complete_statistics(log_path):
                    break
                time.sleep(1.0)

            if rec_done is None:
                _status(f"[auto_fer_eval][adaptive] warn: finalize log not parsable: {log_path}")
            elif not _log_has_complete_statistics(log_path):
                _status(f"[auto_fer_eval][adaptive] warn: finalize log still missing [STATISTICS]: {log_path}")
                # Still record the best effort metrics, but mark as not completed.
                upsert_run(manifest, rec_done)
            else:
                if min_fail > 0 and rec_done.fail_cw is not None and int(rec_done.fail_cw) < min_fail:
                    _status(
                        f"[auto_fer_eval][adaptive] warn: finalize complete but low FAIL_CW: "
                        f"snr={xv} FAIL_CW={rec_done.fail_cw} < {min_fail}"
                    )
                upsert_run(manifest, rec_done)
                completed += 1
                fcw_s = str(rec_done.fail_cw) if rec_done.fail_cw is not None else "?"
                tp_s = str(rec_done.total_packets) if rec_done.total_packets is not None else "?"
                fer_s = f"{rec_done.ldpc_fer:.6g}" if rec_done.ldpc_fer is not None else "?"
                _status(f"[auto_fer_eval][adaptive] done(finalize): snr={xv:.3f} FER={fer_s} FAIL_CW={fcw_s} N={tp_s}")

            job_db.upsert({"key": job_key, "finished_at": time.time()})
            job_db.flush()
            in_flight.pop(xv, None)
            done_any = True

        if not done_any:
            time.sleep(poll_sec)

    return completed


def _load_latest_jobs_by_snr(job_db_path: Path) -> dict[float, dict[str, Any]]:
    """
    Read jobs.json and return the latest job record (by submitted_at/last_update)
    for each SNR.
    """
    if not job_db_path.exists():
        return {}
    try:
        data = json.loads(job_db_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    jobs = data.get("jobs", [])
    if not isinstance(jobs, list):
        return {}

    out: dict[float, dict[str, Any]] = {}
    for rec in jobs:
        if not isinstance(rec, dict):
            continue
        if str(rec.get("axis_type") or "") != "snr":
            continue
        xv = rec.get("axis_value")
        if xv is None:
            continue
        try:
            snr = float(xv)
        except (TypeError, ValueError):
            continue

        t = rec.get("submitted_at")
        if t is None:
            t = rec.get("last_update")
        try:
            t_f = float(t) if t is not None else 0.0
        except (TypeError, ValueError):
            t_f = 0.0

        key = float(round(snr, 9))
        prev = out.get(key)
        prev_t = float(prev.get("_t", 0.0)) if isinstance(prev, dict) else -1.0
        if prev is None or t_f >= prev_t:
            rec2 = dict(rec)
            rec2["_t"] = t_f
            out[key] = rec2

    # Strip internal key
    for k in list(out.keys()):
        out[k].pop("_t", None)
    return out


def export_progress_xlsx(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    xlsx_path: Path,
    job_db_path: Path,
    planned_snrs: list[float],
    pending_snrs: list[float],
    in_flight_snrs: list[float],
) -> None:
    """
    Export a live progress .xlsx for one case.

    - Includes planned/pending/in-flight points even if log is incomplete.
    - Missing fields stay empty.
    - Adds IsComplete flag: 1 if log has [STATISTICS] LDPC FER, else 0.
    """
    # Unique SNR set
    snr_set: dict[float, float] = {}
    for snr in planned_snrs + pending_snrs + in_flight_snrs:
        try:
            snr_set[float(round(float(snr), 9))] = float(snr)
        except Exception:
            pass
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        xv = r.get("axis_value")
        if xv is None:
            continue
        try:
            snr = float(xv)
        except (TypeError, ValueError):
            continue
        snr_set[float(round(snr, 9))] = snr

    snrs = sorted(snr_set.values())

    pending_set = {float(round(float(x), 9)) for x in pending_snrs}
    inflight_set = {float(round(float(x), 9)) for x in in_flight_snrs}

    latest_jobs = _load_latest_jobs_by_snr(job_db_path)

    header = [
        "SNR",
        "RAW_BER",
        "LDPC_FER",
        "FAIL_CW",
        "PACKETS",
        "AvgIter",
        "IsComplete",
        "JobState",
        "Stage",
        "LogPath",
    ]
    rows: list[list[Any]] = [header]

    for snr in snrs:
        key = float(round(float(snr), 9))
        log_path = _log_path_for(c, snr)
        rec = find_run(manifest, "snr", snr)

        parsed = None
        if log_path.exists() and log_path.stat().st_size > 0:
            try:
                parsed = parse_log_text(log_path.read_text(encoding="utf-8", errors="replace"))
            except Exception:
                parsed = None

        def _pick(name: str) -> Any:
            # Prefer parsed log values, then manifest record.
            if parsed is not None:
                v = getattr(parsed, name, None)
                if v is not None:
                    return v
            if rec is not None:
                v2 = rec.get(name) if isinstance(rec, dict) else None
                if v2 is not None:
                    return v2
            return None

        raw_ber = _pick("raw_ber")
        ldpc_fer = _pick("ldpc_fer")
        fail_cw = _pick("fail_cw")
        total_packets = _pick("total_packets")
        avg_iter = _pick("retry_dec_avg_iter")

        is_complete = 1 if _log_has_complete_statistics(log_path) else 0

        job_state = ""
        stage = ""
        job = latest_jobs.get(key)
        if job is not None:
            job_state = str(job.get("last_state") or "")
            stage = str(job.get("stage") or "")

        # If not in job DB, infer from scheduling lists.
        if not job_state:
            if key in inflight_set:
                job_state = "RUN"
                stage = stage or "deep"
            elif key in pending_set:
                job_state = "PENDING"
                stage = stage or "deep"

        rows.append(
            [
                float(snr),
                raw_ber,
                ldpc_fer,
                fail_cw,
                total_packets,
                avg_iter,
                is_complete,
                job_state or None,
                stage or None,
                str(log_path),
            ]
        )

    try:
        write_xlsx(xlsx_path, sheet_name="progress", rows=rows)
    except Exception as e:
        _status(f"[auto_fer_eval][adaptive] WARN: failed to write xlsx: {xlsx_path} ({e})")



def _print_fer_summary(manifest: dict[str, Any], *, case_name: str) -> None:
    """Print a summary table of all completed points sorted by SNR."""
    pts: list[tuple[float, float, Optional[float], Optional[int], Optional[int]]] = []
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        snr = r.get("axis_value")
        fer = r.get("ldpc_fer")
        if snr is None or fer is None:
            continue
        try:
            snr_f = float(snr)
            fer_f = float(fer)
            raw_ber = float(r["raw_ber"]) if r.get("raw_ber") is not None else None
            fcw = int(r["fail_cw"]) if r.get("fail_cw") is not None else None
            tp = int(r["total_packets"]) if r.get("total_packets") is not None else None
        except (TypeError, ValueError):
            continue
        pts.append((snr_f, fer_f, raw_ber, fcw, tp))

    if not pts:
        print()
        _status(f"[auto_fer_eval][adaptive] FER summary for {case_name}: no completed points")
        return

    pts.sort(key=lambda t: t[0])

    print()
    _status(f"[auto_fer_eval][adaptive] FER summary for {case_name}:")
    print("  SNR      FER          RAW_BER      FAIL_CW  PACKETS")
    print("  -------  -----------  -----------  -------  -------")
    for snr, fer, raw_ber, fcw, tp in pts:
        fer_s = f"{fer:.6g}" if fer > 0 else "0"
        raw_s = f"{raw_ber:.6g}" if raw_ber is not None else "-"
        fcw_s = str(fcw) if fcw is not None else "-"
        tp_s = str(tp) if tp is not None else "-"
        print(f"  {snr:<7.3f}  {fer_s:<11}  {raw_s:<11}  {fcw_s:<7}  {tp_s}")
    print()


def adaptive_backfill_gaps(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    executor: Any,
    backfill_step: float,
    decade_threshold: float,
    poll_sec: float,
    fail_fast: bool,
    timeout_log_grace_sec: float,
    max_in_flight: int = 10,
    queue: str = "",
) -> int:
    """
    Backfill gaps where adjacent points have FER jump > decade_threshold.
    Now runs in parallel for efficiency.

    Returns the number of points added.
    """
    if c.axis_type != "snr":
        return 0
    if backfill_step <= 0 or decade_threshold <= 0:
        return 0

    poll_sec = max(0.2, float(poll_sec))
    max_in_flight = max(1, int(max_in_flight))
    job_db = JobDB((Path(c.out_dir) / c.name / "jobs.json").resolve())

    def _job_identity(job: Any) -> tuple[str, str]:
        if hasattr(job, "job_id"):
            return ("lsf", str(getattr(job, "job_id")))
        if hasattr(job, "proc") and hasattr(getattr(job, "proc"), "pid"):
            return ("local", str(getattr(getattr(job, "proc"), "pid")))
        return ("unknown", str(id(job)))

    # Collect all completed points with valid FER.
    pts: list[tuple[float, float]] = []
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        snr = r.get("axis_value")
        fer = r.get("ldpc_fer")
        tp = r.get("total_packets")
        if snr is None or fer is None:
            continue
        try:
            snr_f = float(snr)
            fer_f = float(fer)
            tp_i = int(tp) if tp is not None else None
        except (TypeError, ValueError):
            continue
        fer_eff = fer_effective(fer_f, tp_i)
        if fer_eff is None or fer_eff <= 0:
            continue
        pts.append((snr_f, fer_eff))

    if len(pts) < 2:
        return 0

    pts.sort(key=lambda t: t[0])

    # Find gaps where log10(FER) jump > threshold.
    gaps: list[tuple[float, float, float, float]] = []  # (snr_lo, snr_hi, fer_lo, fer_hi)
    for i in range(len(pts) - 1):
        snr_lo, fer_lo = pts[i]
        snr_hi, fer_hi = pts[i + 1]
        if fer_lo <= 0 or fer_hi <= 0:
            continue
        decade_jump = abs(math.log10(fer_hi) - math.log10(fer_lo))
        if decade_jump > decade_threshold:
            gaps.append((snr_lo, snr_hi, fer_lo, fer_hi))

    if not gaps:
        return 0

    # Collect all backfill points to submit.
    backfill_points: list[float] = []
    for snr_lo, snr_hi, fer_lo, fer_hi in gaps:
        decade_jump = abs(math.log10(fer_hi) - math.log10(fer_lo))
        _status(
            f"[auto_fer_eval][adaptive] backfill gap: snr=[{snr_lo:.3f},{snr_hi:.3f}] "
            f"FER=[{fer_lo:.2g},{fer_hi:.2g}] ({decade_jump:.1f} decades)"
        )

        x = snr_lo + backfill_step
        origin = snr_lo
        while x < snr_hi - 1e-12:
            x = _quantize_axis(x, backfill_step, axis_type="snr", origin=origin)
            if x >= snr_hi - 1e-12:
                break
            # Skip if already exists.
            cached = find_run(manifest, "snr", x)
            if cached and cached.get("ldpc_fer") is not None:
                x += backfill_step
                continue
            backfill_points.append(x)
            x += backfill_step

    if not backfill_points:
        return 0

    # Submit and poll in parallel.
    in_flight: dict[float, dict[str, Any]] = {}
    pending = list(backfill_points)
    added = 0

    while pending or in_flight:
        # Fill pipeline.
        while pending and len(in_flight) < max_in_flight:
            xv = pending.pop(0)
            log_path = _log_path_for(c, xv)
            cmd = build_cmd(c, xv)
            job_name = f"{Path(c.name).name}_{c.log_prefix}_{_format_axis(xv, 'snr')}"
            job = executor.submit(cmd, cwd=c.workdir, log_path=log_path, job_name=job_name, queue=queue)
            backend, jid = _job_identity(job)
            job_key = f"{backend}:{jid}"
            t_submit = time.time()
            job_db.upsert(
                {
                    "key": job_key,
                    "backend": backend,
                    "job_id": jid,
                    "job_name": job_name,
                    "stage": "backfill",
                    "axis_type": "snr",
                    "axis_value": float(xv),
                    "log_path": str(log_path),
                    "cmd": cmd,
                    "submitted_at": t_submit,
                }
            )
            job_db.set_state(job_key, "SUBMITTED", t=t_submit)
            job_db.flush()
            in_flight[xv] = {"job": job, "job_key": job_key, "log_path": log_path, "cmd": cmd}
            _status(f"[auto_fer_eval][adaptive] submit(backfill) snr={xv:.3f} -> {log_path.name}")

        if not in_flight:
            break

        # Poll all in-flight jobs.
        done_any = False
        for xv, info in list(in_flight.items()):
            job = info["job"]
            job_key = str(info["job_key"])
            st = executor.poll(job)
            job_db.set_state(job_key, st.state, t=time.time())

            if not st.done:
                continue

            # DRY_RUN_DONE: skip log parsing, just mark as done
            if st.state == "DRY_RUN_DONE":
                in_flight.pop(xv, None)
                done_any = True
                continue

            log_path = Path(info["log_path"])
            cmd = list(info["cmd"])

            if st.state == "EXIT" and fail_fast:
                raise RuntimeError(f"[auto_fer_eval][adaptive] ERROR: backfill job exited abnormally: snr={xv} log={log_path}")

            rec_done = _parse_done_log(c, xv, log_path, cmd)
            if rec_done is None:
                _status(f"[auto_fer_eval][adaptive] warn: backfill log incomplete: {log_path}")
            else:
                upsert_run(manifest, rec_done)
                added += 1
                fcw_s = str(rec_done.fail_cw) if rec_done.fail_cw is not None else "?"
                tp_s = str(rec_done.total_packets) if rec_done.total_packets is not None else "?"
                fer_s = f"{rec_done.ldpc_fer:.6g}" if rec_done.ldpc_fer is not None else "?"
                _status(f"[auto_fer_eval][adaptive] done(backfill): snr={xv} FER={fer_s} FAIL_CW={fcw_s} N={tp_s}")

            job_db.upsert({"key": job_key, "finished_at": time.time()})
            job_db.flush()
            in_flight.pop(xv, None)
            done_any = True

        if not done_any:
            time.sleep(poll_sec)

    return added


def _parse_killed_job_log(log_path: Path, timeout_grace_sec: float) -> Optional[dict[str, Any]]:
    """
    Try to parse a killed job's log file for partial results.
    Returns a dict with parsed metrics or None if not parsable.
    """
    if not log_path.exists():
        return None

    t_end = time.time() + max(0.0, float(timeout_grace_sec))
    while time.time() < t_end:
        try:
            txt = log_path.read_text(encoding="utf-8", errors="replace")
            m = parse_log_text(txt)
            if m.ldpc_fer is not None:
                return {
                    "ldpc_fer": m.ldpc_fer,
                    "raw_ber": m.raw_ber,
                    "fail_cw": m.fail_cw,
                    "total_packets": m.total_packets,
                }
        except Exception:
            pass
        time.sleep(0.5)
    return None


def gate_plan_start_from_prev_zeros(c: CaseConfig, plan: IntervalPlan, prev_case_dir: Optional[Path]) -> IntervalPlan:
    """
    Loop gating: For deeper loops (L1/L2), avoid wasting budget on high-FER region.

    If previous loop exists and contains a *trailing* run of >=K consecutive points with LDPC_FER==0
    near the "good channel" end, then set this loop's start to the beginning of that run.
    """
    if not c.gate_enable:
        return plan
    if prev_case_dir is None:
        return plan
    prev_manifest_path = prev_case_dir / "manifest.json"
    if not prev_manifest_path.exists():
        return plan

    prev_manifest = load_manifest(prev_manifest_path)
    x_gate = _find_prev_zero_tail_start(prev_manifest, axis_type=c.axis_type, direction=c.direction, k=c.gate_zero_k)
    if x_gate is None:
        return plan

    # Ensure gate start is within current plan bounds.
    if c.direction > 0:
        if x_gate < plan.start - 1e-12 or x_gate > plan.stop + 1e-12:
            return plan
    else:
        if x_gate > plan.start + 1e-12 or x_gate < plan.stop - 1e-12:
            return plan

    if abs(x_gate - plan.start) < 1e-12:
        return plan

    _status(
        f"[auto_fer_eval] gate start from prev zeros: {c.axis_type} start {plan.start} -> {x_gate} (K={c.gate_zero_k})"
    )
    return IntervalPlan(start=float(x_gate), stop=float(plan.stop), axis_type=plan.axis_type, direction=plan.direction)


def _find_prev_zero_tail_start(manifest: dict[str, Any], *, axis_type: str, direction: int, k: int) -> Optional[float]:
    if k <= 0:
        return None

    pts: list[tuple[float, float]] = []
    for r in manifest.get("runs", []):
        if r.get("axis_type") != axis_type:
            continue
        x = r.get("axis_value")
        fer = r.get("ldpc_fer")
        if x is None or fer is None:
            continue
        try:
            x_f = float(x)
            fer_f = float(fer)
        except (TypeError, ValueError):
            continue
        pts.append((x_f, fer_f))

    if not pts:
        return None

    # Sort in the scan direction order: bad -> good.
    pts.sort(key=lambda t: t[0] * float(direction))

    # Scan from the good end backward to find the last contiguous run of zeros with len>=k.
    i = len(pts) - 1
    while i >= 0:
        # Skip non-zero points.
        if abs(pts[i][1]) > 0.0:
            i -= 1
            continue

        end = i
        while i >= 0 and abs(pts[i][1]) == 0.0:
            i -= 1
        start = i + 1

        if (end - start + 1) >= k:
            return pts[start][0]

        # Continue searching earlier runs.
    return None


def _linfit(xs: list[float], ys: list[float]) -> Optional[tuple[float, float]]:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    sxx = 0.0
    sxy = 0.0
    for x, y in zip(xs, ys):
        dx = x - mx
        sxx += dx * dx
        sxy += dx * (y - my)
    if not math.isfinite(sxx) or sxx <= 0.0:
        return None
    a = sxy / sxx
    b = my - a * mx
    if not (math.isfinite(a) and math.isfinite(b)):
        return None
    return a, b


def predict_stop_snr_from_manifest(
    c: CaseConfig,
    manifest: dict[str, Any],
    *,
    min_fail_cw_for_decision: int = 0,
) -> Optional[float]:
    """
    Use the approximate log-log linearity of the waterfall region:
      log10(FER) ≈ a*log10(RAW_BER) + b,  for FER in [fit_fer_lo, fit_fer_hi]

    Then extrapolate to fit_target_fer and map back to SNR via:
      log10(RAW_BER) ≈ c*SNR + d  (local linear approximation)

    Returns a suggested SNR stop (already +fit_stop_margin), or None if insufficient data.
    """
    if not c.fit_enable or c.axis_type != "snr":
        return None

    xs: list[float] = []  # log10(raw_ber)
    ys: list[float] = []  # log10(fer_eff)
    snrs: list[float] = []

    min_fail = max(0, int(min_fail_cw_for_decision))
    for r in manifest.get("runs", []):
        if r.get("axis_type") != "snr":
            continue
        snr = r.get("axis_value")
        raw_ber = r.get("raw_ber")
        fer = r.get("ldpc_fer")
        fcw = r.get("fail_cw")
        total_packets = r.get("total_packets")
        if snr is None or raw_ber is None or fer is None:
            continue
        try:
            snr_f = float(snr)
            raw_f = float(raw_ber)
            fer_f = float(fer)
            fcw_i = int(fcw) if fcw is not None else None
            tp_i = int(total_packets) if total_packets is not None else None
        except (TypeError, ValueError):
            continue

        if min_fail > 0 and fcw_i is not None and fcw_i < min_fail:
            continue

        fer_eff = fer_effective(fer_f, tp_i)
        if fer_eff is None or raw_f <= 0.0 or fer_eff <= 0.0:
            continue
        if fer_eff < c.fit_fer_lo or fer_eff > c.fit_fer_hi:
            continue

        xs.append(math.log10(raw_f))
        ys.append(math.log10(fer_eff))
        snrs.append(snr_f)

    if len(xs) < c.fit_min_points:
        return None

    fit1 = _linfit(xs, ys)  # y = a*x + b
    fit2 = _linfit(snrs, xs)  # x = c*SNR + d
    if fit1 is None or fit2 is None:
        return None
    a, b = fit1
    c1, d1 = fit2

    # Expect RAW_BER decreases with SNR => c1 should be negative.
    if c1 >= 0.0 or abs(c1) < 1e-12:
        return None
    if abs(a) < 1e-12:
        return None

    y_target = math.log10(c.fit_target_fer)
    x_target = (y_target - b) / a
    if not math.isfinite(x_target):
        return None

    snr_target = (x_target - d1) / c1
    if not math.isfinite(snr_target):
        return None

    # The target should typically be beyond the max SNR of the fit window.
    if snr_target <= max(snrs) - 1e-6:
        return None

    # Fit-based stop is only used as a planning hint; snap to 0.01 dB to avoid
    # floating noise like 4.5500000037 in prints/manifests.
    return float(round(snr_target + c.fit_stop_margin, 2))


def pilot_find_interval(c: CaseConfig, manifest: dict[str, Any], dry_run: bool, limiter: Optional[PlanLimiter]) -> IntervalPlan:
    """
    Pilot scan to locate [low, high] interval:
      - low: near first point where FER <= fer_hi
      - high: stop when FER <= fer_pilot_stop, or span cap is reached (snr axis)
    """
    x = c.x_start
    found_low_at: Optional[float] = None
    found_hi_at: Optional[float] = None
    cap_stop: Optional[float] = None
    prev_x: Optional[float] = None
    prev_fer: Optional[float] = None

    # Scan from bad channel (low SNR / high K) towards good channel.
    while _in_range(x, c.x_stop, c.direction):
        if cap_stop is not None:
            # Clamp to cap_stop to avoid running past the requested SNR span.
            if c.direction > 0:
                x = min(x, cap_stop)
            else:
                x = max(x, cap_stop)

        rec = run_point(c, manifest, x, dry_run=dry_run, limiter=limiter)
        fer = rec.ldpc_fer

        if found_low_at is None and fer is not None and fer <= c.fer_hi:
            # We just crossed into FER <= fer_hi. If the previous coarse point
            # still had FER > fer_hi, refine within the coarse interval using
            # step_default (e.g., 0.1) to avoid skipping the waterfall.
            found_low_at = x
            if prev_x is not None and prev_fer is not None and prev_fer > c.fer_hi:
                refine_steps = [c.step_policy.step_default, c.step_policy.step_mid, c.step_policy.step_low]
                found_low_at = refine_pilot_start(
                    c,
                    manifest,
                    x_bad=prev_x,
                    x_good=x,
                    fer_hi=c.fer_hi,
                    fer_too_low=c.fer_pilot_too_low,
                    steps=refine_steps,
                    dry_run=dry_run,
                    limiter=limiter,
                )
            if c.axis_type == "snr" and c.snr_span_cap > 0:
                cap_stop = found_low_at + c.direction * c.snr_span_cap
                if c.direction > 0:
                    cap_stop = min(cap_stop, c.x_stop)
                else:
                    cap_stop = max(cap_stop, c.x_stop)

        # High boundary: either FER reaches the pilot stop target, or we hit the span cap.
        if fer is not None:
            fer_eff = fer_effective(fer, rec.total_packets)
            if fer_eff is not None and fer_eff <= c.fer_pilot_stop:
                found_hi_at = x
                break
        if cap_stop is not None and abs(x - cap_stop) < 1e-12:
            found_hi_at = cap_stop
            break

        prev_x = x
        prev_fer = fer
        nxt = _quantize_axis(x + c.direction * c.pilot_step, c.pilot_step, axis_type=c.axis_type, origin=c.x_start)
        if cap_stop is not None:
            if c.direction > 0 and nxt > cap_stop:
                x = cap_stop
            elif c.direction < 0 and nxt < cap_stop:
                x = cap_stop
            else:
                x = nxt
        else:
            x = nxt

    if found_low_at is None:
        if dry_run:
            # Cannot infer interval without running; return full scan range.
            return IntervalPlan(start=c.x_start, stop=c.x_stop, axis_type=c.axis_type, direction=c.direction)
        raise RuntimeError(f"pilot failed: did not find FER <= {c.fer_hi} within range.")

    start = found_low_at - c.direction * c.low_margin
    # Clamp to start boundary (do not go beyond scan start in the "bad" direction).
    if c.direction > 0:
        start = max(c.x_start, start)
    else:
        start = min(c.x_start, start)

    stop = found_hi_at if found_hi_at is not None else c.x_stop
    return IntervalPlan(start=start, stop=stop, axis_type=c.axis_type, direction=c.direction)


def run_point(
    c: CaseConfig, manifest: dict[str, Any], x: float, dry_run: bool, limiter: Optional[PlanLimiter]
) -> RunRecord:
    axis_value = float(x)
    axis_value = _quantize_axis(axis_value, 1e-12, axis_type=c.axis_type, origin=axis_value)

    cached = find_run(manifest, c.axis_type, axis_value)
    if cached:
        log_path = Path(cached["log_path"])
        if log_path.exists():
            m = parse_log_text(log_path.read_text(encoding="utf-8", errors="replace"))
            if m.ldpc_fer is not None:
                rec = RunRecord(
                    axis_type=c.axis_type,
                    axis_value=axis_value,
                    log_path=str(log_path),
                    cmd=[str(x) for x in cached.get("cmd", [])],
                    raw_ber=m.raw_ber,
                    theo_rber=m.theo_rber,
                    ldpc_fer=m.ldpc_fer,
                    fail_cw=m.fail_cw,
                    total_packets=m.total_packets,
                    retry_dec_avg_iter=m.retry_dec_avg_iter,
                )
                upsert_run(manifest, rec)
                return rec

    log_dir = Path(c.out_dir) / c.name
    log_dir.mkdir(parents=True, exist_ok=True)

    axis_str = _format_axis(axis_value, c.axis_type)
    log_path = log_dir / f"{c.log_prefix}_{axis_str}.log"

    cmd = build_cmd(c, axis_value)
    if log_path.exists():
        m = parse_log_text(log_path.read_text(encoding="utf-8", errors="replace"))
        if m.ldpc_fer is not None:
            rec = RunRecord(
                axis_type=c.axis_type,
                axis_value=axis_value,
                log_path=str(log_path),
                cmd=cmd,
                raw_ber=m.raw_ber,
                theo_rber=m.theo_rber,
                ldpc_fer=m.ldpc_fer,
                fail_cw=m.fail_cw,
                total_packets=m.total_packets,
                retry_dec_avg_iter=m.retry_dec_avg_iter,
            )
            upsert_run(manifest, rec)
            return rec

    if dry_run:
        _status(f"[dry-run] {' '.join(cmd)} > {log_path}")
        rec = RunRecord(axis_type=c.axis_type, axis_value=axis_value, log_path=str(log_path), cmd=cmd)
        upsert_run(manifest, rec)
        if limiter is not None and not log_path.exists():
            limiter.on_new_cmd()
        return rec

    _status(f"[run] {c.axis_type}={axis_value} -> {log_path.name}")
    t0 = time.time()
    with log_path.open("w", encoding="utf-8") as f:
        subprocess.run(cmd, cwd=c.workdir, stdout=f, stderr=subprocess.STDOUT, check=False)
    dt = time.time() - t0
    _status(f"[run] done in {dt:.1f}s")

    txt = log_path.read_text(encoding="utf-8", errors="replace")
    m = parse_log_text(txt)

    rec = RunRecord(
        axis_type=c.axis_type,
        axis_value=axis_value,
        log_path=str(log_path),
        cmd=cmd,
        raw_ber=m.raw_ber,
        theo_rber=m.theo_rber,
        ldpc_fer=m.ldpc_fer,
        fail_cw=m.fail_cw,
        total_packets=m.total_packets,
        retry_dec_avg_iter=m.retry_dec_avg_iter,
    )
    upsert_run(manifest, rec)
    return rec


def _record_from_manifest(d: dict[str, Any]) -> RunRecord:
    return RunRecord(
        axis_type=str(d["axis_type"]),
        axis_value=float(d["axis_value"]),
        log_path=str(d["log_path"]),
        cmd=[str(x) for x in d.get("cmd", [])],
        raw_ber=d.get("raw_ber"),
        theo_rber=d.get("theo_rber"),
        ldpc_fer=d.get("ldpc_fer"),
        fail_cw=d.get("fail_cw"),
        total_packets=d.get("total_packets"),
        retry_dec_avg_iter=d.get("retry_dec_avg_iter"),
    )


def save_interval(path: Path, plan: IntervalPlan) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "axis_type": plan.axis_type,
        "direction": int(plan.direction),
        "start": float(plan.start),
        "stop": float(plan.stop),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def load_interval(path: Path) -> IntervalPlan:
    d = json.loads(path.read_text(encoding="utf-8"))
    return IntervalPlan(
        axis_type=str(d["axis_type"]),
        direction=int(d["direction"]),
        start=float(d["start"]),
        stop=float(d["stop"]),
    )


def build_cmd(c: CaseConfig, x: float) -> list[str]:
    cmd: list[str] = [c.exe, c.sim_mode, c.config, c.ch_model]
    if c.ch_model.upper() != "CLEAN":
        # In this project, ch_para is SNR for AWGN, K for ERR_INJ, etc.
        if c.axis_type == "k":
            cmd.append(str(int(round(x))))
        else:
            cmd.append(str(x))
    cmd.extend(c.cmd_extra_args)
    return cmd


def _format_axis(x: float, axis_type: str) -> str:
    if axis_type == "k":
        return f"k{int(round(x))}"
    # snr: keep up to 3 decimals (for 0.025 steps)
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return f"snr{s}"


def _step_decimals(step: float) -> int:
    if not math.isfinite(step) or step <= 0:
        return 12
    s = f"{step:.12g}"
    if "e" in s or "E" in s:
        return 12
    if "." not in s:
        return 0
    return len(s.split(".")[1].rstrip("0"))


def _quantize_axis(x: float, step: float, axis_type: str, *, origin: float = 0.0) -> float:
    if axis_type == "k":
        return float(int(round(x)))
    if step <= 0:
        return x

    # Quantize relative to origin to avoid snapping to a global 0-grid.
    # Example: x_start=3.2, step=0.5 should produce 3.2, 3.7, 4.2, ...
    q = origin + round((x - origin) / step) * step
    q = round(q, _step_decimals(step))
    # Avoid negative zero and float noise.
    if abs(q) < 1e-12:
        q = 0.0
    return float(q)


def _in_range(x: float, stop: float, direction: int) -> bool:
    if direction > 0:
        return x <= stop + 1e-12
    return x >= stop - 1e-12


def fer_effective(fer: Optional[float], total_packets: Optional[int]) -> Optional[float]:
    """
    Conservative FER estimate used for step selection:
    - if fer>0, use fer
    - if fer==0 and total_packets known, use 95% upper bound ~ 3/N
    """
    if fer is None or not math.isfinite(fer):
        return None
    if fer > 0.0:
        return fer
    if total_packets and total_packets > 0:
        return 3.0 / float(total_packets)
    return fer


def choose_step(
    c: CaseConfig,
    *,
    cur_x: float,
    cur_fer_eff: Optional[float],
    prev_x: Optional[float],
    prev_fer_eff: Optional[float],
) -> float:
    """
    Step selection = segmented step policy + slope guard.

    Segmented policy:
      - uses StepPolicy thresholds (fer_mid/fer_low)

    Slope guard:
      - limits |Δlog10(FER)| per SNR step to <= slope_max_decades (e.g., 0.5 decade),
        so we don't jump over the waterfall with too coarse a step.
    """
    step_seg = c.step_policy.step_for(cur_fer_eff)
    if c.axis_type != "snr":
        return step_seg
    if c.slope_max_decades <= 0:
        return step_seg
    if (
        prev_x is None
        or prev_fer_eff is None
        or cur_fer_eff is None
        or prev_fer_eff <= 0
        or cur_fer_eff <= 0
        or cur_x == prev_x
    ):
        return step_seg

    slope = abs(math.log10(cur_fer_eff) - math.log10(prev_fer_eff)) / abs(cur_x - prev_x)  # decade / dB
    if not math.isfinite(slope) or slope <= 0:
        return step_seg

    step_slope = c.slope_max_decades / slope
    step_slope = max(c.step_policy.step_low, min(c.step_policy.step_default, step_slope))
    return min(step_seg, step_slope)


def refine_first_below(
    c: CaseConfig,
    manifest: dict[str, Any],
    x_bad: float,
    x_good: float,
    fer_hi: float,
    step: float,
    dry_run: bool,
    limiter: Optional[PlanLimiter],
) -> tuple[float, float, RunRecord]:
    """
    Given a bracket where x_bad is worse (FER>fer_hi) and x_good is better (FER<=fer_hi),
    line-scan with a finer step to find the first x where FER<=fer_hi.

    Returns:
      (last_bad_x, first_good_x, first_good_record)
    """
    # Ensure we move from x_bad toward x_good along c.direction.
    last_bad_x = x_bad
    x = _quantize_axis(x_bad + c.direction * step, step, axis_type=c.axis_type, origin=x_bad)
    while _in_range(x, x_good, c.direction):
        rec = run_point(c, manifest, x, dry_run=dry_run, limiter=limiter)
        if rec.ldpc_fer is not None and rec.ldpc_fer <= fer_hi:
            return last_bad_x, x, rec
        if rec.ldpc_fer is not None and rec.ldpc_fer > fer_hi:
            last_bad_x = x
        x = _quantize_axis(x + c.direction * step, step, axis_type=c.axis_type, origin=x_bad)

    # Fallback: accept coarse good point.
    rec_good = run_point(c, manifest, x_good, dry_run=dry_run, limiter=limiter)
    return x_bad, x_good, rec_good


def refine_pilot_start(
    c: CaseConfig,
    manifest: dict[str, Any],
    x_bad: float,
    x_good: float,
    fer_hi: float,
    fer_too_low: float,
    steps: list[float],
    dry_run: bool,
    limiter: Optional[PlanLimiter],
) -> float:
    """
    Refine the first point below fer_hi, but avoid "overshooting" too deep into the waterfall.

    Heuristic:
      - First, refine the coarse bracket with step_default (e.g., 0.1)
      - If the found point is still too good (FER < fer_too_low, e.g., 0.1),
        then keep refining inside the smaller bracket with step_mid / step_low.
    """
    cur_bad = x_bad
    cur_good = x_good
    for step in steps:
        if step <= 0:
            continue
        cur_bad, cur_good, rec_good = refine_first_below(
            c, manifest, cur_bad, cur_good, fer_hi, step, dry_run=dry_run, limiter=limiter
        )

        fer_eff = fer_effective(rec_good.ldpc_fer, rec_good.total_packets)
        if fer_eff is None:
            return cur_good
        if fer_eff >= fer_too_low:
            return cur_good

    return cur_good


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
