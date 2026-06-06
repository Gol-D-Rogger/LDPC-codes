#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from typing import Any, Optional


def _ensure_repo_on_syspath() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    p = str(repo_root)
    if p not in sys.path:
        sys.path.insert(0, p)


_ensure_repo_on_syspath()

try:
    from auto_throughput_eval import auto_throughput_eval as base
except ImportError:
    import auto_throughput_eval as base
from auto_fer_eval.executors import LocalExecutor, LsfExecutor
from auto_fer_eval.job_db import JobDB
from auto_fer_eval.log_parser import parse_log_text
from auto_fer_eval.manifest import RunRecord, find_run, load_manifest, save_manifest, upsert_run


CaseConfig = base.CaseConfig
RunnerConfig = base.RunnerConfig
LsfConfig = base.LsfConfig

HIGH_FER_KILL_PACKETS = 1000
HIGH_FER_KILL_THRESHOLD = 0.99


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(
        prog="auto_throughput_eval_realtime",
        description="Fixed-grid throughput runner with realtime progress csv export.",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run throughput grid for all cases")
    run.add_argument("--config", required=True, help="TOML config path")
    run.add_argument("--dry-run", action="store_true", help="Only print the plan; do not execute commands")

    exp = sub.add_parser("export", help="Refresh progress csv and export final throughput csv if complete")
    exp.add_argument("--config", required=True, help="TOML config path")

    exp_live = sub.add_parser("export-progress", help="Only refresh realtime progress csv from existing files")
    exp_live.add_argument("--config", required=True, help="TOML config path")

    ns = ap.parse_args(argv)
    if ns.cmd == "run":
        cases, runner, lsf = base.load_config(Path(ns.config))
        for c in cases:
            if runner.executor == "local":
                run_case_local(c, runner=runner, dry_run=bool(ns.dry_run))
            else:
                run_case_lsf(c, runner=runner, lsf=lsf, dry_run=bool(ns.dry_run))
            if not bool(ns.dry_run):
                base.export_case_csv(c, lsf=lsf)
        return 0

    if ns.cmd == "export-progress":
        cases, _runner, lsf = base.load_config(Path(ns.config))
        for c in cases:
            export_progress_from_disk(c, lsf=lsf)
        return 0

    if ns.cmd == "export":
        cases, _runner, lsf = base.load_config(Path(ns.config))
        rc = 0
        for c in cases:
            export_progress_from_disk(c, lsf=lsf)
            try:
                base.export_case_csv(c, lsf=lsf)
            except SystemExit as e:
                print(str(e))
                rc = 1
        return rc

    raise SystemExit(f"Unknown cmd: {ns.cmd}")


def _progress_csv_path(c: CaseConfig) -> Path:
    return ((Path(c.out_dir) / c.name).resolve() / "throughput_progress.csv").resolve()


def _load_latest_jobs_by_point(job_db_path: Path, *, axis_type: str) -> dict[str, dict[str, Any]]:
    grouped = base._load_job_records_by_point(job_db_path, axis_type=axis_type)
    out: dict[str, dict[str, Any]] = {}
    for point_key, recs in grouped.items():
        if recs:
            out[point_key] = dict(recs[0])
    return out


def _read_text(log_path: Path) -> Optional[str]:
    if not log_path.exists():
        return None
    try:
        return log_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return None


def _metrics_from_text(text: str, *, source: str) -> Optional[dict[str, Any]]:
    try:
        parsed = parse_log_text(text)
    except Exception:
        return None
    complete = bool(base._text_has_final_statistics(text))
    out = {
        "source": str(source),
        "is_complete": 1 if complete else 0,
    }
    any_metric = False
    for name in ("raw_ber", "ldpc_fer", "fail_cw", "total_packets"):
        v = getattr(parsed, name, None)
        if v is not None:
            out[name] = v
            any_metric = True
    if complete and getattr(parsed, "retry_dec_avg_iter", None) is not None:
        out["retry_dec_avg_iter"] = parsed.retry_dec_avg_iter
        any_metric = True
    return out if any_metric else None


def _first_int(line: str) -> Optional[int]:
    m = re.search(r"^\s*([0-9]+)\b", line)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def read_config_max_iteration(cfg_path: Path) -> Optional[int]:
    try:
        lines = cfg_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except Exception:
        return None

    fast_iter = None
    retry_iter = None
    ibex_iter = None
    decoder = ""
    for line in lines:
        lower = line.lower()
        value = _first_int(line)
        if "ldpc decoder" in lower:
            decoder = line.split("//", 1)[0].strip().upper()
        elif "fast decoder maximum iteration number" in lower and value is not None:
            fast_iter = value
        elif "retry decoder maximum iteration number" in lower and value is not None:
            retry_iter = value
        elif "ibex maximum iteration number" in lower and "normal" not in lower and value is not None:
            ibex_iter = value

    # Fallback to the established cnfg line positions if comments were stripped.
    if fast_iter is None and len(lines) > 13:
        fast_iter = _first_int(lines[13])
    if retry_iter is None and len(lines) > 15:
        retry_iter = _first_int(lines[15])
    if ibex_iter is None and len(lines) > 17:
        ibex_iter = _first_int(lines[17])

    if "IBEX" in decoder and ibex_iter is not None:
        return ibex_iter
    if "FDEC" in decoder and fast_iter is not None:
        return fast_iter
    if any(s in decoder for s in ("RDEC", "MIX", "SKIP")) and retry_iter is not None:
        return retry_iter

    for value in (ibex_iter, retry_iter, fast_iter):
        if value is not None:
            return value
    return None


def _high_fer_kill_candidate(metrics: Optional[dict[str, Any]]) -> bool:
    if not metrics:
        return False
    try:
        packets = int(metrics.get("total_packets"))
        fer = float(metrics.get("ldpc_fer"))
    except (TypeError, ValueError):
        return False
    return packets >= HIGH_FER_KILL_PACKETS and fer >= HIGH_FER_KILL_THRESHOLD


def _high_fer_run_record(
    *,
    c: CaseConfig,
    axis_value: float,
    log_path: Path,
    cmd: list[str],
    metrics: dict[str, Any],
    max_iter: int,
) -> RunRecord:
    def _float_or_none(name: str) -> Optional[float]:
        v = metrics.get(name)
        if v is None:
            return None
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    def _int_or_none(name: str) -> Optional[int]:
        v = metrics.get(name)
        if v is None:
            return None
        try:
            return int(v)
        except (TypeError, ValueError):
            return None

    return RunRecord(
        axis_type=str(c.axis_type),
        axis_value=float(axis_value),
        log_path=str(log_path),
        cmd=list(cmd),
        raw_ber=_float_or_none("raw_ber"),
        theo_rber=_float_or_none("theo_rber"),
        ldpc_fer=_float_or_none("ldpc_fer"),
        fail_cw=_int_or_none("fail_cw"),
        total_packets=_int_or_none("total_packets"),
        retry_dec_avg_iter=float(max_iter),
    )


def _append_high_fer_statistics(log_path: Path, rec: RunRecord) -> None:
    existing = _read_text(log_path) or ""
    if base._text_has_final_statistics(existing):
        return

    def _fmt(v: Any) -> str:
        if v is None:
            return ""
        if isinstance(v, int):
            return str(v)
        try:
            return f"{float(v):.15g}"
        except (TypeError, ValueError):
            return str(v)

    lines = [
        "",
        "[auto_throughput_eval_rt] high FER early stop: synthetic throughput result",
        f"[STATISTICS] Total packets simulated: {_fmt(rec.total_packets)}",
    ]
    if rec.raw_ber is not None:
        lines.append(f"[STATISTICS] RAW BER   : {_fmt(rec.raw_ber)}")
    if rec.theo_rber is not None:
        lines.append(f"[STATISTICS] TheoRBER  : {_fmt(rec.theo_rber)}")
    lines.extend(
        [
            f"[STATISTICS] LDPC FER  : {_fmt(rec.ldpc_fer)}",
            f"[STATISTICS] FAIL CW   : {_fmt(rec.fail_cw)}",
            f"[STATISTICS] Retry Decoder average iterations : {_fmt(rec.retry_dec_avg_iter)}",
        ]
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _partial_metrics_from_log(log_path: Path) -> Optional[dict[str, Any]]:
    txt = _read_text(log_path)
    if not txt:
        return None
    return _metrics_from_text(txt, source="log")


def _live_metrics_score(metrics: Optional[dict[str, Any]]) -> tuple[int, int, int]:
    if not metrics:
        return (-1, -1, -1)
    populated = 0
    for name in ("raw_ber", "ldpc_fer", "fail_cw", "total_packets", "retry_dec_avg_iter"):
        if metrics.get(name) is not None:
            populated += 1
    try:
        packets = int(metrics.get("total_packets")) if metrics.get("total_packets") is not None else -1
    except (TypeError, ValueError):
        packets = -1
    complete = 1 if metrics.get("is_complete") else 0
    return (complete, packets, populated)


def _best_metrics(*candidates: Optional[dict[str, Any]]) -> Optional[dict[str, Any]]:
    picked: Optional[dict[str, Any]] = None
    best_score = (-1, -1, -1)
    for cand in candidates:
        score = _live_metrics_score(cand)
        if score > best_score:
            picked = cand
            best_score = score
    return picked


def _best_live_metrics_from_job(executor: Any, job: Any, log_path: Path) -> Optional[dict[str, Any]]:
    from_log = _partial_metrics_from_log(log_path)

    from_peek = None
    peek_text = getattr(executor, "peek_text", None)
    if callable(peek_text):
        try:
            txt = peek_text(job, timeout_sec=5.0)
        except Exception:
            txt = None
        if txt:
            from_peek = _metrics_from_text(txt, source="bpeek")

    return _best_metrics(from_log, from_peek)


def _job_id_from_info(info: Optional[dict[str, Any]]) -> str:
    if not isinstance(info, dict):
        return ""
    jid = info.get("job_id")
    if jid is None:
        return ""
    return str(jid)


def _job_state_from_info(info: Optional[dict[str, Any]]) -> str:
    if not isinstance(info, dict):
        return ""
    return str(info.get("last_state") or "")


def _local_meta_by_point(local_job_meta: Optional[dict[float, dict[str, Any]]], *, axis_type: str) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for xv, rec in (local_job_meta or {}).items():
        point_key = base._throughput_point_key(axis_type, float(xv))
        out[point_key] = dict(rec)
    return out


def export_progress_csv(
    c: CaseConfig,
    *,
    manifest: dict[str, Any],
    lsf: LsfConfig,
    job_db_path: Optional[Path] = None,
    in_flight: Optional[dict[float, dict[str, Any]]] = None,
    local_job_meta: Optional[dict[float, dict[str, Any]]] = None,
    executor: Any = None,
) -> Path:
    path = _progress_csv_path(c)
    latest_jobs = _load_latest_jobs_by_point(job_db_path, axis_type=c.axis_type) if job_db_path is not None else {}
    local_jobs = _local_meta_by_point(local_job_meta, axis_type=c.axis_type)
    current = in_flight or {}

    point_keys: set[str] = set(latest_jobs.keys()) | set(local_jobs.keys())
    for xv in current.keys():
        point_keys.add(base._throughput_point_key(c.axis_type, float(xv)))
    for rec in manifest.get("runs", []):
        if rec.get("axis_type") != c.axis_type:
            continue
        try:
            point_keys.add(base._throughput_point_key(c.axis_type, float(rec.get("axis_value"))))
        except (TypeError, ValueError):
            continue

    rows: list[list[object]] = []
    axis_header = "SNR" if c.axis_type == "snr" else "K"
    headers = [axis_header, "RBER", "LDPC_FER", "Fail_CW", "Packets", "AvgIter", "JobState", "JobID"]

    axis_values: list[float] = []
    for point_key in point_keys:
        parts = point_key.split(":")
        if len(parts) != 3:
            continue
        tag = parts[-1]
        try:
            if c.axis_type == "k":
                axis_values.append(float(int(tag[1:])))
            else:
                axis_values.append(float(tag[3:]))
        except (TypeError, ValueError):
            continue
    axis_values = sorted(set(axis_values))

    for xv in axis_values:
        point_key = base._throughput_point_key(c.axis_type, xv)
        rec = find_run(manifest, c.axis_type, xv)
        current_info = current.get(float(xv))
        latest_job = latest_jobs.get(point_key)
        local_info = local_jobs.get(point_key)

        log_path_raw = None
        if isinstance(current_info, dict) and current_info.get("log_path"):
            log_path_raw = current_info.get("log_path")
        elif isinstance(local_info, dict) and local_info.get("log_path"):
            log_path_raw = local_info.get("log_path")
        elif isinstance(latest_job, dict) and latest_job.get("log_path"):
            log_path_raw = latest_job.get("log_path")
        elif isinstance(rec, dict) and rec.get("log_path"):
            log_path_raw = rec.get("log_path")
        else:
            log_path_raw = str(base.log_path_for_submission(c, xv, lsf=lsf))
        log_path = Path(str(log_path_raw)).resolve()

        current_live = None
        if executor is not None and isinstance(current_info, dict) and current_info.get("job") is not None:
            current_live = _best_live_metrics_from_job(executor, current_info["job"], log_path)

        partial_log = _partial_metrics_from_log(log_path)
        stored_live = None
        if isinstance(local_info, dict) and isinstance(local_info.get("live_metrics"), dict):
            stored_live = dict(local_info["live_metrics"])
        if stored_live is None and isinstance(latest_job, dict) and isinstance(latest_job.get("live_metrics"), dict):
            stored_live = dict(latest_job["live_metrics"])

        raw_ber = None
        ldpc_fer = None
        fail_cw = None
        total_packets = None
        avg_iter = None
        if rec is not None:
            raw_ber = rec.get("raw_ber")
            ldpc_fer = rec.get("ldpc_fer")
            fail_cw = rec.get("fail_cw")
            total_packets = rec.get("total_packets")
            avg_iter = rec.get("retry_dec_avg_iter")
        else:
            picked = _best_metrics(current_live, partial_log, stored_live)
            if picked is not None:
                raw_ber = picked.get("raw_ber")
                ldpc_fer = picked.get("ldpc_fer")
                fail_cw = picked.get("fail_cw")
                total_packets = picked.get("total_packets")
                avg_iter = picked.get("retry_dec_avg_iter")

        job_state = _job_state_from_info(current_info) or _job_state_from_info(local_info) or _job_state_from_info(latest_job)
        if not job_state and rec is not None:
            job_state = "DONE"
        job_id = _job_id_from_info(current_info) or _job_id_from_info(local_info) or _job_id_from_info(latest_job)

        rows.append([float(xv), raw_ber, ldpc_fer, fail_cw, total_packets, avg_iter, job_state or None, job_id or None])

    base.write_csv(path, headers=headers, rows=rows)
    return path


def export_progress_from_disk(c: CaseConfig, *, lsf: LsfConfig) -> Path:
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)
    job_db_path = (case_dir / "jobs.json").resolve()
    return export_progress_csv(c, manifest=manifest, lsf=lsf, job_db_path=job_db_path)


def _refresh_lsf_live_metrics(*, executor: Any, job_db: JobDB, in_flight: dict[float, dict[str, Any]]) -> None:
    dirty = False
    for info in in_flight.values():
        job_key = str(info.get("job_key") or "")
        if not job_key:
            continue
        live = _best_live_metrics_from_job(executor, info.get("job"), Path(info.get("log_path") or "").resolve())
        if live is None:
            continue
        job_db.upsert({"key": job_key, "live_metrics": live})
        dirty = True
    if dirty:
        job_db.flush()


def _refresh_local_live_metrics(*, executor: Any, in_flight: dict[float, dict[str, Any]], local_job_meta: dict[float, dict[str, Any]]) -> None:
    for xv, info in in_flight.items():
        live = _best_live_metrics_from_job(executor, info.get("job"), Path(info.get("log_path") or "").resolve())
        if live is None:
            continue
        meta = local_job_meta.setdefault(float(xv), {})
        meta["live_metrics"] = live


def run_case_local(c: CaseConfig, *, runner: RunnerConfig, dry_run: bool) -> None:
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)

    cfg_abs = base.resolve_cfg_path(c)
    runtime_cfg = (case_dir / "config_runtime.cnfg").resolve()

    xs = list(base.build_axis_grid(c.x_low, c.x_high, c.x_step, axis_type=c.axis_type))
    todo: list[tuple[float, Path, list[str]]] = []
    cached = 0
    for xv in xs:
        log_path = base.log_path_for(c, xv)
        cmd = base.build_cmd(c, xv, config_arg=str(runtime_cfg))
        rec_done = base._resolve_complete_record(
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
    print(f"\n[auto_throughput_eval_rt] case={c.name} axis={c.axis_type} grid={len(xs)} cached={cached} todo={len(todo)} out={case_dir}")
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
    local_job_meta: dict[float, dict[str, Any]] = {}

    def _flush_progress() -> None:
        _refresh_local_live_metrics(executor=ex, in_flight=in_flight, local_job_meta=local_job_meta)
        export_progress_csv(c, manifest=manifest, lsf=LsfConfig(), in_flight=in_flight, local_job_meta=local_job_meta, executor=ex)

    def cancel_all(reason: str) -> None:
        for xv, info in list(in_flight.items()):
            try:
                ex.cancel(info["job"])
            except Exception:
                pass
            local_job_meta.setdefault(float(xv), {})["last_state"] = "KILLED"
            print(f"[auto_throughput_eval_rt] cancel axis={c.axis_type} x={xv} reason={reason}")
            in_flight.pop(xv, None)

    try:
        if todo:
            base.patch_config_max_sim_num(cfg_abs, runtime_cfg, c.max_sim_num, max_error_num=c.max_err_num)
        high_fer_max_iter = read_config_max_iteration(runtime_cfg if runtime_cfg.exists() else cfg_abs)

        submit_attempts: dict[float, int] = {}
        _flush_progress()
        while todo or in_flight:
            submitted_any = False
            while todo and len(in_flight) < max_in_flight:
                xv, log_path, cmd = todo.pop(0)
                if log_path.exists():
                    try:
                        log_path.unlink()
                    except Exception:
                        pass
                job = ex.submit(cmd, cwd=c.workdir, log_path=log_path)
                _backend, jid = base._job_identity(job)
                meta = local_job_meta.setdefault(float(xv), {})
                meta.update({"job_id": jid, "log_path": str(Path(log_path).resolve()), "last_state": "SUBMITTED"})
                in_flight[xv] = {
                    "job": job,
                    "log_path": Path(log_path).resolve(),
                    "cmd": list(cmd),
                    "job_id": jid,
                    "last_state": "SUBMITTED",
                }
                submit_attempts[xv] = int(submit_attempts.get(xv, 0)) + 1
                submitted_any = True
                print(
                    f"[auto_throughput_eval_rt] submit axis={c.axis_type} x={xv} "
                    f"attempt={submit_attempts[xv]} job_id={jid} -> {Path(log_path).name}"
                )

            if submitted_any:
                _flush_progress()

            if not in_flight:
                break

            time.sleep(poll_sec)
            changed = False
            for xv, info in list(in_flight.items()):
                st = ex.poll(info["job"])
                info["last_state"] = st.state
                local_job_meta.setdefault(float(xv), {})["last_state"] = st.state
                changed = True
                live = _best_live_metrics_from_job(ex, info["job"], Path(info["log_path"]).resolve())
                if live is not None:
                    local_job_meta.setdefault(float(xv), {})["live_metrics"] = live
                if not st.done and _high_fer_kill_candidate(live):
                    if high_fer_max_iter is None:
                        print(
                            f"[auto_throughput_eval_rt] WARN: high FER kill skipped axis={c.axis_type} "
                            f"x={xv}: cannot parse max iteration from {runtime_cfg if runtime_cfg.exists() else cfg_abs}"
                        )
                    else:
                        try:
                            ex.cancel(info["job"])
                        except Exception:
                            pass
                        rec = _high_fer_run_record(
                            c=c,
                            axis_value=xv,
                            log_path=Path(info["log_path"]).resolve(),
                            cmd=list(info["cmd"]),
                            metrics=live,
                            max_iter=high_fer_max_iter,
                        )
                        _append_high_fer_statistics(Path(info["log_path"]).resolve(), rec)
                        upsert_run(manifest, rec)
                        save_manifest(manifest_path, manifest)
                        local_job_meta.setdefault(float(xv), {})["last_state"] = "KILLED_HIGH_FER"
                        in_flight.pop(xv, None)
                        print(
                            f"[auto_throughput_eval_rt] high-fer kill axis={c.axis_type} x={xv} "
                            f"FER={rec.ldpc_fer} packets={rec.total_packets} aver_iter={rec.retry_dec_avg_iter}"
                        )
                        continue
                if not st.done:
                    continue

                in_flight.pop(xv, None)
                rec = base.parse_done_log_with_grace(
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
                    print(f"[auto_throughput_eval_rt] {retry_reason}")
                    if runner.fail_fast:
                        cancel_all(reason=retry_reason)
                        save_manifest(manifest_path, manifest)
                        _flush_progress()
                        raise SystemExit(retry_reason)
                    next_attempt = int(submit_attempts.get(xv, 0)) + 1
                    if next_attempt > runner.max_retries:
                        print(
                            f"[auto_throughput_eval_rt] SKIP axis={c.axis_type} x={xv}: "
                            f"exceeded max_retries={runner.max_retries}"
                        )
                    else:
                        todo.append((float(xv), Path(info["log_path"]).resolve(), list(info["cmd"])))
                        print(
                            f"[auto_throughput_eval_rt] retry axis={c.axis_type} x={xv} "
                            f"next_attempt={next_attempt}"
                        )
                    continue

                upsert_run(manifest, rec)
                save_manifest(manifest_path, manifest)
                local_job_meta.setdefault(float(xv), {})["last_state"] = "DONE"
                print(
                    f"[auto_throughput_eval_rt] done axis={c.axis_type} x={xv} "
                    f"RBER={rec.raw_ber} aver_iter={rec.retry_dec_avg_iter}"
                )

            if changed:
                _flush_progress()
    finally:
        _flush_progress()


def run_case_lsf(c: CaseConfig, *, runner: RunnerConfig, lsf: LsfConfig, dry_run: bool) -> None:
    case_dir = (Path(c.out_dir) / c.name).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = case_dir / "manifest.json"
    manifest = load_manifest(manifest_path)
    job_db_path = (case_dir / "jobs.json").resolve()
    jobs_by_point = base._load_job_records_by_point(job_db_path, axis_type=c.axis_type)

    cfg_abs = base.resolve_cfg_path(c)
    runtime_cfg = (case_dir / "config_runtime.cnfg").resolve()
    matrix_size = base.parse_matrix_size(cfg_abs)
    xs = list(base.build_axis_grid(c.x_low, c.x_high, c.x_step, axis_type=c.axis_type))

    if dry_run:
        todo: list[tuple[float, Path, list[str], str]] = []
        cached = 0
        active_recorded = 0
        for xv in xs:
            log_path = base.log_path_for_submission(c, xv, lsf=lsf)
            cmd = base.build_cmd(c, xv, config_arg=str(runtime_cfg))
            job_name = base.build_job_name(c, matrix_size=matrix_size, axis_value=xv)
            point_key = base._throughput_point_key(c.axis_type, xv)
            rec_done = base._resolve_complete_record(
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
            if any(base._job_record_is_active(rec) for rec in jobs_by_point.get(point_key, [])):
                active_recorded += 1
                continue
            todo.append((xv, log_path, cmd, job_name))

        save_manifest(manifest_path, manifest)
        print(
            f"\n[auto_throughput_eval_rt] case={c.name} axis={c.axis_type} grid={len(xs)} "
            f"cached={cached} active_recorded={active_recorded} todo={len(todo)} out={case_dir}"
        )
        if todo:
            print(
                f"[dry-run] would generate runtime config: {runtime_cfg} (from {cfg_abs}) "
                f"(max_sim_num={c.max_sim_num}, max_err_num={c.max_err_num})"
            )
        cwd_abs = str(Path(c.workdir).resolve())
        for xv, log_path, cmd, job_name in todo:
            bsub_cmd = base.build_bsub_command_line(
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
        log_base_dir="",
    )
    job_db = JobDB(job_db_path)
    max_in_flight = max(1, int(runner.max_in_flight))
    poll_sec = max(0.2, float(runner.poll_sec))

    in_flight: dict[float, dict[str, Any]] = {}
    todo: list[tuple[float, Path, list[str], str, str]] = []
    cached = 0
    adopted = 0

    for xv in xs:
        log_path = base.log_path_for_submission(c, xv, lsf=lsf)
        cmd = base.build_cmd(c, xv, config_arg=str(runtime_cfg))
        job_name = base.build_job_name(c, matrix_size=matrix_size, axis_value=xv)
        point_key = base._throughput_point_key(c.axis_type, xv)

        rec_done = base._resolve_complete_record(
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

            rec_done = base.parse_done_log(axis_type=c.axis_type, axis_value=xv, log_path=job_log_path, cmd=job_cmd, debug=False)
            if rec_done is not None:
                upsert_run(manifest, rec_done)
                cached += 1
                now = time.time()
                if job_key:
                    job_db.set_state(job_key, "DONE", t=now)
                    job_db.upsert({"key": job_key, "finished_at": now})
                reused = True
                break

            if not base._job_record_is_active(job_rec):
                continue

            job = base._resume_executor_job(ex, cmd=job_cmd, cwd=c.workdir, log_path=job_log_path, rec=job_rec)
            if job is None:
                continue

            now = time.time()
            st = ex.poll(job)
            if job_key:
                job_db.set_state(job_key, st.state, t=now)

            if st.done:
                rec_done = base.parse_done_log_with_grace(
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
                "job_id": str(job_rec.get("job_id") or ""),
                "log_path": job_log_path,
                "cmd": job_cmd,
                "job_name": str(job_rec.get("job_name") or job_name),
                "point_key": point_key,
                "last_state": st.state,
            }
            adopted += 1
            reused = True
            print(f"[auto_throughput_eval_rt] adopt axis={c.axis_type} x={xv} job_id={job_rec.get('job_id')} state={st.state}")
            break

        if reused:
            continue

        todo.append((xv, log_path, cmd, job_name, point_key))

    save_manifest(manifest_path, manifest)
    job_db.flush()
    print(
        f"\n[auto_throughput_eval_rt] case={c.name} axis={c.axis_type} grid={len(xs)} "
        f"cached={cached} adopted={adopted} todo={len(todo)} out={case_dir}"
    )

    if todo:
        base.patch_config_max_sim_num(cfg_abs, runtime_cfg, c.max_sim_num, max_error_num=c.max_err_num)

    def _flush_progress() -> None:
        _refresh_lsf_live_metrics(executor=ex, job_db=job_db, in_flight=in_flight)
        export_progress_csv(c, manifest=manifest, lsf=lsf, job_db_path=job_db_path, in_flight=in_flight, executor=ex)

    def cancel_all(reason: str) -> None:
        for xv, info in list(in_flight.items()):
            try:
                ex.cancel(info["job"])
            except Exception:
                pass
            job_key = str(info.get("job_key") or "")
            if job_key:
                job_db.mark_cancel(job_key, reason=reason)
            print(f"[auto_throughput_eval_rt] cancel axis={c.axis_type} x={xv} reason={reason}")
            in_flight.pop(xv, None)
        job_db.flush()

    try:
        submit_attempts: dict[float, int] = {}
        high_fer_max_iter = read_config_max_iteration(runtime_cfg if runtime_cfg.exists() else cfg_abs)
        _flush_progress()
        while todo or in_flight:
            submitted_any = False
            while todo and len(in_flight) < max_in_flight:
                xv, log_path, cmd, job_name, point_key = todo.pop(0)
                if Path(log_path).exists():
                    try:
                        Path(log_path).unlink()
                    except Exception:
                        pass
                job = ex.submit(cmd, cwd=c.workdir, log_path=Path(log_path), job_name=job_name, queue=lsf.queue)
                backend, jid = base._job_identity(job)
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
                    "job_id": jid,
                    "log_path": Path(log_path).resolve(),
                    "cmd": list(cmd),
                    "job_name": job_name,
                    "point_key": point_key,
                    "last_state": "SUBMITTED",
                }
                submit_attempts[xv] = int(submit_attempts.get(xv, 0)) + 1
                submitted_any = True
                print(
                    f"[auto_throughput_eval_rt] bsub axis={c.axis_type} x={xv} "
                    f"attempt={submit_attempts[xv]} job_id={jid} job_name={job_name}"
                )

            if submitted_any:
                _flush_progress()

            if not in_flight:
                break

            time.sleep(poll_sec)
            changed = False
            for xv, info in list(in_flight.items()):
                st = ex.poll(info["job"])
                job_key = str(info.get("job_key") or "")
                now = time.time()
                if job_key:
                    job_db.set_state(job_key, st.state, t=now)
                info["last_state"] = st.state
                changed = True
                live = _best_live_metrics_from_job(ex, info["job"], Path(info["log_path"]).resolve())
                if live is not None and job_key:
                    job_db.upsert({"key": job_key, "live_metrics": live})
                    job_db.flush()
                if not st.done and _high_fer_kill_candidate(live):
                    if high_fer_max_iter is None:
                        print(
                            f"[auto_throughput_eval_rt] WARN: high FER kill skipped axis={c.axis_type} "
                            f"x={xv}: cannot parse max iteration from {runtime_cfg if runtime_cfg.exists() else cfg_abs}"
                        )
                    else:
                        try:
                            ex.cancel(info["job"])
                        except Exception:
                            pass
                        rec = _high_fer_run_record(
                            c=c,
                            axis_value=xv,
                            log_path=Path(info["log_path"]).resolve(),
                            cmd=list(info["cmd"]),
                            metrics=live,
                            max_iter=high_fer_max_iter,
                        )
                        _append_high_fer_statistics(Path(info["log_path"]).resolve(), rec)
                        upsert_run(manifest, rec)
                        save_manifest(manifest_path, manifest)
                        if job_key:
                            job_db.upsert(
                                {
                                    "key": job_key,
                                    "finished_at": now,
                                    "cancel_reason": "high FER >= 0.99 after 1000 packets",
                                    "live_metrics": live,
                                }
                            )
                            job_db.set_state(job_key, "KILLED_HIGH_FER", t=now)
                            job_db.flush()
                        info["last_state"] = "KILLED_HIGH_FER"
                        in_flight.pop(xv, None)
                        print(
                            f"[auto_throughput_eval_rt] high-fer kill axis={c.axis_type} x={xv} "
                            f"FER={rec.ldpc_fer} packets={rec.total_packets} aver_iter={rec.retry_dec_avg_iter}"
                        )
                        continue
                if not st.done:
                    continue

                in_flight.pop(xv, None)
                rec = base.parse_done_log_with_grace(
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
                    print(
                        f"[auto_throughput_eval_rt] done axis={c.axis_type} x={xv} "
                        f"RBER={rec.raw_ber} aver_iter={rec.retry_dec_avg_iter}"
                    )
                    continue

                if job_key:
                    job_db.upsert({"key": job_key, "finished_at": now})
                    job_db.flush()

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

                print(f"[auto_throughput_eval_rt] {retry_reason}")
                if runner.fail_fast:
                    cancel_all(reason=retry_reason)
                    save_manifest(manifest_path, manifest)
                    _flush_progress()
                    raise SystemExit(retry_reason)

                next_attempt = int(submit_attempts.get(xv, 0)) + 1
                if next_attempt > runner.max_retries:
                    print(
                        f"[auto_throughput_eval_rt] SKIP axis={c.axis_type} x={xv}: "
                        f"exceeded max_retries={runner.max_retries}"
                    )
                else:
                    point_key = str(info.get("point_key") or base._throughput_point_key(c.axis_type, float(xv)))
                    todo.append(
                        (
                            float(xv),
                            Path(info["log_path"]).resolve(),
                            list(info["cmd"]),
                            str(info.get("job_name") or base.build_job_name(c, matrix_size=matrix_size, axis_value=float(xv))),
                            point_key,
                        )
                    )
                    print(f"[auto_throughput_eval_rt] retry axis={c.axis_type} x={xv} next_attempt={next_attempt}")

            if changed:
                _flush_progress()
    finally:
        job_db.flush()
        _flush_progress()


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
