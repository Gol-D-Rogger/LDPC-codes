#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]

STAT_LDPC_FER_RE = re.compile(r"^\s*\[STATISTICS\]\s+LDPC\s+FER\s*:\s*([0-9eE+\-\.]+)\s*$", re.M)
STAT_RAW_BER_RE = re.compile(r"^\s*\[STATISTICS\]\s+RAW\s+BER\s*:\s*([0-9eE+\-\.]+)\s*$", re.M)
SIM_LDPC_FER_RE = re.compile(r"^\s*\[SIM\]\s+LDPC\s+FER\s*:\s*([0-9eE+\-\.]+)\s*$", re.M)
SIM_RAW_BER_RE = re.compile(r"^\s*\[SIM\]\s+RAW\s+BER\s*:\s*([0-9eE+\-\.]+)\s*$", re.M)


class SimpleBatchError(RuntimeError):
    pass


@dataclass(frozen=True)
class SizeContext:
    m: int
    n: int
    size_name: str
    output_dir: Path
    config_path: Path
    probe_toml_path: Path
    probe_case_name: str
    probe_out_dir: Path


@dataclass(frozen=True)
class JobSpec:
    matrix_id: int
    snr: float
    job_name: str
    log_path: Path
    cmd: list[str]


def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def log(msg: str) -> None:
    print(f"[{_ts()}] {msg}")


def warn(msg: str) -> None:
    print(f"[{_ts()}] [WARN] {msg}")


def run_cmd(
    cmd: list[str],
    *,
    cwd: Path | None = None,
    capture: bool = True,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        text=True,
        errors="replace",
        check=False,
    )
    if check and p.returncode != 0:
        out = p.stdout or ""
        raise SimpleBatchError(f"Command failed ({p.returncode}): {' '.join(cmd)}\n{out}")
    return p


def parse_m_range(expr: str) -> list[int]:
    vals: set[int] = set()
    for chunk in expr.split(","):
        part = chunk.strip()
        if not part:
            continue
        if "-" in part:
            lo_s, hi_s = part.split("-", 1)
            lo = int(lo_s.strip())
            hi = int(hi_s.strip())
            if lo > hi:
                lo, hi = hi, lo
            vals.update(range(lo, hi + 1))
        else:
            vals.add(int(part))
    if not vals:
        raise SimpleBatchError(f"Invalid --m-range: {expr}")
    return sorted(vals)


def parse_leading_int(line: str, *, field: str) -> int:
    m = re.match(r"^\s*([0-9]+)\b", line)
    if not m:
        raise SimpleBatchError(f"Cannot parse integer for {field}: {line!r}")
    return int(m.group(1))


def patch_leading_int(line: str, value: int) -> str:
    if "//" in line:
        left, right = line.split("//", 1)
        comment = f"//{right.rstrip()}"
    else:
        left, comment = line, ""
    patched, n_sub = re.subn(r"^(\s*)[0-9]+", rf"\g<1>{int(value)}", left, count=1)
    if n_sub == 0:
        raise SimpleBatchError(f"Cannot patch numeric line: {line!r}")
    patched = patched.rstrip()
    if comment:
        return f"{patched} {comment}"
    return patched


def format_snr_value(snr: float) -> str:
    return f"{snr:.6f}".rstrip("0").rstrip(".")


def format_snr_tag(snr: float) -> str:
    return format_snr_value(snr).replace("-", "m").replace(".", "p")


def toml_quote(s: str) -> str:
    return json.dumps(s)


def parse_size_name(size_name: str) -> tuple[int, int] | None:
    m = re.match(r"^([0-9]+)x([0-9]+)$", size_name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def discover_size_name(matrix_root: Path, m: int) -> str | None:
    cands: list[str] = []
    for p in matrix_root.iterdir():
        if not p.is_dir():
            continue
        parsed = parse_size_name(p.name)
        if parsed is None:
            continue
        if parsed[0] == m:
            cands.append(p.name)
    if not cands:
        return None
    cands.sort()
    return cands[0]


def create_size_context(
    *,
    m: int,
    qc: int,
    template_cnfg: Path,
    matrix_root: Path,
    out_root: Path,
) -> SizeContext:
    tpl_lines = template_cnfg.read_text(encoding="utf-8").splitlines()
    if len(tpl_lines) < 5:
        raise SimpleBatchError(f"template too short: {template_cnfg}")

    if qc % 8 != 0:
        raise SimpleBatchError(f"--qc must be divisible by 8, got {qc}")
    row_bytes = qc // 8
    bytes_of_userdata = parse_leading_int(tpl_lines[3], field="bytes_of_userdata")
    n_guess = math.ceil(bytes_of_userdata / row_bytes) + m
    size_name = discover_size_name(matrix_root, m) or f"{m}x{n_guess}"
    parsed = parse_size_name(size_name)
    n_val = parsed[1] if parsed else n_guess

    size_out = (out_root / size_name).resolve()
    for d in (
        size_out / "configs",
        size_out / "toml",
        size_out / "logs",
        size_out / "summary",
        size_out / "summary" / "realtime",
        size_out / "plots",
    ):
        d.mkdir(parents=True, exist_ok=True)

    cnfg_path = size_out / "configs" / f"qc{qc}_M{m}.cnfg"
    tpl_lines[4] = patch_leading_int(tpl_lines[4], row_bytes * m)
    cnfg_path.write_text("\n".join(tpl_lines) + "\n", encoding="utf-8")

    probe_case_name = f"M{m}_probe_id1"
    probe_toml = size_out / "toml" / f"probe_M{m}.toml"
    probe_out = size_out / "probe_autofer"

    log(f"M={m} -> generated cnfg: {cnfg_path}")
    return SizeContext(
        m=m,
        n=n_val,
        size_name=size_name,
        output_dir=size_out,
        config_path=cnfg_path,
        probe_toml_path=probe_toml,
        probe_case_name=probe_case_name,
        probe_out_dir=probe_out,
    )


def write_probe_toml(
    *,
    ctx: SizeContext,
    exe: Path,
    matrix_root: Path,
    queue: str,
    probe_id: int,
    snr_min: float,
    snr_max: float,
    main_step: float,
    pilot_max_sim: int,
    main_max_sim: int,
    probe_max_in_flight: int,
    probe_poll_sec: float,
    probe_timeout_sec: float,
    probe_log_grace_sec: float,
) -> None:
    exe_workdir = str(exe.parent.resolve())
    content = f"""[adaptive]
enable = true
executor = "lsf"
pilot_max_sim_num = {int(pilot_max_sim)}
main_max_sim_num = {int(main_max_sim)}
main_step = {float(main_step)}
trigger_fer = 1e-4
max_in_flight = {int(probe_max_in_flight)}
poll_sec = {float(probe_poll_sec)}
kill_margin = 0.0
startup_ok_required = 1
fail_fast = true
min_fail_cw_for_decision = 1
max_job_runtime_sec = {float(probe_timeout_sec)}
timeout_log_grace_sec = {float(probe_log_grace_sec)}
backfill_step = 0.0

[lsf]
queue = {toml_quote(queue)}
use_cwd = false

[[cases]]
name = {toml_quote(ctx.probe_case_name)}
workdir = {toml_quote(exe_workdir)}
exe = {toml_quote(str(exe))}
sim_mode = "LDPC"
config = {toml_quote(str(ctx.config_path))}
ch_model = "AWGN"
cmd_extra_args = [{toml_quote(str(probe_id))}, {toml_quote(str(matrix_root))}]
out_dir = {toml_quote(str(ctx.probe_out_dir))}
log_prefix = {toml_quote(f"probe_M{ctx.m}_id{probe_id}")}
snr_min = {float(snr_min)}
snr_max = {float(snr_max)}
step_default = {float(main_step)}
step_mid = 0.05
step_low = 0.025
fit_enable = false
"""
    ctx.probe_toml_path.write_text(content, encoding="utf-8")
    log(f"M={ctx.m} -> wrote probe toml: {ctx.probe_toml_path}")


def _probe_main_dir(ctx: SizeContext) -> Path:
    return ctx.probe_out_dir / ctx.probe_case_name / "adaptive" / "main"


def kill_probe_lsf_jobs_best_effort(ctx: SizeContext) -> None:
    jobs_path = _probe_main_dir(ctx) / "jobs.json"
    if not jobs_path.exists():
        return
    try:
        data = json.loads(jobs_path.read_text(encoding="utf-8"))
    except Exception:
        return
    jobs = data.get("jobs")
    if not isinstance(jobs, list):
        return
    for rec in jobs:
        if not isinstance(rec, dict):
            continue
        if str(rec.get("backend", "")).lower() != "lsf":
            continue
        st = str(rec.get("last_state", "")).upper()
        if st in {"DONE", "EXIT", "KILLED"}:
            continue
        jid = rec.get("job_id")
        if jid is None:
            continue
        run_cmd(["bkill", str(jid)], capture=True, check=False)


def find_main_start_with_autofer(
    *,
    ctx: SizeContext,
    python_bin: str,
    auto_fer_eval: Path,
    dry_run: bool,
) -> float:
    if dry_run:
        warn(f"M={ctx.m} dry-run: use main_start=snr_min from toml")
        txt = ctx.probe_toml_path.read_text(encoding="utf-8")
        m = re.search(r"^\s*snr_min\s*=\s*([0-9eE+\-\.]+)\s*$", txt, re.M)
        if not m:
            raise SimpleBatchError("Cannot parse snr_min from probe toml in dry-run mode")
        return float(m.group(1))

    cmd = [python_bin, "-u", str(auto_fer_eval), "adaptive", "--config", str(ctx.probe_toml_path)]
    log(f"M={ctx.m} -> running probe for main_start: {' '.join(cmd)}")
    proc = subprocess.Popen(
        cmd,
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
    )

    main_start: float | None = None
    start_re = re.compile(r"pilot start:\s*axis=snr\s*x[^0-9+\-]*([0-9eE+\-\.]+)")
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line, end="")
        m = start_re.search(line)
        if m:
            try:
                main_start = float(m.group(1))
            except ValueError:
                main_start = None
            if main_start is not None and math.isfinite(main_start):
                break

    if main_start is None:
        rc = proc.wait()
        raise SimpleBatchError(f"M={ctx.m} cannot parse main_start from auto_fer output (rc={rc})")

    if proc.poll() is None:
        try:
            proc.terminate()
            proc.wait(timeout=3.0)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait(timeout=3.0)
    kill_probe_lsf_jobs_best_effort(ctx)
    log(f"M={ctx.m} -> got main_start={format_snr_value(main_start)} and stopped auto_fer_eval")
    return main_start


def parse_log_metrics(text: str) -> tuple[float | None, float | None]:
    raw = _last_float(STAT_RAW_BER_RE, text)
    fer = _last_float(STAT_LDPC_FER_RE, text)
    if raw is None:
        raw = _last_float(SIM_RAW_BER_RE, text)
    if fer is None:
        fer = _last_float(SIM_LDPC_FER_RE, text)
    return raw, fer


def _last_float(regex: re.Pattern[str], text: str) -> float | None:
    ms = regex.findall(text)
    if not ms:
        return None
    try:
        return float(ms[-1])
    except ValueError:
        return None


def has_complete_statistics(log_path: Path) -> bool:
    if not log_path.exists():
        return False
    try:
        txt = log_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False
    return STAT_LDPC_FER_RE.search(txt) is not None


def parse_completed_log(log_path: Path) -> tuple[float | None, float | None] | None:
    if not has_complete_statistics(log_path):
        return None
    txt = log_path.read_text(encoding="utf-8", errors="replace")
    return parse_log_metrics(txt)


def bpeek_metrics(job_id: int) -> tuple[float | None, float | None]:
    p = run_cmd(["bpeek", str(job_id)], capture=True, check=False)
    out = p.stdout or ""
    if not out.strip():
        return None, None
    return parse_log_metrics(out)


def make_job_spec(
    *,
    ctx: SizeContext,
    qc: int,
    matrix_id: int,
    snr: float,
    exe: Path,
    matrix_root: Path,
) -> JobSpec:
    tag = format_snr_tag(snr)
    job_name = f"qc{qc}_M{ctx.m}_id{matrix_id}_snr{tag}"
    log_path = (ctx.output_dir / "logs" / f"{job_name}.log").resolve()
    cmd = [
        str(exe),
        "LDPC",
        str(ctx.config_path),
        "AWGN",
        format_snr_value(snr),
        str(matrix_id),
        str(matrix_root),
    ]
    return JobSpec(
        matrix_id=matrix_id,
        snr=snr,
        job_name=job_name,
        log_path=log_path,
        cmd=cmd,
    )


def poll_bjobs(job_id: int, log_path: Path) -> str:
    p = run_cmd(["bjobs", str(job_id)], capture=True, check=False)
    out = (p.stdout or "").strip()
    if p.returncode != 0 or not out:
        if has_complete_statistics(log_path):
            return "DONE"
        if log_path.exists() and log_path.stat().st_size > 0:
            return "RUN"
        return "UNKNOWN"

    lines = [ln.strip() for ln in out.splitlines() if ln.strip()]
    for ln in reversed(lines):
        if ln.startswith("JOBID"):
            continue
        toks = ln.split()
        if not toks:
            continue
        if toks[0] != str(job_id):
            continue
        return toks[2].upper() if len(toks) >= 3 else "UNKNOWN"
    return "UNKNOWN"


def submit_job(spec: JobSpec, *, queue: str, dry_run: bool) -> int:
    cmd = ["bsub", "-q", queue, "-J", spec.job_name, "-o", str(spec.log_path)] + spec.cmd
    if dry_run:
        log(f"[dry-run] {' '.join(cmd)}")
        return -1
    p = run_cmd(cmd, capture=True, check=True)
    out = p.stdout or ""
    m = re.search(r"Job <([0-9]+)>", out)
    if not m:
        raise SimpleBatchError(f"bsub output not recognized for {spec.job_name}: {out}")
    jid = int(m.group(1))
    log(f"submitted: {spec.job_name} -> job_id={jid}")
    return jid


def snr_key(v: float) -> float:
    return round(float(v), 6)


def build_snr_plan(main_start: float) -> tuple[list[float], list[float]]:
    base = [main_start + 0.1 * i for i in range(5)]  # +0.0 ~ +0.4
    extra = [main_start + 0.35, main_start + 0.45]
    all_pts = sorted({snr_key(x) for x in (base + extra)})
    rank_pts = [snr_key(main_start + 0.35), snr_key(main_start + 0.4), snr_key(main_start + 0.45)]
    return all_pts, rank_pts


def write_snr_plan_json(ctx: SizeContext, main_start: float, all_points: list[float], rank_points: list[float]) -> None:
    payload = {
        "size": ctx.size_name,
        "m": ctx.m,
        "n": ctx.n,
        "main_start": main_start,
        "simulate_points": all_points,
        "ranking_points": rank_points,
    }
    out = ctx.output_dir / "snr_points_simple.json"
    out.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    log(f"M={ctx.m} -> wrote snr plan: {out}")


def ranking_rows_for_snr(
    *,
    ctx: SizeContext,
    qc: int,
    matrix_ids: list[int],
    snr: float,
    exe: Path,
    matrix_root: Path,
    in_flight: dict[tuple[int, float], dict[str, Any]],
    allow_bpeek: bool,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    snr_k = snr_key(snr)
    for mid in matrix_ids:
        spec = make_job_spec(ctx=ctx, qc=qc, matrix_id=mid, snr=snr, exe=exe, matrix_root=matrix_root)
        row: dict[str, Any] = {
            "matrix_id": mid,
            "snr": snr,
            "ldpc_fer": None,
            "raw_ber": None,
            "status": "PENDING",
            "source": "",
            "job_id": "",
            "log_path": str(spec.log_path),
        }

        done = parse_completed_log(spec.log_path)
        if done is not None:
            raw_ber, ldpc_fer = done
            row["raw_ber"] = raw_ber
            row["ldpc_fer"] = ldpc_fer
            row["status"] = "DONE"
            row["source"] = "log"
            rows.append(row)
            continue

        info = in_flight.get((mid, snr_k))
        if info is None:
            rows.append(row)
            continue
        row["job_id"] = info.get("job_id", "")
        row["status"] = "RUNNING"
        if allow_bpeek:
            jid = info.get("job_id")
            if isinstance(jid, int) and jid > 0:
                raw_ber, ldpc_fer = bpeek_metrics(jid)
                if ldpc_fer is not None:
                    row["raw_ber"] = raw_ber
                    row["ldpc_fer"] = ldpc_fer
                    row["source"] = "bpeek"
        rows.append(row)

    def sort_key(r: dict[str, Any]) -> tuple[int, float]:
        fer = r.get("ldpc_fer")
        if fer is None:
            return (1, math.inf)
        try:
            return (0, float(fer))
        except (TypeError, ValueError):
            return (1, math.inf)

    rows.sort(key=sort_key)
    rank = 0
    for r in rows:
        if r.get("ldpc_fer") is not None:
            rank += 1
            r["rank"] = rank
        else:
            r["rank"] = ""
    return rows


def write_ranking_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "matrix_id", "snr", "ldpc_fer", "raw_ber", "status", "source", "job_id", "log_path", "updated_at"])
        now = _ts()
        for r in rows:
            w.writerow(
                [
                    r.get("rank", ""),
                    r.get("matrix_id", ""),
                    r.get("snr", ""),
                    r.get("ldpc_fer", ""),
                    r.get("raw_ber", ""),
                    r.get("status", ""),
                    r.get("source", ""),
                    r.get("job_id", ""),
                    r.get("log_path", ""),
                    now,
                ]
            )


def update_rankings(
    *,
    ctx: SizeContext,
    qc: int,
    matrix_ids: list[int],
    rank_points: list[float],
    exe: Path,
    matrix_root: Path,
    in_flight: dict[tuple[int, float], dict[str, Any]],
    realtime: bool,
    allow_bpeek: bool,
) -> None:
    for snr in rank_points:
        rows = ranking_rows_for_snr(
            ctx=ctx,
            qc=qc,
            matrix_ids=matrix_ids,
            snr=snr,
            exe=exe,
            matrix_root=matrix_root,
            in_flight=in_flight,
            allow_bpeek=allow_bpeek,
        )
        tag = format_snr_tag(snr)
        if realtime:
            out = ctx.output_dir / "summary" / "realtime" / f"ranking_snr{tag}.csv"
        else:
            out = ctx.output_dir / "summary" / f"ranking_snr{tag}.csv"
        write_ranking_csv(out, rows)
        if realtime:
            top = rows[0] if rows else {}
            log(
                f"{ctx.size_name} realtime rank@{format_snr_value(snr)} "
                f"top=id{top.get('matrix_id','?')} fer={top.get('ldpc_fer', '?')}"
            )


def write_points_csv(
    *,
    ctx: SizeContext,
    qc: int,
    matrix_ids: list[int],
    snr_points: list[float],
    exe: Path,
    matrix_root: Path,
) -> tuple[Path, dict[int, list[tuple[float, float, float]]]]:
    all_rows: list[dict[str, Any]] = []
    grouped: dict[int, list[tuple[float, float, float]]] = {}

    for matrix_id in matrix_ids:
        for snr in snr_points:
            spec = make_job_spec(
                ctx=ctx,
                qc=qc,
                matrix_id=matrix_id,
                snr=snr,
                exe=exe,
                matrix_root=matrix_root,
            )
            parsed = parse_completed_log(spec.log_path)
            if parsed is None:
                continue
            raw_ber, ldpc_fer = parsed
            if raw_ber is None or ldpc_fer is None:
                continue
            all_rows.append(
                {
                    "matrix_id": matrix_id,
                    "snr": snr,
                    "raw_ber": raw_ber,
                    "ldpc_fer": ldpc_fer,
                    "log_path": str(spec.log_path),
                }
            )
            grouped.setdefault(matrix_id, []).append((snr, raw_ber, ldpc_fer))

    out_csv = ctx.output_dir / "summary" / "rber_fer_points.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["matrix_id", "snr", "raw_ber", "ldpc_fer", "log_path"])
        for r in sorted(all_rows, key=lambda x: (int(x["matrix_id"]), float(x["snr"]))):
            w.writerow([r["matrix_id"], r["snr"], r["raw_ber"], r["ldpc_fer"], r["log_path"]])
    log(f"{ctx.size_name} -> wrote rber/fer points: {out_csv}")
    return out_csv, grouped


def plot_rber_fer(
    *,
    ctx: SizeContext,
    grouped_points: dict[int, list[tuple[float, float, float]]],
) -> Path | None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        warn("matplotlib not available, skip plotting")
        marker = ctx.output_dir / "plots" / "plot_skipped.txt"
        marker.write_text("matplotlib not available\n", encoding="utf-8")
        return None

    if not grouped_points:
        warn(f"{ctx.size_name} has no complete points for plotting")
        return None

    fig_path = ctx.output_dir / "plots" / f"rber_fer_M{ctx.m}.png"
    plt.figure(figsize=(8, 6))
    for matrix_id in sorted(grouped_points.keys()):
        pts = grouped_points[matrix_id]
        pts_sorted = sorted(pts, key=lambda x: x[1], reverse=True)
        xs = [p[1] for p in pts_sorted]
        ys = [p[2] for p in pts_sorted]
        plt.plot(xs, ys, marker="o", linewidth=1.2, markersize=4, label=f"id{matrix_id}")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("RAW BER")
    plt.ylabel("LDPC FER")
    plt.title(f"{ctx.size_name} RBER-FER")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=180)
    plt.close()
    log(f"{ctx.size_name} -> wrote plot: {fig_path}")
    return fig_path


def submit_and_wait(
    *,
    ctx: SizeContext,
    qc: int,
    matrix_ids: list[int],
    snr_points: list[float],
    rank_points: list[float],
    exe: Path,
    matrix_root: Path,
    queue: str,
    max_in_flight: int,
    poll_sec: float,
    bpeek_interval_sec: float,
    dry_run: bool,
) -> None:
    pending: list[JobSpec] = []
    for mid in matrix_ids:
        for snr in snr_points:
            spec = make_job_spec(ctx=ctx, qc=qc, matrix_id=mid, snr=snr, exe=exe, matrix_root=matrix_root)
            if has_complete_statistics(spec.log_path):
                continue
            pending.append(spec)

    if dry_run:
        for spec in pending:
            submit_job(spec, queue=queue, dry_run=True)
        return

    in_flight_by_jid: dict[int, JobSpec] = {}
    in_flight: dict[tuple[int, float], dict[str, Any]] = {}
    max_in_flight = max(1, int(max_in_flight))
    poll_sec = max(0.2, float(poll_sec))
    bpeek_interval_sec = max(1.0, float(bpeek_interval_sec))
    next_bpeek = time.time() + bpeek_interval_sec

    update_rankings(
        ctx=ctx,
        qc=qc,
        matrix_ids=matrix_ids,
        rank_points=rank_points,
        exe=exe,
        matrix_root=matrix_root,
        in_flight=in_flight,
        realtime=True,
        allow_bpeek=False,
    )

    while pending or in_flight_by_jid:
        while pending and len(in_flight_by_jid) < max_in_flight:
            spec = pending.pop(0)
            jid = submit_job(spec, queue=queue, dry_run=False)
            in_flight_by_jid[jid] = spec
            in_flight[(spec.matrix_id, snr_key(spec.snr))] = {"job_id": jid, "spec": spec}

        done_any = False
        for jid, spec in list(in_flight_by_jid.items()):
            st = poll_bjobs(jid, spec.log_path)
            if st in {"DONE", "EXIT", "ZOMBI"}:
                done_any = True
                in_flight_by_jid.pop(jid, None)
                in_flight.pop((spec.matrix_id, snr_key(spec.snr)), None)
                log(f"job done: {spec.job_name} status={st}")

        now = time.time()
        if now >= next_bpeek:
            update_rankings(
                ctx=ctx,
                qc=qc,
                matrix_ids=matrix_ids,
                rank_points=rank_points,
                exe=exe,
                matrix_root=matrix_root,
                in_flight=in_flight,
                realtime=True,
                allow_bpeek=True,
            )
            next_bpeek = now + bpeek_interval_sec
        elif done_any:
            update_rankings(
                ctx=ctx,
                qc=qc,
                matrix_ids=matrix_ids,
                rank_points=rank_points,
                exe=exe,
                matrix_root=matrix_root,
                in_flight=in_flight,
                realtime=True,
                allow_bpeek=False,
            )

        if not done_any and in_flight_by_jid:
            time.sleep(poll_sec)

    update_rankings(
        ctx=ctx,
        qc=qc,
        matrix_ids=matrix_ids,
        rank_points=rank_points,
        exe=exe,
        matrix_root=matrix_root,
        in_flight={},
        realtime=False,
        allow_bpeek=False,
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Simplified batch simulation with fixed points from main_start")
    ap.add_argument("--template-cnfg", required=True)
    ap.add_argument("--matrix-root", required=True)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--out-dir", default="output/qc192_batch_simple_out")
    ap.add_argument("--m-range", default="19-45")
    ap.add_argument("--n-matrices", type=int, default=5)
    ap.add_argument("--probe-id", type=int, default=1)
    ap.add_argument("--qc", type=int, default=192)
    ap.add_argument("--queue", default="regr_q")

    ap.add_argument("--auto-fer-eval", default=str(REPO_ROOT / "scripts" / "auto_fer_eval" / "auto_fer_eval.py"))
    ap.add_argument("--python-bin", default=sys.executable)
    ap.add_argument("--snr-min", type=float, default=3.0)
    ap.add_argument("--snr-max", type=float, default=7.0)
    ap.add_argument("--main-step", type=float, default=0.1)
    ap.add_argument("--pilot-max-sim", type=int, default=100)
    ap.add_argument("--main-max-sim", type=int, default=1000)
    ap.add_argument("--probe-max-in-flight", type=int, default=20)
    ap.add_argument("--probe-poll-sec", type=float, default=30.0)
    ap.add_argument("--probe-timeout-sec", type=float, default=0.0)
    ap.add_argument("--probe-log-grace-sec", type=float, default=10.0)
    ap.add_argument("--force-probe", action="store_true")

    ap.add_argument("--batch-max-in-flight", type=int, default=20)
    ap.add_argument("--batch-poll-sec", type=float, default=30.0)
    ap.add_argument("--bpeek-interval-sec", type=float, default=3600.0)
    ap.add_argument("--parallel-m", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    return ap.parse_args(argv)


def ensure_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    template_cnfg = Path(args.template_cnfg).expanduser().resolve()
    matrix_root = Path(args.matrix_root).expanduser().resolve()
    exe = Path(args.exe).expanduser().resolve()
    auto_fer_eval = Path(args.auto_fer_eval).expanduser().resolve()

    if not template_cnfg.exists():
        raise SimpleBatchError(f"template cnfg not found: {template_cnfg}")
    if not matrix_root.exists():
        raise SimpleBatchError(f"matrix root not found: {matrix_root}")
    if not exe.exists():
        raise SimpleBatchError(f"exe not found: {exe}")
    if not auto_fer_eval.exists():
        raise SimpleBatchError(f"auto_fer_eval not found: {auto_fer_eval}")
    return template_cnfg, matrix_root, exe, auto_fer_eval


def process_one_m(
    *,
    m: int,
    args: argparse.Namespace,
    template_cnfg: Path,
    matrix_root: Path,
    exe: Path,
    auto_fer_eval: Path,
    out_root: Path,
) -> None:
    ctx = create_size_context(
        m=m,
        qc=int(args.qc),
        template_cnfg=template_cnfg,
        matrix_root=matrix_root,
        out_root=out_root,
    )
    matrix_ids = list(range(1, int(args.n_matrices) + 1))

    write_probe_toml(
        ctx=ctx,
        exe=exe,
        matrix_root=matrix_root,
        queue=args.queue,
        probe_id=int(args.probe_id),
        snr_min=float(args.snr_min),
        snr_max=float(args.snr_max),
        main_step=float(args.main_step),
        pilot_max_sim=int(args.pilot_max_sim),
        main_max_sim=int(args.main_max_sim),
        probe_max_in_flight=int(args.probe_max_in_flight),
        probe_poll_sec=float(args.probe_poll_sec),
        probe_timeout_sec=float(args.probe_timeout_sec),
        probe_log_grace_sec=float(args.probe_log_grace_sec),
    )

    main_start_cache = ctx.output_dir / "main_start.json"
    if args.force_probe or not main_start_cache.exists():
        main_start = find_main_start_with_autofer(
            ctx=ctx,
            python_bin=args.python_bin,
            auto_fer_eval=auto_fer_eval,
            dry_run=bool(args.dry_run),
        )
        main_start_cache.write_text(json.dumps({"main_start": main_start}, indent=2), encoding="utf-8")
    else:
        data = json.loads(main_start_cache.read_text(encoding="utf-8"))
        main_start = float(data["main_start"])
        log(f"M={m} -> reuse cached main_start={format_snr_value(main_start)}")

    snr_points, rank_points = build_snr_plan(main_start)
    write_snr_plan_json(ctx, main_start, snr_points, rank_points)
    log(
        f"M={m} -> simulate points: {', '.join(format_snr_value(x) for x in snr_points)}; "
        f"rank points: {', '.join(format_snr_value(x) for x in rank_points)}"
    )

    submit_and_wait(
        ctx=ctx,
        qc=int(args.qc),
        matrix_ids=matrix_ids,
        snr_points=snr_points,
        rank_points=rank_points,
        exe=exe,
        matrix_root=matrix_root,
        queue=str(args.queue),
        max_in_flight=int(args.batch_max_in_flight),
        poll_sec=float(args.batch_poll_sec),
        bpeek_interval_sec=float(args.bpeek_interval_sec),
        dry_run=bool(args.dry_run),
    )
    if args.dry_run:
        warn(f"M={m} dry-run: skip rber-fer export/plot")
    else:
        _, grouped = write_points_csv(
            ctx=ctx,
            qc=int(args.qc),
            matrix_ids=matrix_ids,
            snr_points=snr_points,
            exe=exe,
            matrix_root=matrix_root,
        )
        plot_rber_fer(ctx=ctx, grouped_points=grouped)
    log(f"M={m} -> done")


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    template_cnfg, matrix_root, exe, auto_fer_eval = ensure_inputs(args)
    m_values = parse_m_range(args.m_range)
    out_root = Path(args.out_dir).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    parallel_m = max(1, int(args.parallel_m))

    log(
        f"start simple-batch: m_range={m_values[0]}..{m_values[-1]} "
        f"queue={args.queue} parallel_m={parallel_m}"
    )

    def _run_one(m: int) -> None:
        process_one_m(
            m=m,
            args=args,
            template_cnfg=template_cnfg,
            matrix_root=matrix_root,
            exe=exe,
            auto_fer_eval=auto_fer_eval,
            out_root=out_root,
        )

    if parallel_m == 1:
        for m in m_values:
            _run_one(m)
    else:
        errs: list[tuple[int, str]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_m) as pool:
            future_to_m = {pool.submit(_run_one, m): m for m in m_values}
            for fut in concurrent.futures.as_completed(future_to_m):
                m = future_to_m[fut]
                try:
                    fut.result()
                except Exception as exc:
                    errs.append((m, str(exc)))
                    warn(f"M={m} failed: {exc}")
        if errs:
            msg = "; ".join([f"M={m}: {e}" for m, e in errs])
            raise SimpleBatchError(f"{len(errs)} M task(s) failed -> {msg}")

    log("all done")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except SimpleBatchError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
