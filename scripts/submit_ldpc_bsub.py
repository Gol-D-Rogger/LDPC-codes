#!/usr/bin/env python3
"""submit_ldpc_bsub.py

Python replacement of submit_ldpc_bsub.sh with clearer structure.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SNR_LIST = ["4.9", "5.0"]


def to_abs_path(path: Optional[str]) -> Optional[Path]:
    if not path:
        return None
    p = Path(path)
    if not p.is_absolute():
        p = WORKSPACE_ROOT / p
    return p.resolve()


def run_command(cmd: List[str], dry_run: bool = False, capture: bool = False,
                allow_failure: bool = False) -> str:
    if dry_run:
        print("[DRY-RUN]", " ".join(cmd))
        return ""
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        text=True,
        check=False,
    )
    if not allow_failure and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    if capture and result.stdout:
        return result.stdout
    return ""


@dataclass
class GenerationResult:
    matrix_dir: Path
    sim_root: Path
    new_ids: List[str]
    generated: bool


@dataclass
class Paths:
    exec_path: Path
    config_path: Path
    gen_bin: Path
    gen_out_base: Path
    matrix_dir: Path
    sim_matrix_root: Path
    out_dir: Path
    rename_script: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Submit LDPC simulations via LSF")
    parser.add_argument("--exec", default="gen4_ldpc_sim/ssd_fc_dq")
    parser.add_argument("--config", default="gen4_ldpc_sim/config/sdec_dq.cnfg")
    parser.add_argument("--gen-job-name")
    parser.add_argument("--matrix-dir")
    parser.add_argument("--out-dir", default="runs/ldpc_batch")
    parser.add_argument("--queue", default=os.environ.get("LSF_QUEUE", "regr_q"))
    parser.add_argument("--job-prefix", default="ldpc")
    parser.add_argument("--cwd")
    parser.add_argument("--snr")
    parser.add_argument("--snr-seq")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--gen-bin", default="GenLDPC/ldpc_gen")
    parser.add_argument("--gen-n", type=int)
    parser.add_argument("--gen-k", type=int)
    parser.add_argument("--gen-count", type=int)
    parser.add_argument("--gen-queue")
    parser.add_argument("--gen-mode", choices=["lsf", "local"], default="lsf")
    parser.add_argument("--gen-out-base")
    parser.add_argument("--rename-enable", action="store_true", default=True)
    parser.add_argument("--no-rename", action="store_true")
    parser.add_argument("--sim-matrix-dir")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def ensure_build(paths: Paths, build: bool, dry_run: bool) -> None:
    if not build:
        return
    run_command(["make", "-C", str(WORKSPACE_ROOT / "GenLDPC")], dry_run)
    run_command(["make", "-C", str(WORKSPACE_ROOT / "gen4_ldpc_sim")], dry_run)


def list_files(directory: Path, pattern: str) -> List[Path]:
    if not directory.exists():
        return []
    return sorted(directory.glob(pattern))


def wait_for_job_completion(job_name: str, dry_run: bool) -> None:
    if dry_run:
        return
    if shutil.which("bwait"):
        run_command(["bwait", "-w", f"ended({job_name})"], dry_run)
        return
    print("Waiting for generator jobs to finish (polling bjobs)...")
    while True:
        if shutil.which("bjobs"):
            out = run_command(["bjobs", "-J", job_name], dry_run=False, capture=True, allow_failure=True)
            rem = sum(1 for line in out.splitlines() if line.strip()) - 1
            if rem <= 0:
                break
        else:
            print("[WARN] bjobs not found; skip polling.", file=sys.stderr)
            break
        time.sleep(5)


def run_local_generation(gen_bin: Path, gen_out_base: Path, gen_n: int, gen_k: int,
                         gen_count: int, dry_run: bool) -> None:
    for phase in range(1, gen_count + 1):
        cmd = (
            f'GENLDPC_OUT_DIR="{gen_out_base}" "{gen_bin}" {gen_n} {gen_k} {phase} "{gen_out_base}"'
        )
        run_command(["bash", "-lc", cmd], dry_run)


def gather_new_ids(matrix_dir: Path, files_before: List[Path]) -> List[str]:
    files_after = [p.resolve() for p in list_files(matrix_dir, "*_QC_H_*.txt")]
    before_set = {p.resolve() for p in files_before}
    new_files = [str(p) for p in files_after if p not in before_set]

    new_ids: List[str] = []
    pattern_new = re.compile(r"_([0-9]+)\.txt$")
    for file_path in new_files:
        match = pattern_new.search(Path(file_path).name)
        if match:
            new_ids.append(match.group(1))
    if not new_ids:
        old_pattern = re.compile(r"_([0-9]+)_([0-9]+)\.txt$")
        for p in list_files(matrix_dir, "*_QC_H_*_*.txt"):
            match = old_pattern.search(p.name)
            if match:
                new_ids.append(f"{match.group(1)}_{match.group(2)}")
    return new_ids


def generate_matrices(args: argparse.Namespace, paths: Paths, dry_run: bool) -> GenerationResult:
    required = [args.gen_n, args.gen_k, args.gen_count]
    if any(v is None for v in required):
        raise SystemExit("生成矩阵需同时指定 --gen-n --gen-k --gen-count")

    mk = args.gen_n - args.gen_k
    matrix_dir = paths.gen_out_base / f"{mk}x{args.gen_n}" / "matrix"
    mask_dir = matrix_dir.parent / "mask_matrix"
    if not dry_run:
        matrix_dir.mkdir(parents=True, exist_ok=True)
        mask_dir.mkdir(parents=True, exist_ok=True)

    files_before = list_files(matrix_dir, "*_QC_H_*.txt")

    if args.gen_mode == "lsf":
        queue = args.gen_queue or args.queue
        job_name = args.gen_job_name or f"{args.job_prefix}_gen_{args.gen_n}x{args.gen_k}"
        log_dir = paths.out_dir / "logs_gen"
        if not dry_run:
            log_dir.mkdir(parents=True, exist_ok=True)
        gen_inner_cmd = (
            f'export GENLDPC_OUT_DIR="{paths.gen_out_base}"; '
            f'"{paths.gen_bin}" {args.gen_n} {args.gen_k} $LSB_JOBINDEX "{paths.gen_out_base}"'
        )
        bsub_cmd = [
            "bsub",
            "-q", queue,
            "-J", f"{job_name}[1-{args.gen_count}]",
            "-o", str(log_dir / "gen_%I.out"),
            "-e", str(log_dir / "gen_%I.err"),
            "--",
            "bash", "-lc", gen_inner_cmd,
        ]
        run_command(bsub_cmd, dry_run)
        wait_for_job_completion(job_name, dry_run)
    else:
        run_local_generation(paths.gen_bin, paths.gen_out_base, args.gen_n, args.gen_k, args.gen_count, dry_run)

    if args.rename_enable and not args.no_rename:
        rename_cmd = [
            "bash",
            str(paths.rename_script),
            str(mk),
            str(args.gen_n),
            "--apply",
            "--base-dir",
            str(paths.gen_out_base),
        ]
        run_command(rename_cmd, dry_run)

    new_ids = gather_new_ids(matrix_dir, files_before)
    generated = bool(new_ids)
    if generated:
        print(f"新增矩阵: {len(new_ids)} 个 -> {' '.join(new_ids)}")
    else:
        print("未检测到新增矩阵文件")
    return GenerationResult(matrix_dir=matrix_dir, sim_root=paths.gen_out_base, new_ids=new_ids, generated=generated)


def parse_snr_list(args: argparse.Namespace) -> List[str]:
    if args.snr:
        return args.snr.split()
    if args.snr_seq:
        try:
            start, step, end = map(float, args.snr_seq.split(":"))
        except ValueError:
            raise SystemExit("--snr-seq 需要格式 start:step:end")
        values = []
        x = start
        while x <= end + 1e-9:
            values.append(f"{x:.6f}")
            x += step
        return values
    return DEFAULT_SNR_LIST


def gather_matrix_ids(matrix_dir: Path, new_ids: List[str], generated: bool, allow_missing: bool = False) -> List[str]:
    def sort_key(value: str) -> List[int]:
        return [int(part) for part in value.split("_")]

    if generated and new_ids:
        return sorted(new_ids, key=sort_key)

    if not matrix_dir.exists():
        if allow_missing:
            return []
        raise SystemExit(f"矩阵目录不存在: {matrix_dir}")

    ids: List[str] = []
    seen = set()
    pattern_numeric = re.compile(r"_([0-9]+)(\.[^.]+)?$")
    pattern_trailing = re.compile(r"\.([0-9]+)$")
    for file_path in matrix_dir.glob("*.txt"):
        base = file_path.name
        match = pattern_numeric.search(base)
        if match:
            val = match.group(1)
        else:
            match = pattern_trailing.search(base)
            if not match:
                continue
            val = match.group(1)
        if val not in seen:
            ids.append(val)
            seen.add(val)
    if not ids:
        if allow_missing:
            return []
        raise SystemExit(f"矩阵目录为空：{matrix_dir}")
    return sorted(ids, key=sort_key)


def submit_simulations(args: argparse.Namespace, paths: Paths, matrix_ids: List[str],
                       snr_list: List[str], sim_root: Path) -> None:
    if not args.dry_run:
        paths.out_dir.mkdir(parents=True, exist_ok=True)

    exec_dir = paths.exec_path.parent
    exec_bin = paths.exec_path.name
    config_abs = str(paths.config_path)

    for matrix_id in matrix_ids:
        for snr in snr_list:
            log_dir = paths.out_dir / f"matrix{matrix_id}"
            if not args.dry_run:
                log_dir.mkdir(parents=True, exist_ok=True)
            log_o = log_dir / f"snr{snr}.out"
            log_e = log_dir / f"snr{snr}.err"

            if sim_root:
                inner_cmd = (
                    f'cd "{exec_dir}" && export LDPC_MATRIX_DIR="{sim_root}"; '
                    f'./"{exec_bin}" LDPC "{config_abs}" AWGN {snr} {matrix_id} "{sim_root}"'
                )
            else:
                inner_cmd = (
                    f'cd "{exec_dir}" && ./"{exec_bin}" LDPC "{config_abs}" AWGN {snr} {matrix_id}'
                )

            bsub_cmd = [
                "bsub",
                "-q", args.queue,
                "-J", f"{args.job_prefix}_M{matrix_id}_S{snr}",
                "-o", str(log_o),
                "-e", str(log_e),
            ]
            if args.cwd:
                bsub_cmd.extend(["-cwd", args.cwd])
            bsub_cmd.extend(["--", "bash", "-lc", inner_cmd])
            run_command(bsub_cmd, args.dry_run)


def prepare_paths(args: argparse.Namespace) -> Paths:
    exec_path = to_abs_path(args.exec)
    config_path = to_abs_path(args.config)
    gen_bin = to_abs_path(args.gen_bin)

    matrix_dir_arg = to_abs_path(args.matrix_dir) if args.matrix_dir else None
    matrix_dir = matrix_dir_arg if matrix_dir_arg else (WORKSPACE_ROOT / "gen4_ldpc_sim/matrix").resolve()

    if args.gen_out_base:
        gen_out_base = to_abs_path(args.gen_out_base)
    elif matrix_dir_arg is not None:
        gen_out_base = matrix_dir.parent
    else:
        gen_out_base = (WORKSPACE_ROOT / "GenLDPC/output").resolve()
    sim_matrix_root = to_abs_path(args.sim_matrix_dir) if args.sim_matrix_dir else matrix_dir.parent
    out_dir = to_abs_path(args.out_dir) or (WORKSPACE_ROOT / "runs/ldpc_batch").resolve()
    rename_script = WORKSPACE_ROOT / "scripts" / "rename_ldpc_matrices.sh"

    return Paths(
        exec_path=exec_path,
        config_path=config_path,
        gen_bin=gen_bin,
        gen_out_base=gen_out_base,
        matrix_dir=matrix_dir,
        sim_matrix_root=sim_matrix_root,
        out_dir=out_dir,
        rename_script=rename_script,
    )


def main() -> None:
    args = parse_args()
    if args.no_rename:
        args.rename_enable = False

    paths = prepare_paths(args)

    ensure_build(paths, args.build, args.dry_run)

    generation_result: Optional[GenerationResult] = None
    if any(v is not None for v in (args.gen_n, args.gen_k, args.gen_count)):
        generation_result = generate_matrices(args, paths, args.dry_run)
        matrix_dir = generation_result.matrix_dir
        sim_root = generation_result.sim_root
    else:
        matrix_dir = paths.matrix_dir
        sim_root = paths.sim_matrix_root

    snr_list = parse_snr_list(args)

    new_ids = generation_result.new_ids if generation_result else []
    generated = generation_result.generated if generation_result else False
    matrix_ids = gather_matrix_ids(matrix_dir, new_ids, generated, allow_missing=args.dry_run)
    if args.dry_run and not matrix_ids:
        fallback_dir = paths.matrix_dir
        if fallback_dir.exists():
            matrix_ids = gather_matrix_ids(fallback_dir, [], False, allow_missing=True)
        if not matrix_ids:
            matrix_ids = ["1"]

    print(f"Detected matrix IDs: {' '.join(matrix_ids)}")
    if matrix_ids:
        print(f"MATRIX_END = {matrix_ids[-1]}")
    print(f"SNRs: {' '.join(snr_list)}")
    print(f"Queue: {args.queue}")
    print(f"Exec:  {paths.exec_path}")
    print(f"Config:{paths.config_path}")
    print(f"Matrix:{matrix_dir}")
    print(f"Out:   {paths.out_dir}")
    print()

    submit_simulations(args, paths, matrix_ids, snr_list, sim_root)
    print(f"Done. Logs under: {paths.out_dir}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(f"Command failed: {' '.join(exc.cmd)}", file=sys.stderr)
        if exc.stderr:
            print(exc.stderr, file=sys.stderr)
        sys.exit(exc.returncode)
