"""
Job submission helpers for LDPC automation.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .commands import CommandRunner


@dataclass
class SimulationJob:
    matrix_id: str
    snr: float
    model: str
    queue: str
    job_name: str
    log_path: Path
    err_path: Optional[Path]
    command: str


def build_bsub_command(job: SimulationJob) -> List[str]:
    cmd = [
        "bsub",
        "-q",
        job.queue,
        "-J",
        job.job_name,
        "-o",
        str(job.log_path),
    ]
    if job.err_path:
        cmd.extend(["-e", str(job.err_path)])
    cmd.extend(["--", "bash", "-lc", job.command])
    return cmd


def submit_jobs(runner: CommandRunner, jobs: Iterable[SimulationJob]) -> None:
    for job in jobs:
        cmd = build_bsub_command(job)
        runner.run(cmd, capture_output=False, check=True)


def build_matrix_job_name(prefix: str, matrix_id: str, model: str, snr: Optional[float] = None) -> str:
    parts = [prefix, f"M{matrix_id}", model]
    if snr is not None:
        parts.append(f"S{snr}")
    return "_".join(parts)


def format_inner_command(
    exec_path: Path,
    config_path: Path,
    model: str,
    snr: Optional[float],
    matrix_id: str,
    matrix_dir: Optional[Path] = None,
    extra_args: Optional[Sequence[str]] = None,
) -> str:
    exec_dir = exec_path.parent
    exec_bin = exec_path.name
    cmd_parts = [f'cd "{exec_dir}"']
    env_cmd = ""
    if matrix_dir:
        env_cmd = f'export LDPC_MATRIX_DIR="{matrix_dir}"; '
    args = []
    if model.upper() == "AWGN":
        if snr is None:
            raise ValueError("AWGN 模式需要提供 SNR。")
        args = ["LDPC", str(config_path), "AWGN", str(snr), matrix_id]
        if matrix_dir:
            args.append(str(matrix_dir))
    else:
        args = ["LDPC", str(config_path), model, matrix_id]
        if snr is not None:
            args.insert(-1, str(snr))
    if extra_args:
        args.extend(extra_args)
    inner = " ".join(f'"{part}"' if " " in part else part for part in args)
    cmd_parts.append(f"{env_cmd}./\"{exec_bin}\" {inner}")
    return " && ".join(cmd_parts)


def format_custom_command(template: str, **kwargs: str) -> str:
    return template.format(**kwargs)
