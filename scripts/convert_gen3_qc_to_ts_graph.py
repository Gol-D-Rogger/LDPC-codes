#!/usr/bin/env python3
"""
Convert a standard Gen3 QC-LDPC base matrix to TS_enum `.graph` format.

Input matrix format:
  - whitespace-separated integers
  - `-1` means an all-zero ZxZ submatrix
  - `0..Z-1` means a ZxZ circulant permutation with the given shift

Output `.graph` format for TS_enum:
  - line 1: "<num_vn> <num_cn>"
  - line 2..: for each expanded CN row, a space-separated list of 1-based VN indices

Gen3 uses standard QC-LDPC structure only. There is no Gen4-style fade/mask logic here.
"""

from __future__ import annotations

import argparse
import re
import shlex
import subprocess
import sys
from pathlib import Path


def read_qc_matrix(path: Path) -> list[list[int]]:
    rows: list[list[int]] = []
    with path.open() as f:
        for line_no, line in enumerate(f, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                row = [int(token) for token in stripped.split()]
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: non-integer token") from exc
            rows.append(row)

    if not rows:
        raise ValueError(f"{path}: matrix file is empty")

    expected_cols = len(rows[0])
    for row_idx, row in enumerate(rows, start=1):
        if len(row) != expected_cols:
            raise ValueError(
                f"{path}:{row_idx}: expected {expected_cols} columns, got {len(row)}"
            )

    return rows


def infer_ex_factor(matrix_path: Path, explicit_value: int | None) -> int:
    if explicit_value is not None:
        return explicit_value

    match = re.search(r"ex(\d+)", matrix_path.name)
    if match:
        return int(match.group(1))

    raise ValueError(
        "Unable to infer ex-factor from filename. Pass --ex-factor explicitly."
    )


def default_graph_path(matrix_path: Path) -> Path:
    graph_name = matrix_path.stem.replace("_QC_H", "_graph") + ".graph"
    return matrix_path.with_name(graph_name)


def validate_matrix_values(matrix: list[list[int]], ex_factor: int, matrix_path: Path) -> None:
    for row_idx, row in enumerate(matrix, start=1):
        for col_idx, value in enumerate(row, start=1):
            if value < -1:
                raise ValueError(
                    f"{matrix_path}:{row_idx}:{col_idx}: invalid shift {value}; expected -1 or [0, {ex_factor - 1}]"
                )
            if value >= ex_factor:
                raise ValueError(
                    f"{matrix_path}:{row_idx}:{col_idx}: shift {value} exceeds ex-factor {ex_factor}"
                )


def write_ts_graph(matrix: list[list[int]], ex_factor: int, out_path: Path) -> tuple[int, int]:
    base_rows = len(matrix)
    base_cols = len(matrix[0])
    num_cn = base_rows * ex_factor
    num_vn = base_cols * ex_factor

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write(f"{num_vn} {num_cn}\n")
        for base_row in matrix:
            active_cols = [(col_idx, shift) for col_idx, shift in enumerate(base_row) if shift >= 0]
            for local_row in range(ex_factor):
                neighbors = [
                    str(col_idx * ex_factor + ((local_row + shift) % ex_factor) + 1)
                    for col_idx, shift in active_cols
                ]
                f.write(" ".join(neighbors) + "\n")

    return num_vn, num_cn


def run_ts_enum(
    graph_path: Path,
    trap_path: Path,
    ts_enum_path: Path,
    wine_cmd: str,
    maxit: int,
    qc: int,
    fast: bool,
) -> None:
    cmd = [
        wine_cmd,
        str(ts_enum_path),
        "-maxit",
        str(maxit),
        "-qc",
        str(qc),
    ]
    if fast:
        cmd.append("-fast")
    cmd.extend([str(graph_path), str(trap_path)])

    print("TS_enum command:")
    print("  " + shlex.join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a Gen3 QC-H matrix to TS_enum .graph format."
    )
    parser.add_argument("matrix", help="Path to Gen3 QC_H matrix file")
    parser.add_argument(
        "--ex-factor",
        type=int,
        default=None,
        help="Circulant size Z. If omitted, infer from filename like ex256.",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Output .graph path. Default: next to input matrix, with _QC_H -> _graph.",
    )
    parser.add_argument(
        "--run-ts",
        action="store_true",
        help="Run TS_enum immediately after generating the .graph file.",
    )
    parser.add_argument(
        "--ts-enum",
        default="TS/TS_enum.exe",
        help="Path to TS_enum.exe when --run-ts is used.",
    )
    parser.add_argument(
        "--wine-cmd",
        default="wine",
        help="Wine command used to run TS_enum.exe.",
    )
    parser.add_argument(
        "--trap-out",
        default=None,
        help="Output .trap path when --run-ts is used. Default: same basename as .graph.",
    )
    parser.add_argument(
        "--maxit",
        type=int,
        default=25,
        help="TS_enum -maxit value when --run-ts is used.",
    )
    parser.add_argument(
        "--qc",
        type=int,
        default=None,
        help="TS_enum -qc value when --run-ts is used. Default: ex-factor.",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Pass -fast to TS_enum when --run-ts is used.",
    )
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    if not matrix_path.is_file():
        raise SystemExit(f"[ERROR] matrix file not found: {matrix_path}")

    try:
        matrix = read_qc_matrix(matrix_path)
        ex_factor = infer_ex_factor(matrix_path, args.ex_factor)
        validate_matrix_values(matrix, ex_factor, matrix_path)
    except ValueError as exc:
        raise SystemExit(f"[ERROR] {exc}") from exc

    out_path = Path(args.out) if args.out else default_graph_path(matrix_path)
    num_vn, num_cn = write_ts_graph(matrix, ex_factor, out_path)

    print(f"Input matrix: {matrix_path}")
    print(f"Base shape: {len(matrix)} x {len(matrix[0])}")
    print(f"Ex-factor: {ex_factor}")
    print(f"Expanded graph: {num_cn} CN x {num_vn} VN")
    print(f"Graph file: {out_path}")

    if not args.run_ts:
        return

    ts_enum_path = Path(args.ts_enum)
    if not ts_enum_path.is_file():
        raise SystemExit(f"[ERROR] TS_enum not found: {ts_enum_path}")

    trap_path = Path(args.trap_out) if args.trap_out else out_path.with_suffix(".trap")
    qc = args.qc if args.qc is not None else ex_factor

    try:
        run_ts_enum(
            graph_path=out_path,
            trap_path=trap_path,
            ts_enum_path=ts_enum_path,
            wine_cmd=args.wine_cmd,
            maxit=args.maxit,
            qc=qc,
            fast=args.fast,
        )
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"[ERROR] TS_enum failed with exit code {exc.returncode}") from exc

    print(f"Trap file: {trap_path}")


if __name__ == "__main__":
    main()
