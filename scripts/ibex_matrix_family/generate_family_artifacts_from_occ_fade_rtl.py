#!/usr/bin/env python3
"""
Generate WrapBaseDelta and FirstMaskShift from:

1. a fixed WrapBase table provided by the caller, and
2. `ldpc_matrix_assignments_maxk67.svh` containing payload/parity occupied/fade bitmaps.

This script does NOT search WrapBase. Instead it treats WrapBase as a fixed input
and rebuilds each row's shift chain from left to right:

    first_active_shift = (wrap_base[row] + delta[row]) mod Z
    next_active_shift  = (current_shift + delta[row]) mod Z

for every active CPM encountered in the cropped `(M, K)` matrix view.

After rebuilding the shift view, it emits:
  - WrapBase: the provided fixed table
  - WrapBaseDelta: Definition A
      min k >= 0 such that (last_element + k * delta) mod Z == 0
  - FirstMaskShift:
      (Z + first_element(last_row) - first_parity_CPM(last_row)) mod Z

Input assignment format is the SVH emitted by scripts such as
`generate_family_artifacts.py`, for example:

  /tmp/origin_dv_rtl_param/lut/ldpc_matrix_assignments_maxk67.svh

Supports both 13-row and 17-row family views. When the assignment file contains
a 17-rate family, `--max-m 13` can be used to clamp to the 13-row view.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any

from generate_family_artifacts import IBEX_DELTA, format_concat


DEFAULT_WRAP_BASE_13 = [0, 270, 415, 337, 182, 301, 445, 70, 490, 65, 85, 75, 216]
ASSIGNMENT_RE = re.compile(
    r"assign\s+ldpc_matrix_(payload|parity)_(occupied|fade)\[\s*(\d+)\]\[\s*(\d+)\]\s*=\s*(\d+)'b([01]+)\s*;"
)


@dataclass
class AssignmentRate:
    m: int
    n: int
    payload_cols: int
    z: int
    occupied: list[list[int]]
    fade: list[list[int]]


def parse_int_list(raw: str) -> list[int]:
    values = [int(token.strip()) for token in raw.split(",") if token.strip()]
    if not values:
        raise ValueError("Expected a non-empty comma-separated integer list")
    return values


def bitmap_bit(bitmap: int, row: int) -> int:
    return (bitmap >> row) & 1


def parse_assignments_svh(
    assignments_path: Path,
    payload_cols: int,
    z: int,
    min_m: int,
    max_m: int | None,
) -> list[AssignmentRate]:
    raw_maps: dict[int, dict[str, dict[str, dict[int, int]]]] = {}
    with assignments_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            match = ASSIGNMENT_RE.match(line)
            if not match:
                continue
            column_group, value_kind, m_raw, index_raw, _width_raw, bits = match.groups()
            m = int(m_raw)
            index = int(index_raw)
            if m < min_m:
                continue
            if max_m is not None and m > max_m:
                continue
            raw_maps.setdefault(
                m,
                {
                    "payload": {"occupied": {}, "fade": {}},
                    "parity": {"occupied": {}, "fade": {}},
                },
            )
            raw_maps[m][column_group][value_kind][index] = int(bits, 2)

    if not raw_maps:
        raise ValueError(f"No ldpc_matrix_assignments found in {assignments_path}")

    rates: list[AssignmentRate] = []
    for m in sorted(raw_maps):
        occupied = [[0 for _ in range(payload_cols + m)] for _ in range(m)]
        fade = [[0 for _ in range(payload_cols + m)] for _ in range(m)]
        payload_occ = raw_maps[m]["payload"]["occupied"]
        payload_fade = raw_maps[m]["payload"]["fade"]
        parity_occ = raw_maps[m]["parity"]["occupied"]
        parity_fade = raw_maps[m]["parity"]["fade"]

        for col in range(payload_cols):
            occ_bitmap = payload_occ.get(col, 0)
            fade_bitmap = payload_fade.get(col, 0)
            for row in range(m):
                occupied[row][col] = bitmap_bit(occ_bitmap, row)
                fade[row][col] = bitmap_bit(fade_bitmap, row)

        for parity_idx in range(m):
            col = payload_cols + parity_idx
            occ_bitmap = parity_occ.get(parity_idx, 0)
            fade_bitmap = parity_fade.get(parity_idx, 0)
            for row in range(m):
                occupied[row][col] = bitmap_bit(occ_bitmap, row)
                fade[row][col] = bitmap_bit(fade_bitmap, row)

        rates.append(
            AssignmentRate(
                m=m,
                n=payload_cols + m,
                payload_cols=payload_cols,
                z=z,
                occupied=occupied,
                fade=fade,
            )
        )
    return rates


def crop_rate_shape(rate: AssignmentRate, payload_cols: int) -> tuple[list[list[int]], list[list[int]]]:
    cols = payload_cols + rate.m
    occupied = [[0 for _ in range(cols)] for _ in range(rate.m)]
    fade = [[0 for _ in range(cols)] for _ in range(rate.m)]
    for row in range(rate.m):
        for col in range(cols):
            if col < payload_cols:
                src_col = col
            else:
                src_col = rate.n - cols + col
            occupied[row][col] = rate.occupied[row][src_col]
            fade[row][col] = rate.fade[row][src_col]
    return occupied, fade


def rebuild_shifts_from_wrap_base(
    rate: AssignmentRate,
    payload_cols: int,
    occupied: list[list[int]],
    fade: list[list[int]],
    delta: list[int],
    wrap_base: list[int],
    z: int,
) -> dict[str, Any]:
    cols = payload_cols + rate.m
    shifts = [[-1 for _ in range(cols)] for _ in range(rate.m)]
    first_element = [0 for _ in range(rate.m)]
    last_element = [0 for _ in range(rate.m)]
    first_parity_shift = [0 for _ in range(rate.m)]
    row_weight = [0 for _ in range(rate.m)]

    for row in range(rate.m):
        cursor = (wrap_base[row] + delta[row]) % z
        first_seen = False
        parity_seen = False
        for col in range(cols):
            if not (occupied[row][col] or fade[row][col]):
                continue
            shifts[row][col] = cursor
            if not first_seen:
                first_element[row] = cursor
                first_seen = True
            last_element[row] = cursor
            if col >= payload_cols and not parity_seen:
                first_parity_shift[row] = cursor
                parity_seen = True
            if occupied[row][col]:
                row_weight[row] += 1
            cursor = (cursor + delta[row]) % z
        if not first_seen:
            raise ValueError(f"No active CPM in row {row} for M={rate.m},K={payload_cols}")

    return {
        "occupied": occupied,
        "fade": fade,
        "shifts": shifts,
        "first_element": first_element,
        "last_element": last_element,
        "first_parity_shift": first_parity_shift,
        "row_weight": row_weight,
    }


def compute_wrap_base_delta_from_last_element(last_element: int, delta: int, z: int) -> int:
    if delta == 0:
        if last_element % z != 0:
            raise ValueError(f"delta is zero but last_element={last_element} does not wrap to 0")
        return 0

    for k in range(z):
        if (last_element + k * delta) % z == 0:
            return k
    raise ValueError(
        f"Unable to find wrap_base_delta for last_element={last_element}, delta={delta}, z={z}"
    )


def build_summary(
    assignments_path: Path,
    payload_cols: int,
    z: int,
    lut_k_min: int,
    lut_k_max: int,
    min_m: int,
    max_m: int | None,
    wrap_base: list[int],
) -> tuple[dict[str, Any], list[AssignmentRate], dict[tuple[int, int], dict[str, Any]]]:
    rates = parse_assignments_svh(
        assignments_path=assignments_path,
        payload_cols=payload_cols,
        z=z,
        min_m=min_m,
        max_m=max_m,
    )
    if not rates:
        raise ValueError("No eligible rates discovered after M filtering")

    discovered_max_m = max(rate.m for rate in rates)
    discovered_min_m = min(rate.m for rate in rates)
    family_max_m = discovered_max_m if max_m is None else max_m
    if len(wrap_base) < family_max_m:
        raise ValueError(
            f"wrap_base length {len(wrap_base)} is shorter than required max_m={family_max_m}"
        )

    delta = list(IBEX_DELTA[:family_max_m])
    fixed_wrap_base = list(wrap_base[:family_max_m])

    configs: dict[tuple[int, int], dict[str, Any]] = {}
    for rate in rates:
        for k in range(lut_k_min, lut_k_max + 1):
            occupied, fade = crop_rate_shape(rate, k)
            configs[(rate.m, k)] = rebuild_shifts_from_wrap_base(
                rate=rate,
                payload_cols=k,
                occupied=occupied,
                fade=fade,
                delta=delta,
                wrap_base=fixed_wrap_base,
                z=z,
            )

    wrap_base_delta_rows: list[list[int]] = []
    wrap_base_delta_debug: list[list[dict[str, Any]]] = []
    lut_labels_by_row: list[list[str]] = []
    for row in range(family_max_m):
        row_values: list[int] = []
        row_debug: list[dict[str, Any]] = []
        row_labels: list[str] = []
        start_m = max(min_m, row + 1)
        for m in range(start_m, family_max_m + 1):
            for k in range(lut_k_min, lut_k_max + 1):
                cfg = configs[(m, k)]
                delta_num = compute_wrap_base_delta_from_last_element(
                    last_element=cfg["last_element"][row],
                    delta=delta[row],
                    z=z,
                )
                row_values.append(delta_num)
                row_labels.append(f"M{m}_K{k}")
                row_debug.append(
                    {
                        "m": m,
                        "payload_cols": k,
                        "label": f"M{m}_K{k}",
                        "row": row,
                        "first_element": cfg["first_element"][row],
                        "last_element": cfg["last_element"][row],
                        "wrap_base_delta": delta_num,
                        "delta": delta[row],
                    }
                )
        wrap_base_delta_rows.append(row_values)
        wrap_base_delta_debug.append(row_debug)
        lut_labels_by_row.append(row_labels)

    first_mask_shift_by_m: dict[int, list[int]] = {}
    first_mask_shift_debug: dict[int, list[dict[str, Any]]] = {}
    for m in range(min_m, family_max_m + 1):
        row = m - 1
        values: list[int] = []
        debug_items: list[dict[str, Any]] = []
        for k in range(lut_k_min, lut_k_max + 1):
            cfg = configs[(m, k)]
            first_element = cfg["first_element"][row]
            parity_first = cfg["first_parity_shift"][row]
            first_mask_shift = (z + first_element - parity_first) % z
            values.append(first_mask_shift)
            debug_items.append(
                {
                    "m": m,
                    "payload_cols": k,
                    "label": f"M{m}_K{k}",
                    "row": row,
                    "first_element": first_element,
                    "parity_first": parity_first,
                    "first_mask_shift": first_mask_shift,
                }
            )
        first_mask_shift_by_m[m] = values
        first_mask_shift_debug[m] = debug_items

    return (
        {
            "assignments_svh": str(assignments_path.resolve()),
            "source_semantics": (
                "assignments svh only; occupied/fade bitmaps reconstructed from ldpc_matrix_* assignments"
            ),
            "rebuild_semantics": (
                "fixed WrapBase input; left-to-right rebuild with first_active = wrap_base + delta"
            ),
            "wrap_base_delta_semantics": (
                "Definition A: min k >= 0 such that (last_element + k*delta) mod Z == 0"
            ),
            "first_mask_shift_semantics": (
                "(Z + first_element(last_row) - first_parity_CPM(last_row)) mod Z"
            ),
            "requested_min_m": min_m,
            "requested_max_m": max_m if max_m is not None else discovered_max_m,
            "discovered_min_m": discovered_min_m,
            "discovered_max_m": discovered_max_m,
            "min_m": min_m,
            "max_m": family_max_m,
            "lut_k_min": lut_k_min,
            "lut_k_max": lut_k_max,
            "z": z,
            "delta": delta,
            "wrap_base": fixed_wrap_base,
            "wrap_base_delta": wrap_base_delta_rows,
            "wrap_base_delta_debug": wrap_base_delta_debug,
            "wrap_base_delta_labels_by_row": lut_labels_by_row,
            "first_mask_shift_by_m": first_mask_shift_by_m,
            "first_mask_shift_debug": first_mask_shift_debug,
        },
        rates,
        configs,
    )


def write_reconstructed_matrices(
    summary: dict[str, Any],
    configs: dict[tuple[int, int], dict[str, Any]],
    output_root: Path,
) -> list[Path]:
    matrix_dir = output_root / "lut" / "reconstructed_matrix"
    matrix_dir.mkdir(parents=True, exist_ok=True)

    output_paths: list[Path] = []
    for m in range(summary["min_m"], summary["max_m"] + 1):
        for k in range(summary["lut_k_min"], summary["lut_k_max"] + 1):
            cfg = configs[(m, k)]
            matrix_path = matrix_dir / f"M{m:02d}_K{k:02d}_matrix.csv"
            with matrix_path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerows(cfg["shifts"])
            output_paths.append(matrix_path)
    return output_paths


def build_wrap_base_delta_csv_layout(summary: dict[str, Any]) -> tuple[list[str], list[list[int]]]:
    labels = [
        f"M{m}_K{k}"
        for m in range(summary["min_m"], summary["max_m"] + 1)
        for k in range(summary["lut_k_min"], summary["lut_k_max"] + 1)
    ]
    rows: list[list[int]] = []
    for row in range(summary["max_m"]):
        values_by_label = {
            item["label"]: item["wrap_base_delta"]
            for item in summary["wrap_base_delta_debug"][row]
        }
        rows.append([values_by_label.get(label, 0) for label in labels])
    return labels, rows


def write_outputs(
    summary: dict[str, Any],
    configs: dict[tuple[int, int], dict[str, Any]],
    output_root: Path,
) -> list[Path]:
    lut_dir = output_root / "lut"
    lut_dir.mkdir(parents=True, exist_ok=True)

    json_path = lut_dir / "matrix_lut_summary_fixed_wrap_base.json"
    csv_path = lut_dir / "wrap_base_delta_table_fixed_wrap_base.csv"
    mask_csv_path = lut_dir / "first_mask_shift_table_fixed_wrap_base.csv"
    svh_path = lut_dir / "ldpc_matrix_lut_fixed_wrap_base.svh"
    max_wrap_delta_len = max(len(values) for values in summary["wrap_base_delta"])

    json_path.write_text(json.dumps(summary, indent=2) + "\n")

    csv_labels, csv_rows = build_wrap_base_delta_csv_layout(summary)
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row_slot"] + csv_labels)
        for row, values in enumerate(csv_rows):
            writer.writerow([f"row_slot{row}", *values])

    with mask_csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["m", "k64", "k65", "k66", "k67"])
        for m in range(summary["min_m"], summary["max_m"] + 1):
            writer.writerow([m, *summary["first_mask_shift_by_m"][m]])

    lines = [
        "// Auto-generated by scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py",
        "// Fixed-WrapBase, left-to-right rebuild semantics.",
        "// WrapBase is caller-provided and not searched by this script.",
        "// Shift chain per row: first_active = WrapBase[row] + Delta[row], then +Delta for each next active CPM.",
        "// WrapBaseDelta uses Definition A:",
        "// min k >= 0 such that (last_element + k*delta) mod Z == 0.",
        "// FirstMaskShift uses: (Z + first_element(last_row) - first_parity_CPM(last_row)) mod Z.",
        "",
        f"localparam int LDPC_MATRIX_MIN_M = {summary['min_m']};",
        f"localparam int LDPC_MATRIX_MAX_M = {summary['max_m']};",
        f"localparam int LDPC_MATRIX_LUT_K_MIN = {summary['lut_k_min']};",
        f"localparam int LDPC_MATRIX_LUT_K_MAX = {summary['lut_k_max']};",
        f"localparam logic [9*{summary['max_m']}-1:0] ldpc_matrix_delta_tmp_adj_ord = "
        f"{format_concat(summary['delta'], 9, reverse=False)};",
        f"localparam logic [9*{summary['max_m']}-1:0] ldpc_matrix_wrap_base_tmp_adj_ord = "
        f"{format_concat(summary['wrap_base'], 9, reverse=False)};",
        "",
    ]

    for m in range(summary["min_m"], summary["max_m"] + 1):
        values = summary["first_mask_shift_by_m"][m]
        lines.append(
            f"wire [{len(values)}*9-1:0] ldpc_matrix_first_mask_shift{m:02d} = "
            f"{format_concat(values, 9, reverse=False)};"
        )

    lines.append("")
    for row, values in enumerate(summary["wrap_base_delta"]):
        padded_values = values + [0] * (max_wrap_delta_len - len(values))
        lines.append(
            f"wire [6*{max_wrap_delta_len}-1:0] ldpc_matrix_wrap_base_deltas{row}_adj_ord = "
            f"{format_concat(padded_values, 6, reverse=False)};"
        )

    svh_path.write_text("\n".join(lines) + "\n")
    matrix_paths = write_reconstructed_matrices(summary, configs, output_root)
    return [json_path, csv_path, mask_csv_path, svh_path, *matrix_paths]


def write_manifest(output_root: Path, summary: dict[str, Any], output_paths: list[Path]) -> Path:
    manifest_path = output_root / "manifest.json"
    manifest = {
        "assignments_svh": summary["assignments_svh"],
        "requested_min_m": summary["requested_min_m"],
        "requested_max_m": summary["requested_max_m"],
        "discovered_min_m": summary["discovered_min_m"],
        "discovered_max_m": summary["discovered_max_m"],
        "source_semantics": summary["source_semantics"],
        "rebuild_semantics": summary["rebuild_semantics"],
        "wrap_base_delta_semantics": summary["wrap_base_delta_semantics"],
        "first_mask_shift_semantics": summary["first_mask_shift_semantics"],
        "outputs": [str(path) for path in output_paths],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Generate WrapBaseDelta and FirstMaskShift from fixed WrapBase + "
            "ldpc_matrix_assignments_maxk67.svh."
        ),
        epilog=(
            "wrap_base format:\n"
            "  Comma-separated integers, with or without spaces.\n"
            "  Example for 13-row view:\n"
            "    --wrap-base \"0,270,415,337,182,301,445,70,490,65,85,75,216\"\n"
            "\n"
            "Examples:\n"
            "  13-row view (uses default wrap_base if --wrap-base is omitted):\n"
            "    python3 scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py \\\n"
            "      --assignments-svh /tmp/ibex_family_stage1/lut/ldpc_matrix_assignments_maxk67.svh \\\n"
            "      --output-root /tmp/ibex_family_stage2 \\\n"
            "      --max-m 13\n"
            "\n"
            "  17-row view (must pass a 17-element wrap_base):\n"
            "    python3 scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py \\\n"
            "      --assignments-svh /tmp/ibex_family_stage1/lut/ldpc_matrix_assignments_maxk67.svh \\\n"
            "      --output-root /tmp/ibex_family_stage2_17 \\\n"
            "      --max-m 17 \\\n"
            "      --wrap-base \"b0,b1,...,b16\""
        ),
    )
    parser.add_argument(
        "--assignments-svh",
        type=Path,
        default=Path("/tmp/origin_dv_rtl_param/lut/ldpc_matrix_assignments_maxk67.svh"),
        help="Input assignment SVH containing ldpc_matrix_payload/parity occupied/fade assigns.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("IBEX/ibex_matrix_fixed_wrap_base_lut"),
        help="Output directory for generated artifacts.",
    )
    parser.add_argument("--payload-cols", type=int, default=67, help="Full-width payload column count.")
    parser.add_argument("--circulant-size", type=int, default=512, help="Circulant size Z.")
    parser.add_argument("--lut-k-min", type=int, default=64)
    parser.add_argument("--lut-k-max", type=int, default=67)
    parser.add_argument("--min-m", type=int, default=5)
    parser.add_argument(
        "--max-m",
        type=int,
        default=None,
        help="Optional M clamp. Use --max-m 13 to derive the 13-row family view from a 17-row input file.",
    )
    parser.add_argument(
        "--wrap-base",
        default=",".join(str(value) for value in DEFAULT_WRAP_BASE_13),
        help=(
            "Comma-separated WrapBase table. Spaces are allowed between entries. "
            "Length must be at least max_m. "
            "Default is the 13-row RTL WrapBase provided in the thread."
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    wrap_base = parse_int_list(args.wrap_base)
    summary, _, configs = build_summary(
        assignments_path=args.assignments_svh.resolve(),
        payload_cols=args.payload_cols,
        z=args.circulant_size,
        lut_k_min=args.lut_k_min,
        lut_k_max=args.lut_k_max,
        min_m=args.min_m,
        max_m=args.max_m,
        wrap_base=wrap_base,
    )

    output_root = args.output_root.resolve()
    output_paths = write_outputs(summary, configs, output_root)
    manifest_path = write_manifest(output_root, summary, output_paths)

    print(f"[fixed-wrap-base] assignments_svh: {summary['assignments_svh']}")
    print(
        f"[fixed-wrap-base] family view: M={summary['min_m']}..{summary['max_m']}, "
        f"K={summary['lut_k_min']}..{summary['lut_k_max']}"
    )
    print(f"[fixed-wrap-base] wrap_base: {summary['wrap_base']}")
    for path in output_paths:
        print(f"[fixed-wrap-base] wrote {path}")
    print(f"[fixed-wrap-base] wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
