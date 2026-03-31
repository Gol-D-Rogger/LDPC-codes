#!/usr/bin/env python3
"""
Auto-search WrapBase from assignments SVH so Definition-A wrap_base_deltas remain small.

This is a separate variant of generate_family_artifacts_from_occ_fade_rtl.py:

1. Parse `ldpc_matrix_assignments_maxk67.svh`
2. For each row-slot, search wrap_base in [0, Z) to minimize the family-wide
   maximum Definition-A wrap_base_delta over all `(M, K)` views that use the row.
3. Rebuild shifts with the chosen wrap_base and emit:
   - ldpc_matrix_wrap_base_tmp_adj_ord
   - ldpc_matrix_wrap_base_deltasX_adj_ord
   - ldpc_matrix_first_mask_shiftMM
   - reconstructed matrices

The search is motivated by the RTL logic in IBEX/src/ldpc_codec.cpp::ldpc_config(),
which shows that an overall row shift can be absorbed into wrap_base / delta_num.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from generate_family_artifacts import IBEX_DELTA, format_concat
from generate_family_artifacts_from_occ_fade_rtl import (
    AssignmentRate,
    build_wrap_base_delta_csv_layout,
    build_summary,
    compute_wrap_base_delta_from_last_element,
    crop_rate_shape,
    parse_assignments_svh,
    write_reconstructed_matrices,
)


def compute_rtl_shift_dist_from_ldpc_config(
    first_element: int,
    last_element: int,
    delta: int,
    z: int,
    max_delta_num: int = 50,
) -> tuple[int, int]:
    """
    Mirror IBEX/src/ldpc_codec.cpp::ldpc_config() RTL search:

      init_base = (first_element - delta + bits) % bits
      for j in [0, bits):
        shift_base = init_base + j
        for k in [0, 50):
          if ((last + shift_base + k*delta + j) % bits) == shift_base:
            shift_dist = j
            wrap_num_deltas = k
    """
    init_base = (first_element - delta + z) % z
    for j in range(z):
        shift_base = init_base + j
        for k in range(max_delta_num):
            if ((last_element + shift_base + k * delta + j) % z) == shift_base:
                return j, k
    return -1, -1


def search_wrap_base_by_row(
    rates: list[AssignmentRate],
    z: int,
    lut_k_min: int,
    lut_k_max: int,
    min_m: int,
    max_m: int,
    width_limit: int,
) -> tuple[list[int], list[dict[str, Any]]]:
    delta = list(IBEX_DELTA[:max_m])

    active_counts: dict[tuple[int, int], list[int]] = {}
    for rate in rates:
        for k in range(lut_k_min, lut_k_max + 1):
            occupied, fade = crop_rate_shape(rate, k)
            row_counts = []
            for row in range(rate.m):
                row_counts.append(
                    sum(1 for col in range(k + rate.m) if occupied[row][col] or fade[row][col])
                )
            active_counts[(rate.m, k)] = row_counts

    wrap_base: list[int] = []
    debug_rows: list[dict[str, Any]] = []
    for row in range(max_m):
        start_m = max(min_m, row + 1)
        cases: list[tuple[int, int, int]] = []
        for m in range(start_m, max_m + 1):
            for k in range(lut_k_min, lut_k_max + 1):
                count = active_counts[(m, k)][row]
                cases.append((m, k, count))

        if delta[row] == 0:
            fixed_values: list[dict[str, int]] = []
            for m, k, count in cases:
                fixed_values.append(
                    {
                        "m": m,
                        "k": k,
                        "active_count": count,
                        "last_element": 0,
                        "wrap_base_delta": 0,
                    }
                )
            wrap_base.append(0)
            debug_rows.append(
                {
                    "row": row,
                    "chosen_wrap_base": 0,
                    "score": {
                        "fits_width_limit": True,
                        "max_wrap_base_delta": 0,
                        "sum_wrap_base_delta": 0,
                        "candidate": 0,
                    },
                    "cases": fixed_values,
                }
            )
            continue

        best_tuple: tuple[int, int, int, int] | None = None
        best_base = 0
        best_values: list[dict[str, int]] = []

        for candidate in range(z):
            delta_values: list[dict[str, int]] = []
            max_delta = 0
            sum_delta = 0
            fits_limit = 1
            for m, k, count in cases:
                last_element = (candidate + count * delta[row]) % z
                delta_num = compute_wrap_base_delta_from_last_element(last_element, delta[row], z)
                delta_values.append(
                    {
                        "m": m,
                        "k": k,
                        "active_count": count,
                        "last_element": last_element,
                        "wrap_base_delta": delta_num,
                    }
                )
                max_delta = max(max_delta, delta_num)
                sum_delta += delta_num
                if delta_num > width_limit:
                    fits_limit = 0

            score = (1 - fits_limit, max_delta, sum_delta, candidate)
            if best_tuple is None or score < best_tuple:
                best_tuple = score
                best_base = candidate
                best_values = delta_values

        wrap_base.append(best_base)
        debug_rows.append(
            {
                "row": row,
                "chosen_wrap_base": best_base,
                "score": {
                    "fits_width_limit": best_tuple[0] == 0,
                    "max_wrap_base_delta": best_tuple[1],
                    "sum_wrap_base_delta": best_tuple[2],
                    "candidate": best_tuple[3],
                },
                "cases": best_values,
            }
        )

    return wrap_base, debug_rows


def write_outputs_autobase(
    summary: dict[str, Any],
    configs: dict[tuple[int, int], dict[str, Any]],
    output_root: Path,
) -> list[Path]:
    lut_dir = output_root / "lut"
    lut_dir.mkdir(parents=True, exist_ok=True)

    json_path = lut_dir / "matrix_lut_summary_autobase.json"
    csv_path = lut_dir / "wrap_base_delta_table_autobase.csv"
    shift_dist_csv_path = lut_dir / "rtl_shift_dist_table_autobase.csv"
    mask_csv_path = lut_dir / "first_mask_shift_table_autobase.csv"
    svh_path = lut_dir / "ldpc_matrix_lut_autobase.svh"
    max_wrap_delta_len = max(len(values) for values in summary["wrap_base_delta"])

    json_path.write_text(json.dumps(summary, indent=2) + "\n")

    csv_labels, csv_rows = build_wrap_base_delta_csv_layout(summary)
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row_slot"] + csv_labels)
        for row, values in enumerate(csv_rows):
            writer.writerow([f"row_slot{row}", *values])

    shift_dist_rows: list[list[int]] = []
    for row in range(summary["max_m"]):
        values_by_label: dict[str, int] = {}
        for item in summary["wrap_base_delta_debug"][row]:
            shift_dist, _delta_num = compute_rtl_shift_dist_from_ldpc_config(
                first_element=item["first_element"],
                last_element=item["last_element"],
                delta=item["delta"],
                z=summary["z"],
            )
            values_by_label[item["label"]] = shift_dist
        shift_dist_rows.append([values_by_label.get(label, 0) for label in csv_labels])

    with shift_dist_csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row_slot"] + csv_labels)
        for row, values in enumerate(shift_dist_rows):
            writer.writerow([f"row_slot{row}", *values])

    with mask_csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["m", "k64", "k65", "k66", "k67"])
        for m in range(summary["min_m"], summary["max_m"] + 1):
            writer.writerow([m, *summary["first_mask_shift_by_m"][m]])

    lines = [
        "// Auto-generated by scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl_autobase.py",
        "// WrapBase is auto-searched per row-slot to minimize family-wide Definition-A wrap_base_delta.",
        "// Shift chain per row: first_active = WrapBase[row] + Delta[row], then +Delta for each next active CPM.",
        "// WrapBaseDelta uses Definition A:",
        "// min k >= 0 such that (last_element + k*delta) mod Z == 0.",
        f"// Width padding is applied so every wrap_base_deltas row is [6*{max_wrap_delta_len}-1:0].",
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
    return [json_path, csv_path, shift_dist_csv_path, mask_csv_path, svh_path, *matrix_paths]


def write_manifest_autobase(output_root: Path, summary: dict[str, Any], output_paths: list[Path]) -> Path:
    manifest_path = output_root / "manifest.json"
    manifest = {
        "assignments_svh": summary["assignments_svh"],
        "source_semantics": summary["source_semantics"],
        "rebuild_semantics": summary["rebuild_semantics"],
        "wrap_base_delta_semantics": summary["wrap_base_delta_semantics"],
        "first_mask_shift_semantics": summary["first_mask_shift_semantics"],
        "wrap_base_search": summary["wrap_base_search"],
        "outputs": [str(path) for path in output_paths],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Auto-search wrap_base from assignments SVH to reduce Definition-A wrap_base_delta width."
    )
    parser.add_argument(
        "--assignments-svh",
        type=Path,
        default=Path("/tmp/origin_dv_rtl_param/lut/ldpc_matrix_assignments_maxk67.svh"),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("IBEX/ibex_matrix_autobase_lut"),
    )
    parser.add_argument("--payload-cols", type=int, default=67)
    parser.add_argument("--circulant-size", type=int, default=512)
    parser.add_argument("--lut-k-min", type=int, default=64)
    parser.add_argument("--lut-k-max", type=int, default=67)
    parser.add_argument("--min-m", type=int, default=5)
    parser.add_argument("--max-m", type=int, default=None)
    parser.add_argument(
        "--width-limit",
        type=int,
        default=63,
        help="Preferred maximum wrap_base_delta value. Search minimizes overflow first, then max value.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rates = parse_assignments_svh(
        assignments_path=args.assignments_svh.resolve(),
        payload_cols=args.payload_cols,
        z=args.circulant_size,
        min_m=args.min_m,
        max_m=args.max_m,
    )
    family_max_m = max(rate.m for rate in rates) if args.max_m is None else args.max_m

    wrap_base, search_debug = search_wrap_base_by_row(
        rates=rates,
        z=args.circulant_size,
        lut_k_min=args.lut_k_min,
        lut_k_max=args.lut_k_max,
        min_m=args.min_m,
        max_m=family_max_m,
        width_limit=args.width_limit,
    )

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
    summary["wrap_base_search"] = {
        "strategy": "per-row brute force search over wrap_base in [0, Z)",
        "width_limit": args.width_limit,
        "rows": search_debug,
    }

    output_root = args.output_root.resolve()
    output_paths = write_outputs_autobase(summary, configs, output_root)
    manifest_path = write_manifest_autobase(output_root, summary, output_paths)

    max_delta = max(max(row) for row in summary["wrap_base_delta"])
    print(f"[autobase] assignments_svh: {summary['assignments_svh']}")
    print(
        f"[autobase] family view: M={summary['min_m']}..{summary['max_m']}, "
        f"K={summary['lut_k_min']}..{summary['lut_k_max']}"
    )
    print(f"[autobase] wrap_base: {summary['wrap_base']}")
    print(f"[autobase] max wrap_base_delta: {max_delta}")
    for path in output_paths:
        print(f"[autobase] wrote {path}")
    print(f"[autobase] wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
