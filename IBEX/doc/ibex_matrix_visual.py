#!/usr/bin/env python3
"""
IBEX QC-LDPC Matrix Structure Visualizer
=========================================
Reads real matrix data from ibex_matrix/ and ibex_matrix_flat_13rate/,
generates publication-quality figures explaining:

  Fig 1 - 80-Column Internal Matrix Layout (13x80 full matrix heatmap)
  Fig 2 - 5x5 Parity Block Structure
  Fig 3 - Trimming Mechanism (how 80-col maps to different (M,K) configs)
  Fig 4 - Shift Assignment & Delta Stepping (per-row shift sequences)
  Fig 5 - LUT Structure (36-entry table indexed by extra_rows x extra_userdata_cols)
  Fig 6 - generate_family_artifacts.py Pipeline Overview
  Fig 7 - Wraparound Identity Verification
  Fig 8 - RTL Two-Stage Pipeline & Left-to-Right Shift Chain Semantics
  Fig 9 - WrapBaseDelta "Ring Track" Intuition & AutoBase Search
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent  # -> LDPC-codes
IBEX_ROOT = REPO_ROOT / "IBEX"
MATRIX_DIR = IBEX_ROOT / "ibex_matrix"
FLAT_DIR = IBEX_ROOT / "ibex_matrix_flat_13rate"
OUTPUT_DIR = IBEX_ROOT / "doc" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

Z = 512
BASE_PAYLOAD_COLS = 67
MAX_INTERNAL_COLS = 80
PARITY_5X5_START = 75
MIN_M = 5
MAX_M_13 = 13


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------
def read_int_matrix(path: Path) -> list[list[int]]:
    rows = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append([int(t) for t in line.split()])
    return rows


def load_bm_schematic(m: int, n: int) -> dict:
    path = FLAT_DIR / "bm_schematic" / f"bm_schematic_{m}x{n}ex{Z}_w4.txt"
    occupied = []
    fade = []
    in_fade_section = False
    in_merged_section = False
    with path.open() as f:
        for line in f:
            if "Fade matrix" in line:
                in_fade_section = True
                continue
            if "fade shown as" in line:
                in_merged_section = True
                in_fade_section = False
                continue
            if not line.startswith("ROW"):
                continue
            data_part = line.split(":")[1].strip().replace("|", "")
            if in_merged_section:
                break
            if in_fade_section:
                fade.append([1 if c == "1" else 0 for c in data_part])
            else:
                occupied.append([1 if c == "1" else 0 for c in data_part])
    return {"occupied": occupied, "fade": fade, "m": m, "n": n}


def load_lut_summary() -> dict:
    path = FLAT_DIR / "lut" / "matrix_lut_summary.json"
    with path.open() as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Figure 1: 80-Column Internal Matrix Layout
# ---------------------------------------------------------------------------
def fig1_matrix_layout():
    m, n = 13, 80
    matrix_path = MATRIX_DIR / f"{m}x{n}" / "matrix"
    shift_file = sorted(matrix_path.glob("*_QC_H_*.txt"))[0]
    occ_path = MATRIX_DIR / f"{m}x{n}" / "occupied_matrix"
    occ_file = sorted(occ_path.glob("*_occupied_*.txt"))[0]
    fade_path = MATRIX_DIR / f"{m}x{n}" / "fade_matrix"
    fade_file = sorted(fade_path.glob("*_fade_*.txt"))[0]

    shifts = read_int_matrix(shift_file)
    occ = read_int_matrix(occ_file)
    fade = read_int_matrix(fade_file)

    # Build color matrix: 0=empty, 1=occupied, 2=fade, 3=parity-5x5-occupied
    color_mat = np.zeros((m, n), dtype=float)
    for i in range(m):
        for j in range(n):
            if occ[i][j] == 1:
                if j >= PARITY_5X5_START:
                    color_mat[i][j] = 3   # 5x5 block
                elif j >= BASE_PAYLOAD_COLS:
                    color_mat[i][j] = 4   # diagonal parity
                else:
                    color_mat[i][j] = 1   # payload occupied
            elif fade[i][j] == 1:
                color_mat[i][j] = 2       # fade
            else:
                color_mat[i][j] = 0       # empty

    fig, ax = plt.subplots(figsize=(20, 5))
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#F0F0F0", "#2196F3", "#FF9800", "#E91E63", "#4CAF50"])
    ax.imshow(color_mat, cmap=cmap, aspect="auto", interpolation="nearest")

    # Zone boundaries
    for x in [BASE_PAYLOAD_COLS - 0.5, PARITY_5X5_START - 0.5]:
        ax.axvline(x, color="black", linewidth=2, linestyle="--")

    # Zone labels
    ax.text(33, -1.2, "Payload (col 0-66)\ncol_wt=4, 67 cols",
            ha="center", fontsize=10, fontweight="bold", color="#2196F3")
    ax.text(71, -1.2, "Diag Parity\n(col 67-74)",
            ha="center", fontsize=10, fontweight="bold", color="#4CAF50")
    ax.text(77.5, -1.2, "5x5 Block\n(col 75-79)",
            ha="center", fontsize=10, fontweight="bold", color="#E91E63")

    ax.set_xlabel("Internal Column Index", fontsize=12)
    ax.set_ylabel("Row Index", fontsize=12)
    ax.set_title("Figure 1: IBEX 13x80 Internal Matrix Layout (srand(1) generated)", fontsize=14, fontweight="bold")
    ax.set_xticks(range(0, 80, 5))
    ax.set_yticks(range(13))

    legend_elements = [
        mpatches.Patch(facecolor="#2196F3", label="Payload Occupied"),
        mpatches.Patch(facecolor="#FF9800", label="Fade"),
        mpatches.Patch(facecolor="#4CAF50", label="Diagonal Parity Occupied"),
        mpatches.Patch(facecolor="#E91E63", label="5x5 Parity Block"),
        mpatches.Patch(facecolor="#F0F0F0", edgecolor="gray", label="Empty"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9, ncol=3)

    fig.tight_layout()
    path = OUTPUT_DIR / "fig1_matrix_layout.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 2: 5x5 Parity Block Structure
# ---------------------------------------------------------------------------
def fig2_5x5_block():
    block = np.zeros((5, 5), dtype=int)
    labels = [["" for _ in range(5)] for _ in range(5)]
    for i in range(5):
        for j in range(5):
            k = 4 - j  # k = 79 - (75+j)
            use = (k % 5) != i
            if i == 4 and j == 4:
                use = False  # invertibility exclusion
            block[i][j] = 1 if use else 0
            if not use:
                if k % 5 == i:
                    labels[i][j] = "diag\nexcl"
                else:
                    labels[i][j] = "invert\nexcl"
            else:
                labels[i][j] = "1"

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [1, 1.5]})

    # Left: 5x5 block heatmap
    ax = axes[0]
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#FFCDD2", "#C8E6C9"])
    ax.imshow(block, cmap=cmap, aspect="equal", interpolation="nearest")
    for i in range(5):
        for j in range(5):
            color = "green" if block[i][j] else "red"
            ax.text(j, i, labels[i][j], ha="center", va="center",
                    fontsize=9, fontweight="bold", color=color)

    ax.set_xticks(range(5))
    ax.set_xticklabels([f"col {75+j}" for j in range(5)])
    ax.set_yticks(range(5))
    ax.set_yticklabels([f"row {i}" for i in range(5)])
    ax.set_title("5x5 Parity Block Occupied Pattern", fontsize=12, fontweight="bold")

    # Right: entry count per row
    ax2 = axes[1]
    entry_counts = [sum(block[i]) for i in range(5)]
    colors = ["#4CAF50" if c == 4 else "#FF9800" for c in entry_counts]
    bars = ax2.barh(range(5), entry_counts, color=colors, edgecolor="black")
    ax2.set_yticks(range(5))
    ax2.set_yticklabels([f"row {i}" for i in range(5)])
    ax2.set_xlabel("Entry Count in 5x5 Block")
    ax2.set_title("Entries per Row in 5x5 Block", fontsize=12, fontweight="bold")
    ax2.invert_yaxis()
    for i, (bar, cnt) in enumerate(zip(bars, entry_counts)):
        ax2.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
                 str(cnt), va="center", fontsize=11, fontweight="bold")
    ax2.set_xlim(0, 6)

    # Delta as identity label (separate table, no causal link to entry count)
    delta_info = "Row identity:  " + "  |  ".join(
        f"row {i}: delta={[0,13,19,29,41][i]}" for i in range(5)
    )
    note = ("Rule: occupied[i][j] = ((79-j) % 5) != i\n"
            "Exception: row 4, col 79 excluded for invertibility\n"
            f"{delta_info}\n"
            "Note: entry count is determined by the 5x5 geometry, NOT by delta")
    fig.text(0.5, -0.02, note, ha="center", fontsize=9, style="italic", color="gray")

    fig.suptitle("Figure 2: 5x5 Parity Block Structure", fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    path = OUTPUT_DIR / "fig2_5x5_block.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 3: Trimming Mechanism
# ---------------------------------------------------------------------------
def fig3_trimming():
    fig, axes = plt.subplots(4, 1, figsize=(20, 10), sharex=True)
    configs = [
        (8, 75, "M=8, N=75 (K=67)"),
        (8, 73, "M=8, N=73 (K=65)"),
        (11, 78, "M=11, N=78 (K=67)"),
        (11, 75, "M=11, N=75 (K=64)"),
    ]

    for ax, (m, n, title) in zip(axes, configs):
        k = n - m
        # Build column map: target col j -> internal col
        col_map = {}
        for j in range(n):
            if j < k:
                col_map[j] = j  # payload direct
            else:
                col_map[j] = 80 - n + j  # parity from right

        # Color the 80-column bar
        bar = np.zeros(80)
        for j, internal in col_map.items():
            if internal < BASE_PAYLOAD_COLS:
                bar[internal] = 1  # kept payload
            elif internal < PARITY_5X5_START:
                bar[internal] = 3  # kept diag parity
            else:
                bar[internal] = 4  # kept 5x5

        # Mark skipped
        for c in range(80):
            if bar[c] == 0:
                if c < BASE_PAYLOAD_COLS:
                    bar[c] = 0.3  # skipped payload
                elif c < PARITY_5X5_START:
                    bar[c] = 0.3  # skipped diag
                else:
                    bar[c] = 0.3  # skipped 5x5

        from matplotlib.colors import Normalize
        colors = []
        for c in range(80):
            v = bar[c]
            if v == 1:
                colors.append("#2196F3")    # kept payload
            elif v == 3:
                colors.append("#4CAF50")    # kept diag parity
            elif v == 4:
                colors.append("#E91E63")    # kept 5x5
            else:
                colors.append("#E0E0E0")    # skipped

        ax.bar(range(80), [1]*80, color=colors, edgecolor="white", linewidth=0.3)

        # Zone boundaries
        for x in [BASE_PAYLOAD_COLS - 0.5, PARITY_5X5_START - 0.5]:
            ax.axvline(x, color="black", linewidth=1.5, linestyle="--", alpha=0.5)

        # Annotations
        skipped_payload = [c for c in range(BASE_PAYLOAD_COLS) if bar[c] < 1]
        skipped_mid = [c for c in range(BASE_PAYLOAD_COLS, PARITY_5X5_START) if bar[c] < 1]
        kept_payload = k
        kept_parity = m

        info = f"{title}  |  Keep payload col 0..{kept_payload - 1}, parity (rightmost {m} cols)"
        skip_parts = []
        if skipped_payload:
            skip_parts.append(f"payload [{skipped_payload[0]}..{skipped_payload[-1]}]")
        if skipped_mid:
            skip_parts.append(f"diag [{skipped_mid[0]}..{skipped_mid[-1]}]")
        if skip_parts:
            info += f"  |  SKIP internal " + " + ".join(skip_parts)
        ax.set_title(info, fontsize=11, fontweight="bold", loc="left")
        ax.set_ylim(0, 1.2)
        ax.set_yticks([])

        # Arrow showing trimming
        if skipped_payload or skipped_mid:
            skip_all = skipped_payload + skipped_mid
            for sc in skip_all:
                ax.text(sc, 0.5, "X", ha="center", va="center",
                        fontsize=7, color="red", fontweight="bold")

    axes[-1].set_xlabel("Internal 80-Column Index", fontsize=12)
    axes[-1].set_xticks(range(0, 80, 5))

    legend_elements = [
        mpatches.Patch(facecolor="#2196F3", label="Kept Payload"),
        mpatches.Patch(facecolor="#4CAF50", label="Kept Diag Parity"),
        mpatches.Patch(facecolor="#E91E63", label="Kept 5x5 Block"),
        mpatches.Patch(facecolor="#E0E0E0", edgecolor="gray", label="Skipped (trimmed)"),
    ]
    axes[0].legend(handles=legend_elements, loc="upper right", fontsize=8, ncol=4)

    fig.suptitle("Figure 3: Trimming Mechanism - 80-Col Internal -> Target (M,K)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = OUTPUT_DIR / "fig3_trimming.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 4: Shift Assignment & Delta Stepping
# ---------------------------------------------------------------------------
def _read_origin_matrix() -> list[list[int]]:
    """Read origin_matrix.csv and return as list of lists."""
    import csv
    origin_path = MATRIX_DIR / "origin_matrix.csv"
    with origin_path.open() as f:
        return [[int(x) for x in row] for row in csv.reader(f)]


DELTAS = [0, 13, 19, 29, 41, 67, 73, 79, 91, 97, 103, 111, 119]
DELTAS_17 = [0, 13, 19, 29, 41, 67, 73, 79, 91, 97, 103, 111, 119, 127, 131, 137, 149]
RTL_WRAP_BASE_13 = [0, 270, 415, 337, 182, 301, 445, 70, 490, 65, 85, 75, 216]
WRAP_BASE = RTL_WRAP_BASE_13
WND_M13_K67 = [31, 32, 33, 33, 32, 22, 17, 12, 6, 3, 2, 0, 0]


def _family_matrix_dir(m: int) -> Path:
    return MATRIX_DIR / f"{m}x{m + BASE_PAYLOAD_COLS}"


def _read_family_triplet(m: int) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    rate_dir = _family_matrix_dir(m)
    shift_file = sorted((rate_dir / "matrix").glob("*_QC_H_*.txt"))[0]
    occ_file = sorted((rate_dir / "occupied_matrix").glob("*_occupied_*.txt"))[0]
    fade_file = sorted((rate_dir / "fade_matrix").glob("*_fade_*.txt"))[0]
    return read_int_matrix(shift_file), read_int_matrix(occ_file), read_int_matrix(fade_file)


def _crop_occ_fade(occ: list[list[int]], fade: list[list[int]], m: int,
                   payload_cols: int) -> tuple[list[list[int]], list[list[int]]]:
    full_cols = len(occ[0])
    parity_start = full_cols - m
    keep_cols = list(range(payload_cols)) + list(range(parity_start, full_cols))

    occ_crop = [[row[c] for c in keep_cols] for row in occ]
    fade_crop = [[row[c] for c in keep_cols] for row in fade]
    return occ_crop, fade_crop


def _rebuild_shifts_from_wrap_base(
    occ: list[list[int]],
    fade: list[list[int]],
    deltas: list[int],
    wrap_bases: list[int],
    z: int = Z,
) -> list[list[int]]:
    rows = len(occ)
    cols = len(occ[0])
    shifts = [[-1 for _ in range(cols)] for _ in range(rows)]

    for r in range(rows):
        delta = deltas[r]
        cursor = (wrap_bases[r] + delta) % z
        for c in range(cols):
            if occ[r][c] or fade[r][c]:
                shifts[r][c] = cursor
                cursor = (cursor + delta) % z
    return shifts


def _active_cols(occ_row: list[int], fade_row: list[int]) -> list[int]:
    return [idx for idx, (o, f) in enumerate(zip(occ_row, fade_row)) if o or f]


def _first_active_shift(row_shifts: list[int]) -> int:
    return next(v for v in row_shifts if v >= 0)


def _last_active_shift(row_shifts: list[int]) -> int:
    return next(v for v in reversed(row_shifts) if v >= 0)


def _first_parity_shift(row_shifts: list[int], occ_row: list[int], fade_row: list[int],
                        payload_cols: int) -> int:
    for col in range(payload_cols, len(row_shifts)):
        if occ_row[col] or fade_row[col]:
            return row_shifts[col]
    raise ValueError("No parity CPM found in row")


def _wrap_base_delta_from_last(last_element: int, delta: int, z: int = Z) -> int:
    if delta == 0:
        return 0
    for k in range(z):
        if (last_element + k * delta) % z == 0:
            return k
    raise ValueError(f"No wrap_base_delta found for last={last_element}, delta={delta}")


def _reconstructed_family_view(m: int, payload_cols: int,
                               wrap_bases: list[int]) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    _, occ_full, fade_full = _read_family_triplet(m)
    occ, fade = _crop_occ_fade(occ_full, fade_full, m, payload_cols)
    shifts = _rebuild_shifts_from_wrap_base(occ, fade, DELTAS_17[:m], wrap_bases[:m])
    return shifts, occ, fade


def fig4_shift_delta():
    origin = _read_origin_matrix()
    m, n = len(origin), len(origin[0])

    fig, axes = plt.subplots(3, 1, figsize=(18, 12))

    # Panel A: Delta chain formula verification S_r(n) vs origin_matrix
    ax = axes[0]
    sample_rows = [1, 5, 10]
    colors_rows = ["#2196F3", "#4CAF50", "#E91E63"]
    for row, color in zip(sample_rows, colors_rows):
        shifts = [origin[row][c] for c in range(n) if origin[row][c] != -1]
        T = len(shifts)
        # Formula predictions
        predicted = [(WRAP_BASE[row] + (i + 1) * DELTAS[row]) % Z for i in range(T)]
        ax.plot(range(T), shifts, "o", color=color, markersize=4, alpha=0.7,
                label=f"Row {row} origin (T={T})")
        ax.plot(range(T), predicted, "-", color=color, linewidth=1.5, alpha=0.5,
                label=f"Row {row} S_r(n) formula")

    ax.set_xlabel("Chain Position n")
    ax.set_ylabel("Shift Value (mod 512)")
    ax.set_title(
        r"A. Delta Chain Verification: $S_r(n) = (\mathrm{WRAP\_BASE}[r] + (n+1) \times \Delta[r])\ \mathrm{mod}\ 512$",
        fontsize=12, fontweight="bold",
    )
    ax.legend(fontsize=8, ncol=3)
    ax.set_ylim(-10, 520)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
    ax.axhline(Z, color="gray", linewidth=0.5, linestyle=":")

    # Panel B: Row 1 shift bar chart with delta arrows
    ax = axes[1]
    row = 1
    delta = DELTAS[row]
    shifts = [origin[row][c] for c in range(n) if origin[row][c] != -1]
    x = np.arange(len(shifts))
    ax.bar(x, shifts, color="#2196F3", alpha=0.7, edgecolor="white")
    for i in range(1, min(len(shifts), 20)):
        diff = (shifts[i] - shifts[i - 1]) % Z
        ax.annotate(
            "", xy=(i, shifts[i]), xytext=(i - 1, shifts[i - 1]),
            arrowprops=dict(arrowstyle="->", color="red", lw=1.5),
        )
        if i <= 8:
            label = f"+{delta}" if diff == delta else f"+{diff}(wrap)"
            ax.text(
                i - 0.5, max(shifts[i - 1], shifts[i]) + 15,
                label, ha="center", fontsize=7, color="red",
            )

    wb = WRAP_BASE[row]
    ax.set_xlabel("Chain Position n")
    ax.set_ylabel("Shift Value")
    ax.set_title(
        f"B. Row 1 Delta Chain: delta={delta}, WRAP_BASE={wb}, "
        f"S(0)=({wb}+{delta})%512={shifts[0]}",
        fontsize=10, fontweight="bold",
    )
    ax.set_ylim(-10, 560)

    # Panel C: Delta values per row
    ax = axes[2]
    row_indices = list(range(MAX_M_13))
    bars = ax.bar(
        row_indices, DELTAS,
        color=["#FF9800" if d > 0 else "#9E9E9E" for d in DELTAS],
        edgecolor="black",
    )
    ax.set_xlabel("Row Slot Index")
    ax.set_ylabel("Delta Value")
    ax.set_title(
        "C. Per-Row Delta Values (structural constants, shared by all (M,K) configs)",
        fontsize=12, fontweight="bold",
    )
    ax.set_xticks(row_indices)
    for i, d in enumerate(DELTAS):
        ax.text(i, d + 2, str(d), ha="center", fontsize=8, fontweight="bold")

    fig.suptitle(
        r"Figure 4: Delta Chain $S_r(n)$ — Verified Against origin_matrix.csv",
        fontsize=14, fontweight="bold",
    )
    fig.tight_layout()
    path = OUTPUT_DIR / "fig4_shift_delta.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 5: LUT Structure (per-slot 36-entry table)
# ---------------------------------------------------------------------------
def fig5_lut_structure():
    lut = load_lut_summary()
    min_m_lut = lut["min_m"]  # 5
    wbd = lut["wrap_base_delta"]  # 17 rows, each 52 entries (13 M-values x 4 K-values)
    k_cols = lut["lut_payload_cols"]  # [64,65,66,67]
    m_values = lut["m_values"]  # [5..17]

    # For display, show first 13 rows (slot 0-12) with M=5..13 scope
    # Reshape wbd[row] into 2D: (extra_rows, extra_userdata_cols) within M=5..13
    num_k = len(k_cols)
    num_m_13 = 9  # M=5..13

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    sample_slots = [0, 1, 4, 7, 10, 12]
    for ax, slot in zip(axes.flat, sample_slots):
        table = np.zeros((num_m_13, num_k), dtype=int)
        valid = np.ones((num_m_13, num_k), dtype=bool)
        for mi, m_val in enumerate(range(5, 14)):
            for ki, kval in enumerate(k_cols):
                idx = mi * num_k + ki
                if slot >= m_val:
                    valid[mi][ki] = False
                    table[mi][ki] = -1
                else:
                    table[mi][ki] = wbd[slot][idx]

        display = np.where(valid, table, np.nan)
        im = ax.imshow(display, aspect="auto", cmap="YlOrRd", interpolation="nearest",
                       vmin=0, vmax=max(1, np.nanmax(display)))

        for mi in range(num_m_13):
            for ki in range(num_k):
                if valid[mi][ki]:
                    ax.text(ki, mi, str(table[mi][ki]), ha="center", va="center",
                            fontsize=8, fontweight="bold")
                else:
                    ax.text(ki, mi, "-", ha="center", va="center",
                            fontsize=8, color="gray")

        ax.set_xticks(range(num_k))
        ax.set_xticklabels([f"K={k}" for k in k_cols], fontsize=8)
        ax.set_yticks(range(num_m_13))
        ax.set_yticklabels([f"M={m}" for m in range(5, 14)], fontsize=8)
        ax.set_title(f"Slot {slot}  (delta={lut['delta'][slot]})", fontsize=11, fontweight="bold")

        # Mark invalid region
        for mi in range(num_m_13):
            for ki in range(num_k):
                if not valid[mi][ki]:
                    ax.add_patch(plt.Rectangle((ki-0.5, mi-0.5), 1, 1,
                                               fill=True, facecolor="#E0E0E0",
                                               edgecolor="gray", linewidth=0.5))
                    ax.text(ki, mi, "-", ha="center", va="center",
                            fontsize=8, color="gray")

    fig.suptitle("Figure 5: Script-Derived Per-Slot first_step LUT (from generate_family_artifacts.py)\n"
                 "NOT the legacy DV wrap_base_deltaXX_adj_ord  |  Index = (M-5)*4 + (K-64), gray = slot absent",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    path = OUTPUT_DIR / "fig5_lut_structure.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 6: Pipeline Overview
# ---------------------------------------------------------------------------
def fig6_pipeline():
    fig, ax = plt.subplots(figsize=(18, 10))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def box(x, y, w, h, text, color="#E3F2FD", border="#1565C0", fontsize=9):
        rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                                        facecolor=color, edgecolor=border, linewidth=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", wrap=True)

    def arrow(x1, y1, x2, y2, text="", color="#1565C0"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=2))
        if text:
            mx, my = (x1+x2)/2, (y1+y2)/2
            ax.text(mx, my + 0.2, text, ha="center", fontsize=7, color=color, style="italic")

    # Input sources
    box(0.3, 8.5, 3.7, 1.0,
        "ibex_matrix/\n5x72 .. 17x84\nmax-K matrices + occupied/fade",
        "#FFF3E0", "#E65100")
    box(4.8, 8.5, 3.5, 1.0,
        "WrapBase source\nfixed RTL list or\nAutoBase search target",
        "#FFF3E0", "#E65100")

    # Stage 1
    box(0.5, 6.4, 4.0, 1.3,
        "Stage 1\n"
        "generate_family_artifacts.py\n"
        "copy_normalized()\nwrite_bm_schematic()\nwrite_assignments_svh()\nwrite_ens_files()",
        "#E3F2FD")
    arrow(2.15, 8.5, 2.15, 7.7)

    box(5.3, 6.4, 3.8, 1.3,
        "Stage 1 Outputs\nmatrix / occupied / fade\nENS / assignments SVH",
        "#E8F5E9", "#2E7D32")
    arrow(4.5, 7.05, 5.3, 7.05)

    # Stage 2a / 2b
    box(0.8, 4.0, 4.2, 1.45,
        "Stage 2a\n"
        "generate_family_artifacts_from_occ_fade_rtl.py\n"
        "read assignments SVH + fixed WrapBase\n"
        "rebuild shifts left->right",
        "#FFEBEE", "#C62828")
    arrow(7.2, 6.4, 3.0, 5.45, "assignments SVH")
    arrow(4.8, 8.5, 3.9, 5.45, "fixed WrapBase")

    box(5.8, 4.0, 4.2, 1.45,
        "Stage 2b\n"
        "generate_family_artifacts_from_occ_fade_rtl_autobase.py\n"
        "search WrapBase\n"
        "minimize max WrapBaseDelta",
        "#F3E5F5", "#6A1B9A")
    arrow(7.2, 6.4, 7.9, 5.45, "assignments SVH")
    arrow(6.6, 8.5, 7.9, 5.45, "candidate base")

    box(11.0, 6.25, 3.2, 1.45,
        "RDEC_schedule\nseparate path\nuses existing QC entries\nno semantic change here",
        "#E8F5E9", "#2E7D32")
    arrow(9.1, 7.05, 11.0, 7.05)

    # Outputs
    box(0.7, 1.8, 4.3, 1.35,
        "Fixed-Base Outputs\nwrap_base_delta / first_mask_shift\nSVH + JSON + CSV\nreconstructed matrices",
        "#E8F5E9", "#2E7D32")
    arrow(2.9, 4.0, 2.9, 3.15)

    box(5.8, 1.8, 4.4, 1.35,
        "AutoBase Outputs\nautobase LUT + CSV\nshift_dist / chosen WrapBase\nreconstructed matrices",
        "#E8F5E9", "#2E7D32")
    arrow(7.9, 4.0, 7.9, 3.15)

    box(12.0, 1.8, 4.8, 1.35,
        "Key semantics\nStage 1 writes bitmaps\nStage 2 rebuilds shifts\nleft to right from WrapBase",
        "#FFFDE7", "#F57F17", fontsize=8)

    # Final manifest
    box(3.4, 0.2, 6.5, 1.05,
        "Family artifact set = Stage 1 outputs + Stage 2 parameter tables + reconstructed matrix views",
        "#F3E5F5", "#6A1B9A")
    arrow(5.0, 1.8, 6.65, 1.25)
    arrow(8.0, 1.8, 6.65, 1.25)

    ax.set_title("Figure 6: Two-Stage Offline Pipeline (Stage 1 bitmaps, Stage 2 RTL-aligned parameters)\n"
                 "RDEC_schedule remains a parallel artifact path derived from existing QC entries",
                 fontsize=13, fontweight="bold", pad=20)

    path = OUTPUT_DIR / "fig6_pipeline.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 7: Wraparound Identity Verification
# ---------------------------------------------------------------------------
def fig7_wraparound_identity():
    origin = _read_origin_matrix()
    m, n = len(origin), len(origin[0])

    fig, axes = plt.subplots(2, 1, figsize=(16, 10))

    # Panel A: S_r(wnd) vs measured wraparound for all 13 rows
    ax = axes[0]
    rows_idx = list(range(m))
    wraparounds = []
    s_wnd_vals = []
    for r in range(m):
        shifts = [origin[r][c] for c in range(n) if origin[r][c] != -1]
        first, last = shifts[0], shifts[-1]
        wrap = (Z + first - last) % Z
        wraparounds.append(wrap)
        wnd = WND_M13_K67[r]
        s_wnd = (WRAP_BASE[r] + (wnd + 1) * DELTAS[r]) % Z
        s_wnd_vals.append(s_wnd)

    bar_width = 0.35
    x = np.arange(m)
    ax.bar(x - bar_width / 2, wraparounds, bar_width, label="Measured wraparound",
           color="#2196F3", edgecolor="black", linewidth=0.5)
    ax.bar(x + bar_width / 2, s_wnd_vals, bar_width, label=r"$S_r(\mathrm{wnd})$ formula",
           color="#FF9800", edgecolor="black", linewidth=0.5, alpha=0.8)
    for i in range(m):
        match = "=" if wraparounds[i] == s_wnd_vals[i] else "X"
        ax.text(i, max(wraparounds[i], s_wnd_vals[i]) + 10, match,
                ha="center", fontsize=9, fontweight="bold",
                color="green" if match == "=" else "red")

    ax.set_xlabel("Row Index")
    ax.set_ylabel("Value (mod 512)")
    ax.set_title(
        r"A. Wraparound Identity: $S_r(\mathrm{wnd}) = (\mathrm{WRAP\_BASE} + (\mathrm{wnd}+1) \times \Delta)\ \mathrm{mod}\ 512 = (Z + \mathrm{FE} - \mathrm{LE})\ \mathrm{mod}\ Z$",
        fontsize=11, fontweight="bold",
    )
    ax.set_xticks(rows_idx)
    ax.legend(fontsize=10)
    ax.set_ylim(0, 560)

    # Panel B: Decoder rotation = wrap_base + wnd * delta
    ax = axes[1]
    decoder_rot = [(WRAP_BASE[r] + WND_M13_K67[r] * DELTAS[r]) % Z for r in range(m)]
    expected_rot = [(Z + wraparounds[r] - DELTAS[r]) % Z for r in range(m)]

    ax.bar(x - bar_width / 2, decoder_rot, bar_width,
           label="WRAP_BASE + wnd*DELTA", color="#4CAF50", edgecolor="black", linewidth=0.5)
    ax.bar(x + bar_width / 2, expected_rot, bar_width,
           label="wraparound - DELTA", color="#E91E63", edgecolor="black", linewidth=0.5, alpha=0.8)

    for i in range(m):
        match = "=" if decoder_rot[i] == expected_rot[i] else "X"
        ax.text(i, max(decoder_rot[i], expected_rot[i]) + 10, match,
                ha="center", fontsize=9, fontweight="bold",
                color="green" if match == "=" else "red")

    ax.set_xlabel("Row Index")
    ax.set_ylabel("Decoder Rotation (mod 512)")
    ax.set_title(
        "B. Decoder Rotation = WRAP_BASE + wnd * DELTA (M=13, K=67)",
        fontsize=11, fontweight="bold",
    )
    ax.set_xticks(rows_idx)
    ax.legend(fontsize=10)
    ax.set_ylim(0, 560)

    fig.suptitle(
        "Figure 7: Wraparound Identity Verification (origin_matrix.csv, M=13, K=67)",
        fontsize=14, fontweight="bold",
    )
    fig.tight_layout()
    path = OUTPUT_DIR / "fig7_wraparound_identity.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 8: G/E Sub-Matrix Structure & RTL Two-Stage Derivation
# ---------------------------------------------------------------------------
def fig8_constant_derivation():
    g = 5  # G sub-matrix size

    def last_element_for_row(row: int) -> int:
        if row < g - 1:
            return (4 * DELTAS_17[row]) % Z
        if row == g - 1:
            return (3 * DELTAS_17[row]) % Z
        return 0

    fig = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1.2, 1, 1.2], hspace=0.35, wspace=0.3)

    # --- Panel A: Parity sub-matrix structure for M=10 (10x77) ---
    ax_a = fig.add_subplot(gs[0, :])
    m_demo = 10
    n_demo = m_demo + 67
    parity_cols = m_demo
    parity_start = n_demo - parity_cols  # col 67

    # Build parity block occupied pattern
    parity_mat = np.zeros((m_demo, parity_cols), dtype=float)
    parity_labels = [["" for _ in range(parity_cols)] for _ in range(m_demo)]

    # G sub-matrix: rows 0..4, rightmost 5 parity cols
    g_parity_start = parity_cols - g  # col offset 5 within parity region
    for i in range(g):
        g_ind = 0
        for j in range(g):
            real_col = g_parity_start + j
            k = g - 1 - j
            use = (k % g) != i
            if i == g - 1 and j == g - 1:
                use = False
            if use:
                g_ind += 1
                shift = (g_ind * DELTAS_17[i]) % Z
                parity_mat[i][real_col] = 3  # G block
                parity_labels[i][real_col] = str(shift)

    # E sub-matrix: rows 5..9, diagonal entries
    for i in range(g, m_demo):
        diag_col = parity_cols - 1 - i  # relative to parity start
        if 0 <= diag_col < parity_cols:
            parity_mat[i][diag_col] = 2  # E block
            parity_labels[i][diag_col] = "0"

    from matplotlib.colors import ListedColormap
    cmap_p = ListedColormap(["#F5F5F5", "#2196F3", "#FF9800", "#E91E63"])
    ax_a.imshow(parity_mat, cmap=cmap_p, aspect="auto", interpolation="nearest")

    for i in range(m_demo):
        for j in range(parity_cols):
            if parity_labels[i][j]:
                color = "white" if parity_mat[i][j] >= 2 else "black"
                ax_a.text(j, i, parity_labels[i][j], ha="center", va="center",
                          fontsize=8, fontweight="bold", color=color)

    # Zone boundary between E and G
    ax_a.axvline(g_parity_start - 0.5, color="black", linewidth=2, linestyle="--")

    # Zone labels
    mid_e = (g_parity_start - 1) / 2
    mid_g = g_parity_start + (g - 1) / 2
    ax_a.text(mid_e, -1.0, "E sub-matrix\n(diagonal, shift=0)",
              ha="center", fontsize=10, fontweight="bold", color="#FF9800")
    ax_a.text(mid_g, -1.0, "G sub-matrix\n(shift = g_ind × delta)",
              ha="center", fontsize=10, fontweight="bold", color="#E91E63")

    ax_a.set_xticks(range(parity_cols))
    ax_a.set_xticklabels([f"p{j}" for j in range(parity_cols)], fontsize=8)
    ax_a.set_yticks(range(m_demo))
    ax_a.set_yticklabels([f"row {i}" for i in range(m_demo)], fontsize=9)
    ax_a.set_title(f"A. Parity Region Structure (M={m_demo}, {parity_cols} parity cols) — "
                   "Shift values shown inside cells",
                   fontsize=12, fontweight="bold", pad=20)

    legend_elements = [
        mpatches.Patch(facecolor="#E91E63", label="G sub-matrix (rows 0-4)"),
        mpatches.Patch(facecolor="#FF9800", label="E sub-matrix (rows 5+, identity)"),
        mpatches.Patch(facecolor="#F5F5F5", edgecolor="gray", label="Empty"),
    ]
    ax_a.legend(handles=legend_elements, loc="lower right", fontsize=9, ncol=3)

    # --- Panel B: last_element per row (bar chart) ---
    ax_b = fig.add_subplot(gs[1, 0])
    num_rows = 17
    last_elements = [last_element_for_row(r) for r in range(num_rows)]
    colors_le = ["#E91E63" if r < g - 1 else "#FF5722" if r == g - 1 else "#9E9E9E"
                 for r in range(num_rows)]
    bars = ax_b.bar(range(num_rows), last_elements, color=colors_le, edgecolor="black", linewidth=0.5)
    for i, (bar, le) in enumerate(zip(bars, last_elements)):
        if le > 0:
            ax_b.text(i, le + 3, str(le), ha="center", fontsize=7, fontweight="bold")
        else:
            ax_b.text(i, le + 3, "0", ha="center", fontsize=7, color="gray")
    ax_b.set_xlabel("Row Index", fontsize=10)
    ax_b.set_ylabel("last_element Value", fontsize=10)
    ax_b.set_title("B. Fixed last_element[r] (Parity Anchor)", fontsize=11, fontweight="bold")
    ax_b.set_xticks(range(num_rows))

    le_legend = [
        mpatches.Patch(facecolor="#E91E63", label="G rows 0-3: 4×δ"),
        mpatches.Patch(facecolor="#FF5722", label="G row 4: 3×δ"),
        mpatches.Patch(facecolor="#9E9E9E", label="E rows 5+: 0"),
    ]
    ax_b.legend(handles=le_legend, fontsize=8, loc="upper left")

    # --- Panel C: Delta per row (extended to 17 rows) ---
    ax_c = fig.add_subplot(gs[1, 1])
    colors_d = ["#FF9800" if d > 0 else "#9E9E9E" for d in DELTAS_17]
    bars_d = ax_c.bar(range(num_rows), DELTAS_17, color=colors_d, edgecolor="black", linewidth=0.5)
    for i, (bar, d) in enumerate(zip(bars_d, DELTAS_17)):
        ax_c.text(i, d + 2, str(d), ha="center", fontsize=7, fontweight="bold")
    ax_c.set_xlabel("Row Index", fontsize=10)
    ax_c.set_ylabel("Delta Value", fontsize=10)
    ax_c.set_title("C. Fixed delta[r] (Step Size Constants, 17 Rows)", fontsize=11, fontweight="bold")
    ax_c.set_xticks(range(num_rows))

    # --- Panel D: RTL two-stage derivation pipeline ---
    ax_d = fig.add_subplot(gs[2, :])
    ax_d.set_xlim(0, 20)
    ax_d.set_ylim(0, 4)
    ax_d.axis("off")

    def draw_box(x, y, w, h, text, color="#E3F2FD", border="#1565C0", fs=9):
        rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                        facecolor=color, edgecolor=border, linewidth=2)
        ax_d.add_patch(rect)
        ax_d.text(x + w/2, y + h/2, text, ha="center", va="center",
                  fontsize=fs, fontweight="bold")

    def draw_arrow(x1, y1, x2, y2, text="", color="#1565C0"):
        ax_d.annotate("", xy=(x2, y2), xytext=(x1, y1),
                      arrowprops=dict(arrowstyle="-|>", color=color, lw=2))
        if text:
            mx, my = (x1+x2)/2, (y1+y2)/2
            ax_d.text(mx, my + 0.15, text, ha="center", fontsize=7, color=color, style="italic")

    # Stage-1 inputs
    draw_box(0.2, 2.5, 3.9, 1.2,
             "Stage 1: assignments SVH\n(occupied/fade bitmaps)",
             "#FFF3E0", "#E65100")
    draw_box(4.8, 2.5, 3.7, 1.2,
             "WrapBase[r]\n(fixed input or\nAutoBase search)",
             "#FFF3E0", "#E65100")

    # Process
    draw_box(9.4, 2.5, 4.3, 1.2,
             "rebuild_shifts_from_wrap_base()\n"
             "cursor = WrapBase[r] + delta[r]\n"
             "cursor += delta[r]\n"
             "(left to right)",
             "#E3F2FD", "#1565C0")
    draw_arrow(4.1, 3.1, 9.4, 3.1, "occupied/fade bitmaps")
    draw_arrow(8.5, 3.1, 9.4, 3.1, "WrapBase + delta")

    # Outputs
    draw_box(0.4, 0.3, 3.2, 1.2,
             "WrapBase\n(transparent input,\nnot derived)",
             "#E8F5E9", "#2E7D32")
    draw_box(4.3, 0.3, 4.0, 1.2,
             "WrapBaseDelta\n= min k: LE+k·δ≡0\nper (row, M, K)",
             "#E8F5E9", "#2E7D32")
    draw_box(9.0, 0.3, 4.0, 1.2,
             "FirstMaskShift\n= (Z+FE-parity_1st)%Z\nper (M, K)",
             "#E8F5E9", "#2E7D32")
    draw_box(13.8, 0.3, 5.0, 1.2,
             "RTL Aligned\nmatches ldpc_codec.cpp\nshift semantics",
             "#FCE4EC", "#880E4F")

    draw_box(14.4, 2.5, 4.1, 1.2,
             "SVH + JSON + CSV +\nreconstructed matrices",
             "#E8F5E9", "#2E7D32")

    draw_arrow(11.55, 2.5, 2.0, 1.5, "FE, LE")
    draw_arrow(11.55, 2.5, 6.3, 1.5)
    draw_arrow(11.55, 2.5, 11.0, 1.5)
    draw_arrow(13.7, 3.1, 14.4, 3.1)

    ax_d.set_title("D. RTL Two-Stage Parameter Derivation Pipeline",
                   fontsize=12, fontweight="bold")

    fig.suptitle("Figure 8: G/E Sub-Matrix Structure & RTL Left-to-Right Shift Reconstruction\n"
                 "Panel D follows the current RTL semantics: assignments SVH + WrapBase -> rebuild shifts -> WrapBaseDelta / FirstMaskShift",
                 fontsize=14, fontweight="bold", y=0.99)

    path = OUTPUT_DIR / "fig8_constant_derivation.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# Figure 9: WrapBaseDelta Ring-Track Intuition
# ---------------------------------------------------------------------------
def fig9_wrap_base_delta_ring():
    fig = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(2, 2, hspace=0.32, wspace=0.22)

    # --- Panel A: ring-track intuition on a single row ---
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_aspect("equal")
    ax_a.axis("off")

    demo_row = 1
    demo_m = 13
    demo_k = 67
    shifts_13x67, occ_13x67, fade_13x67 = _reconstructed_family_view(
        demo_m, demo_k, RTL_WRAP_BASE_13[:demo_m]
    )
    demo_row_shifts = [v for v in shifts_13x67[demo_row] if v >= 0]
    wrap_base = RTL_WRAP_BASE_13[demo_row]
    delta = DELTAS_17[demo_row]
    first_element = demo_row_shifts[0]
    last_element = demo_row_shifts[-1]
    wrap_base_delta = _wrap_base_delta_from_last(last_element, delta)

    theta = np.linspace(0, 2 * np.pi, 512)
    ax_a.plot(np.cos(theta), np.sin(theta), color="#37474F", linewidth=2)

    def point_on_ring(value: int, radius: float = 1.0) -> tuple[float, float]:
        ang = np.pi / 2 - 2 * np.pi * value / Z
        return radius * np.cos(ang), radius * np.sin(ang)

    x0, y0 = point_on_ring(0)
    ax_a.scatter([x0], [y0], s=90, color="#000000", zorder=5)
    ax_a.text(x0, y0 + 0.12, "0", ha="center", fontsize=10, fontweight="bold")

    xb, yb = point_on_ring(wrap_base)
    ax_a.scatter([xb], [yb], s=70, color="#FFB300", zorder=5)
    ax_a.text(xb - 0.16, yb - 0.08, f"WrapBase={wrap_base}", fontsize=9,
              color="#E65100", fontweight="bold")

    xf, yf = point_on_ring(first_element)
    ax_a.scatter([xf], [yf], s=85, color="#1E88E5", zorder=5)
    ax_a.text(xf + 0.08, yf + 0.08, f"FE={first_element}", fontsize=9,
              color="#1565C0", fontweight="bold")

    sample_points = demo_row_shifts[:4]
    for idx, value in enumerate(sample_points[1:], start=1):
        xp, yp = point_on_ring(value)
        ax_a.scatter([xp], [yp], s=55, color="#64B5F6", zorder=4)
        ax_a.text(xp + 0.05, yp, str(value), fontsize=8, color="#1565C0")
    ax_a.text(0, 0, "active CPMs\nstep by +delta", ha="center", va="center",
              fontsize=10, color="#546E7A")

    xl, yl = point_on_ring(last_element)
    ax_a.scatter([xl], [yl], s=85, color="#E53935", zorder=5)
    ax_a.text(xl + 0.08, yl - 0.12, f"LE={last_element}", fontsize=9,
              color="#B71C1C", fontweight="bold")

    ax_a.annotate(
        "",
        xy=(x0, y0),
        xytext=(xl, yl),
        arrowprops=dict(arrowstyle="-|>", color="#B71C1C", lw=2,
                        linestyle=(0, (4, 3))),
    )
    ax_a.text(0.0, -1.28,
              f"Remaining arc = WrapBaseDelta × delta = {wrap_base_delta} × {delta}",
              ha="center", fontsize=10, color="#B71C1C", fontweight="bold")
    ax_a.set_title(
        f"A. Ring Track Example (row {demo_row}, delta={delta}, M={demo_m}, K={demo_k})",
        fontsize=12, fontweight="bold"
    )

    # --- Panel B: K variation vs WrapBaseDelta ---
    ax_b = fig.add_subplot(gs[0, 1])
    demo_rows = [1, 2, 4, 7, 10]
    colors = ["#1E88E5", "#43A047", "#E53935", "#8E24AA", "#FB8C00"]
    k_values = [64, 65, 66, 67]

    for row, color in zip(demo_rows, colors):
        y_vals = []
        for k in k_values:
            shifts, _, _ = _reconstructed_family_view(13, k, RTL_WRAP_BASE_13)
            last_shift = _last_active_shift(shifts[row])
            y_vals.append(_wrap_base_delta_from_last(last_shift, DELTAS_17[row]))
        ax_b.plot(k_values, y_vals, marker="o", linewidth=2, color=color,
                  label=f"row {row} (δ={DELTAS_17[row]})")

    ax_b.set_xlabel("K (payload columns)")
    ax_b.set_ylabel("WrapBaseDelta")
    ax_b.set_title("B. Smaller K Usually Pushes WrapBaseDelta Up",
                   fontsize=12, fontweight="bold")
    ax_b.set_xticks(k_values)
    ax_b.grid(alpha=0.25, linestyle="--")
    ax_b.legend(fontsize=8, ncol=2)

    # --- Panel C: AutoBase search curve for one row ---
    ax_c = fig.add_subplot(gs[1, 0])
    search_row = 4
    search_delta = DELTAS_17[search_row]
    candidate_wrap_bases = np.arange(Z)
    max_wrap_base_delta = []

    occ_full, fade_full = _read_family_triplet(13)[1:]
    for candidate in candidate_wrap_bases:
        per_k = []
        for k in k_values:
            occ_k, fade_k = _crop_occ_fade(occ_full, fade_full, 13, k)
            row_occ = occ_k[search_row]
            row_fade = fade_k[search_row]
            row_shifts = _rebuild_shifts_from_wrap_base(
                [row_occ], [row_fade], [search_delta], [int(candidate)]
            )[0]
            per_k.append(_wrap_base_delta_from_last(_last_active_shift(row_shifts), search_delta))
        max_wrap_base_delta.append(max(per_k))

    max_wrap_base_delta = np.array(max_wrap_base_delta)
    best_base = int(candidate_wrap_bases[np.argmin(max_wrap_base_delta)])
    best_val = int(max_wrap_base_delta.min())
    official_base = RTL_WRAP_BASE_13[search_row]
    official_val = int(max_wrap_base_delta[official_base])

    ax_c.plot(candidate_wrap_bases, max_wrap_base_delta, color="#3949AB", linewidth=1.8)
    ax_c.scatter([best_base], [best_val], color="#E53935", s=70, zorder=5,
                 label=f"best={best_base}, maxΔ={best_val}")
    ax_c.scatter([official_base], [official_val], color="#FFB300", s=70, zorder=5,
                 label=f"RTL WrapBase={official_base}, maxΔ={official_val}")
    ax_c.axvline(official_base, color="#FFB300", linewidth=1.2, linestyle="--", alpha=0.7)
    ax_c.set_xlabel("Candidate WrapBase (0..511)")
    ax_c.set_ylabel("max WrapBaseDelta over K=64..67")
    ax_c.set_title(
        f"C. AutoBase Search Curve (row {search_row}, delta={search_delta})",
        fontsize=12, fontweight="bold"
    )
    ax_c.grid(alpha=0.25, linestyle="--")
    ax_c.legend(fontsize=8)

    # --- Panel D: first-mask-shift gate distance ---
    ax_d = fig.add_subplot(gs[1, 1])
    demo_m2 = 5
    demo_k2 = 66
    wrap_bases_5 = RTL_WRAP_BASE_13[:demo_m2]
    shifts_5x66, occ_5x66, fade_5x66 = _reconstructed_family_view(demo_m2, demo_k2, wrap_bases_5)
    last_row = demo_m2 - 1
    row_shifts = shifts_5x66[last_row]
    row_occ = occ_5x66[last_row]
    row_fade = fade_5x66[last_row]
    payload_cols = demo_k2

    payload_pts = [(col, shift) for col, shift in enumerate(row_shifts[:payload_cols]) if shift >= 0]
    parity_pts = [
        (col, shift)
        for col, shift in enumerate(row_shifts[payload_cols:], start=payload_cols)
        if shift >= 0
    ]
    fe = _first_active_shift(row_shifts)
    first_parity = _first_parity_shift(row_shifts, row_occ, row_fade, payload_cols)
    first_mask_shift = (Z + fe - first_parity) % Z

    ax_d.axhline(0, color="#546E7A", linewidth=1.5)
    ax_d.scatter([shift for _, shift in payload_pts], [0] * len(payload_pts),
                 color="#1E88E5", s=60, label="payload active CPMs")
    ax_d.scatter([shift for _, shift in parity_pts], [0] * len(parity_pts),
                 color="#43A047", s=70, marker="s", label="parity active CPMs")
    ax_d.scatter([fe], [0], color="#E53935", s=100, zorder=6)
    ax_d.scatter([first_parity], [0], color="#6A1B9A", s=110, marker="D", zorder=6)
    ax_d.annotate(f"FE={fe}", xy=(fe, 0), xytext=(fe, 0.22),
                  ha="center", fontsize=9, fontweight="bold", color="#B71C1C")
    ax_d.annotate(f"first parity={first_parity}", xy=(first_parity, 0), xytext=(first_parity, -0.28),
                  ha="center", fontsize=9, fontweight="bold", color="#6A1B9A")
    ax_d.annotate(
        "",
        xy=(first_parity, 0.12),
        xytext=(fe, 0.12),
        arrowprops=dict(arrowstyle="<->", color="#455A64", lw=2),
    )
    ax_d.text((fe + first_parity) / 2.0, 0.18,
              f"FirstMaskShift = (512 + {fe} - {first_parity}) % 512 = {first_mask_shift}",
              ha="center", fontsize=9, color="#37474F", fontweight="bold")
    ax_d.set_xlim(-10, 522)
    ax_d.set_ylim(-0.45, 0.45)
    ax_d.set_yticks([])
    ax_d.set_xlabel("Shift value on the rebuilt row-local delta chain")
    ax_d.set_title(
        f"D. FirstMaskShift as a Gate Distance (M={demo_m2}, K={demo_k2}, last row)",
        fontsize=12, fontweight="bold"
    )
    ax_d.grid(axis="x", alpha=0.2, linestyle="--")
    ax_d.legend(fontsize=8, loc="upper left")

    fig.suptitle(
        'Figure 9: WrapBaseDelta "Ring Track" Intuition & AutoBase Search',
        fontsize=14, fontweight="bold"
    )
    path = OUTPUT_DIR / "fig9_wrap_base_delta_ring.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    print("Generating IBEX Matrix Visualizations...")
    fig1_matrix_layout()
    fig2_5x5_block()
    fig3_trimming()
    fig4_shift_delta()
    fig5_lut_structure()
    fig6_pipeline()
    fig7_wraparound_identity()
    fig8_constant_derivation()
    fig9_wrap_base_delta_ring()
    print(f"\nAll figures saved to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
