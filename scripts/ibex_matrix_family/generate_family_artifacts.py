#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
# IBEX QC-LDPC fixed design constants (from GenLDPC_T4/src/gen_ldpc.c)
# ---------------------------------------------------------------------------
IBEX_DELTA = (0, 13, 19, 29, 41, 67, 73, 79, 91, 97, 103, 111, 119, 127, 131, 137, 149)
IBEX_G_SIZE = 5
IBEX_Z = 512


def ibex_last_element(row: int, delta: int, z: int = IBEX_Z) -> int:
    """Compute the fixed last_element (rightmost active parity CPM shift) for a row.

    Rows 0..3 (G sub-matrix): 4 * delta  (4 active parity entries)
    Row 4 (G last row):       3 * delta  (3 entries, one removed for invertibility)
    Rows 5+ (E sub-matrix):   0          (identity shift on diagonal)
    """
    if row < IBEX_G_SIZE - 1:
        return (4 * delta) % z
    if row == IBEX_G_SIZE - 1:
        return (3 * delta) % z
    return 0


def ibex_last_element_per_row(z: int = IBEX_Z) -> list[int]:
    """Return fixed last_element for all 17 rows."""
    return [ibex_last_element(row, IBEX_DELTA[row], z) for row in range(len(IBEX_DELTA))]


RATE_DIR_RE = re.compile(r"^(?P<m>\d+)x(?P<n>\d+)$")
FLAT_MATRIX_RE = re.compile(
    r"^LDPC_(?P<m>\d+)x(?P<n>\d+)ex(?P<z>\d+)_w(?P<wt>\d+)_dense(?P<dense>\d+)_QC_H_\d+(?:_\d+)?\.txt$"
)
FLAT_OCCUPIED_RE = re.compile(
    r"^LDPC_(?P<m>\d+)x(?P<n>\d+)ex(?P<z>\d+)_w(?P<wt>\d+)_dense(?P<dense>\d+)_occupied_\d+(?:_\d+)?\.txt$"
)


@dataclass
class Entry:
    row: int
    col: int
    shift: int


@dataclass
class RateMatrix:
    m: int
    n: int
    payload_cols: int
    z: int
    parity_bytes: int
    prefix: str
    occupied_path: Path
    fade_path: Path
    occupied: list[list[int]]
    fade: list[list[int]]
    matrix_path: Path | None = None
    shifts: list[list[int]] | None = None

    @property
    def extra_user_cols(self) -> int:
        return max(0, self.payload_cols - 64)

    @property
    def extra_bits_of_parity(self) -> int:
        bytes_per_col = self.z // 8
        unused_bytes = self.m * bytes_per_col - self.parity_bytes
        if unused_bytes < 0:
            raise ValueError(
                f"Invalid parity bytes for {self.m}x{self.n}: {self.parity_bytes}"
            )
        if unused_bytes == 0:
            return 0
        return (bytes_per_col - unused_bytes) * 8


@dataclass(frozen=True)
class LutConfig:
    m: int
    payload_cols: int
    rate: RateMatrix

    @property
    def n(self) -> int:
        return self.m + self.payload_cols

    @property
    def extra_user_cols(self) -> int:
        return max(0, self.payload_cols - 64)

    @property
    def label(self) -> str:
        return f"M{self.m}_K{self.payload_cols}"


def read_int_matrix(path: Path) -> list[list[int]]:
    rows: list[list[int]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            rows.append([int(token) for token in line.split()])
    return rows


def pick_matrix_file(directory: Path, pattern: str) -> Path:
    candidates = sorted(directory.glob(pattern))
    if not candidates:
        raise FileNotFoundError(f"No files matching {pattern} under {directory}")
    if len(candidates) == 1:
        return candidates[0]

    preferred_suffixes = (
        "_1_1.txt",
        "_1.txt",
        "_0.txt",
    )
    for suffix in preferred_suffixes:
        for candidate in candidates:
            if candidate.name.endswith(suffix):
                return candidate
    return candidates[0]


def infer_prefix(path: Path, kind: str) -> str:
    patterns = {
        "matrix": r"^(?P<prefix>.+)_QC_H_\d+(?:_\d+)?\.txt$",
        "occupied": r"^(?P<prefix>.+)_occupied_\d+(?:_\d+)?\.txt$",
        "fade": r"^(?P<prefix>.+)_fade_\d+(?:_\d+)?\.txt$",
    }
    match = re.match(patterns[kind], path.name)
    if not match:
        raise ValueError(f"Cannot infer prefix from {path}")
    return match.group("prefix")


def normalize_output_name(prefix: str, kind: str) -> str:
    if kind == "matrix":
        return f"{prefix}_QC_H_1.txt"
    if kind == "occupied":
        return f"{prefix}_occupied_1.txt"
    if kind == "fade":
        return f"{prefix}_fade_1.txt"
    raise ValueError(f"Unsupported kind {kind}")


def copy_normalized(src: Path, dst_dir: Path, filename: str) -> Path:
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst_path = dst_dir / filename
    shutil.copyfile(src, dst_path)
    return dst_path


def parity_bytes_for_rate(m: int, z: int, m5_parity_bytes: int) -> int:
    if m == 5:
        return m5_parity_bytes
    return m * (z // 8) - 1


def discover_rates(
    input_root: Path, payload_cols: int, z: int, m5_parity_bytes: int,
    skip_shifts: bool = False,
) -> list[RateMatrix]:
    rate_dirs = [child for child in sorted(input_root.iterdir()) if child.is_dir() and RATE_DIR_RE.match(child.name)]
    if rate_dirs:
        return discover_rates_from_nested_root(input_root, payload_cols, z, m5_parity_bytes, skip_shifts=skip_shifts)

    matrix_dir = input_root / "matrix"
    occupied_dir = input_root / "occupied_matrix"
    fade_dir = input_root / "fade_matrix"
    if skip_shifts and occupied_dir.is_dir() and fade_dir.is_dir():
        return discover_rates_from_flat_root(input_root, payload_cols, z, m5_parity_bytes, skip_shifts=True)
    if matrix_dir.is_dir() and occupied_dir.is_dir() and fade_dir.is_dir():
        return discover_rates_from_flat_root(input_root, payload_cols, z, m5_parity_bytes, skip_shifts=skip_shifts)

    raise FileNotFoundError(
        f"Input root {input_root} is neither a nested per-rate tree nor a flat matrix/fade/occupied tree"
    )


def discover_rates_from_nested_root(
    input_root: Path, payload_cols: int, z: int, m5_parity_bytes: int,
    skip_shifts: bool = False,
) -> list[RateMatrix]:
    rates: list[RateMatrix] = []
    for child in sorted(input_root.iterdir()):
        if not child.is_dir():
            continue
        match = RATE_DIR_RE.match(child.name)
        if not match:
            continue
        m = int(match.group("m"))
        n = int(match.group("n"))
        if n - m != payload_cols:
            continue

        occupied_path = pick_matrix_file(child / "occupied_matrix", "*_occupied_*.txt")
        fade_path = pick_matrix_file(child / "fade_matrix", "*_fade_*.txt")

        if skip_shifts:
            prefix = infer_prefix(occupied_path, "occupied")
            matrix_path = None
            shifts = None
        else:
            matrix_path = pick_matrix_file(child / "matrix", "*_QC_H_*.txt")
            prefix = infer_prefix(matrix_path, "matrix")
            shifts = read_int_matrix(matrix_path)

        occupied = read_int_matrix(occupied_path)
        fade = read_int_matrix(fade_path)
        if len(occupied) != m or len(fade) != m:
            raise ValueError(f"Matrix row count mismatch under {child}")
        if shifts is not None:
            if len(shifts) != m:
                raise ValueError(f"Shift matrix row count mismatch under {child}")
            if any(len(row) != n for row in shifts):
                raise ValueError(f"Shift matrix column count mismatch in {matrix_path}")
        if any(len(row) != n for row in occupied):
            raise ValueError(f"Occupied matrix column count mismatch in {occupied_path}")
        if any(len(row) != n for row in fade):
            raise ValueError(f"Fade matrix column count mismatch in {fade_path}")

        rates.append(
            RateMatrix(
                m=m,
                n=n,
                payload_cols=payload_cols,
                z=z,
                parity_bytes=parity_bytes_for_rate(m, z, m5_parity_bytes),
                prefix=prefix,
                occupied_path=occupied_path,
                fade_path=fade_path,
                occupied=occupied,
                fade=fade,
                matrix_path=matrix_path,
                shifts=shifts,
            )
        )

    if not rates:
        raise FileNotFoundError(
            f"No rate directories matching payload_cols={payload_cols} found under {input_root}"
        )
    return sorted(rates, key=lambda rate: rate.m)


def discover_rates_from_flat_root(
    input_root: Path, payload_cols: int, z: int, m5_parity_bytes: int,
    skip_shifts: bool = False,
) -> list[RateMatrix]:
    rates: list[RateMatrix] = []
    matrix_dir = input_root / "matrix"
    occupied_dir = input_root / "occupied_matrix"
    fade_dir = input_root / "fade_matrix"

    # When skip_shifts and matrix_dir is missing, discover rates from occupied files.
    if skip_shifts and not matrix_dir.is_dir():
        for occ_path in sorted(occupied_dir.glob("*_occupied_*.txt")):
            match = FLAT_OCCUPIED_RE.match(occ_path.name)
            if not match:
                continue
            m = int(match.group("m"))
            n = int(match.group("n"))
            matrix_z = int(match.group("z"))
            if matrix_z != z or n - m != payload_cols:
                continue
            prefix = infer_prefix(occ_path, "occupied")
            fade_path = fade_dir / normalize_output_name(prefix, "fade")
            if not fade_path.exists():
                fade_path = pick_matrix_file(fade_dir, f"{prefix}_fade_*.txt")
            occupied = read_int_matrix(occ_path)
            fade = read_int_matrix(fade_path)
            if len(occupied) != m or len(fade) != m:
                raise ValueError(f"Matrix row count mismatch for prefix {prefix}")
            if any(len(row) != n for row in occupied):
                raise ValueError(f"Occupied matrix column count mismatch in {occ_path}")
            if any(len(row) != n for row in fade):
                raise ValueError(f"Fade matrix column count mismatch in {fade_path}")
            rates.append(
                RateMatrix(
                    m=m, n=n, payload_cols=payload_cols, z=z,
                    parity_bytes=parity_bytes_for_rate(m, z, m5_parity_bytes),
                    prefix=prefix,
                    occupied_path=occ_path, fade_path=fade_path,
                    occupied=occupied, fade=fade,
                )
            )
        if not rates:
            raise FileNotFoundError(
                f"No occupied files matching payload_cols={payload_cols} found under {occupied_dir}"
            )
        return sorted(rates, key=lambda rate: rate.m)

    # Default path: discover from matrix file names.
    for matrix_path in sorted(matrix_dir.glob("*_QC_H_*.txt")):
        match = FLAT_MATRIX_RE.match(matrix_path.name)
        if not match:
            continue
        m = int(match.group("m"))
        n = int(match.group("n"))
        matrix_z = int(match.group("z"))
        if matrix_z != z:
            continue
        if n - m != payload_cols:
            continue

        prefix = infer_prefix(matrix_path, "matrix")
        occupied_path = occupied_dir / normalize_output_name(prefix, "occupied")
        fade_path = fade_dir / normalize_output_name(prefix, "fade")
        if not occupied_path.exists():
            occupied_path = pick_matrix_file(occupied_dir, f"{prefix}_occupied_*.txt")
        if not fade_path.exists():
            fade_path = pick_matrix_file(fade_dir, f"{prefix}_fade_*.txt")

        if skip_shifts:
            shifts = None
        else:
            shifts = read_int_matrix(matrix_path)

        occupied = read_int_matrix(occupied_path)
        fade = read_int_matrix(fade_path)
        if len(occupied) != m or len(fade) != m:
            raise ValueError(f"Matrix row count mismatch for prefix {prefix}")
        if shifts is not None:
            if len(shifts) != m:
                raise ValueError(f"Shift matrix row count mismatch for prefix {prefix}")
            if any(len(row) != n for row in shifts):
                raise ValueError(f"Shift matrix column count mismatch in {matrix_path}")
        if any(len(row) != n for row in occupied):
            raise ValueError(f"Occupied matrix column count mismatch in {occupied_path}")
        if any(len(row) != n for row in fade):
            raise ValueError(f"Fade matrix column count mismatch in {fade_path}")

        rates.append(
            RateMatrix(
                m=m, n=n, payload_cols=payload_cols, z=z,
                parity_bytes=parity_bytes_for_rate(m, z, m5_parity_bytes),
                prefix=prefix,
                occupied_path=occupied_path, fade_path=fade_path,
                occupied=occupied, fade=fade,
                matrix_path=matrix_path if not skip_shifts else None,
                shifts=shifts,
            )
        )

    if not rates:
        raise FileNotFoundError(
            f"No flat matrix files matching payload_cols={payload_cols} found under {input_root}"
        )
    return sorted(rates, key=lambda rate: rate.m)


def build_qc_entries(rate: RateMatrix) -> tuple[list[list[Entry]], list[list[Entry]]]:
    row_entries: list[list[Entry]] = [[] for _ in range(rate.m)]
    col_entries: list[list[Entry]] = [[] for _ in range(rate.n)]
    use_all_nonzero = rate.extra_bits_of_parity > 0

    for row in range(rate.m):
        for col in range(rate.n):
            shift = rate.shifts[row][col]
            include = shift >= 0 if use_all_nonzero else rate.occupied[row][col] == 1
            if not include:
                continue
            entry = Entry(row=row, col=col, shift=shift)
            row_entries[row].append(entry)
            col_entries[col].append(entry)

    return row_entries, col_entries


def prev_in_col(entry: Entry, col_entries: list[list[Entry]]) -> Entry:
    entries = col_entries[entry.col]
    index = next(i for i, candidate in enumerate(entries) if candidate.row == entry.row)
    if index == 0:
        return entries[-1]
    return entries[index - 1]


def write_ens_files(rate: RateMatrix, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    file_stem = f"{rate.m}x{rate.n}ex{rate.z}_w4"
    enc1_path = out_dir / f"enc1_sched_{file_stem}.txt"
    alias_path = out_dir / f"ENS_{file_stem}.txt"

    lines: list[str] = []
    for col in range(rate.n):
        occ_shift: list[int] = []
        occ_row: list[int] = []
        fade_rows: list[int] = []
        fade_shift = 0x1FF

        for row in range(rate.m):
            if rate.occupied[row][col] == 1 and rate.shifts[row][col] >= 0:
                occ_shift.append(rate.shifts[row][col])
                occ_row.append(row)
            if rate.fade[row][col] == 1:
                fade_rows.append(row)
                fade_shift = rate.shifts[row][col]

        if len(occ_shift) > 4:
            raise ValueError(
                f"Column {col} in {rate.m}x{rate.n} has more than 4 occupied entries"
            )

        while len(occ_shift) < 4:
            occ_shift.append(0x1FF)
            occ_row.append(0x1F)

        fade_row = fade_rows[-1] if fade_rows else 0x1F
        word = 0
        word |= occ_shift[0] & 0x1FF
        word |= (occ_shift[1] & 0x1FF) << 9
        word |= (occ_shift[2] & 0x1FF) << 18
        word |= (occ_shift[3] & 0x1FF) << 27
        word |= (fade_shift & 0x1FF) << 36
        word |= (occ_row[0] & 0x1F) << 45
        word |= (occ_row[1] & 0x1F) << 50
        word |= (occ_row[2] & 0x1F) << 55
        word |= (occ_row[3] & 0x1F) << 60
        word |= (fade_row & 0x1F) << 65
        lines.append(f"{word:018X}")

    for output_path in (enc1_path, alias_path):
        output_path.write_text("\n".join(lines) + "\n")
    return [enc1_path, alias_path]


def write_bm_schematic(rate: RateMatrix, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / f"bm_schematic_{rate.m}x{rate.n}ex{rate.z}_w4.txt"
    base_payload_cols = min(rate.payload_cols, 64)

    def is_nonzero(row: int, col: int) -> bool:
        if rate.extra_bits_of_parity > 0:
            return rate.shifts[row][col] >= 0
        return rate.occupied[row][col] == 1

    lines = [
        f"# Base-matrix schematic (bm_m={rate.m}, bm_n={rate.n}, Z={rate.z})",
        "# Legend: 1=non-zero CPM, 0=zero CPM, X=zero CPM in extra-userdata col (placeholder/skip)",
        "# Groups: base_userdata_cols|extra_userdata_cols|parity_cols",
        f"# base_userdata_cols={base_payload_cols}, extra_userdata_cols={rate.payload_cols - base_payload_cols}, parity_cols={rate.m}",
        f"# NOTE: extra_bits_of_userdata=0, extra_bits_of_parity={rate.extra_bits_of_parity}",
        "",
    ]

    for row in range(rate.m):
        chars: list[str] = []
        for col in range(base_payload_cols):
            chars.append("1" if is_nonzero(row, col) else "0")
        chars.append("|")
        for col in range(base_payload_cols, rate.payload_cols):
            chars.append("1" if is_nonzero(row, col) else "X")
        chars.append("|")
        for col in range(rate.payload_cols, rate.n):
            chars.append("1" if is_nonzero(row, col) else "0")
        lines.append(f"ROW {row:2d}: {''.join(chars)}")

    lines.extend(
        [
            "",
            f"# Fade matrix (bm_m={rate.m}, bm_n={rate.n}, Z={rate.z})",
            "# Legend: 1=fade CPM, .=non-fade",
            "# Groups: base_userdata_cols|extra_userdata_cols|parity_cols",
            "",
        ]
    )
    for row in range(rate.m):
        chars = []
        for col in range(base_payload_cols):
            chars.append("1" if rate.fade[row][col] == 1 else ".")
        chars.append("|")
        for col in range(base_payload_cols, rate.payload_cols):
            chars.append("1" if rate.fade[row][col] == 1 else ".")
        chars.append("|")
        for col in range(rate.payload_cols, rate.n):
            chars.append("1" if rate.fade[row][col] == 1 else ".")
        lines.append(f"ROW {row:2d}: {''.join(chars)}")

    lines.extend(
        [
            "",
            "# Base-matrix schematic (fade shown as '*')",
            "# Legend: 1=non-zero CPM, 0=zero CPM, X=zero CPM in extra-userdata col (placeholder/skip), *=fade CPM",
            "# Groups: base_userdata_cols|extra_userdata_cols|parity_cols",
            "",
        ]
    )
    for row in range(rate.m):
        chars = []
        for col in range(base_payload_cols):
            chars.append("*" if rate.fade[row][col] == 1 else ("1" if is_nonzero(row, col) else "0"))
        chars.append("|")
        for col in range(base_payload_cols, rate.payload_cols):
            if rate.fade[row][col] == 1:
                chars.append("*")
            else:
                chars.append("1" if is_nonzero(row, col) else "X")
        chars.append("|")
        for col in range(rate.payload_cols, rate.n):
            chars.append("*" if rate.fade[row][col] == 1 else ("1" if is_nonzero(row, col) else "0"))
        lines.append(f"ROW {row:2d}: {''.join(chars)}")

    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def pack_sched_word(
    rate: RateMatrix,
    row: int,
    entry: Entry,
    pre_entry: Entry,
    emitted_real_count: int,
    row_weight: int,
) -> int:
    shift_delta = (entry.shift - pre_entry.shift + rate.z) % rate.z
    last_in_row = 1 if emitted_real_count == (row_weight - 1) else 0
    is_extra_userdata = 1 if 64 <= entry.col < rate.payload_cols else 0

    mask_flag = 0
    delta_to_last = 0
    if rate.extra_bits_of_parity > 0:
        if rate.fade[row][entry.col] == 1:
            mask_flag = 1
        elif rate.occupied[row][entry.col] == 1 and row == rate.m - 1:
            mask_flag = 1

        last_row_shift = rate.shifts[rate.m - 1][entry.col]
        if last_row_shift >= 0:
            delta_to_last = (last_row_shift - entry.shift + rate.z) % rate.z

    word = 0
    word |= entry.col & 0x7F
    word |= (entry.shift & 0x1FF) << 7
    word |= (pre_entry.row & 0x1F) << 16
    word |= (shift_delta & 0x1FF) << 21
    word |= (last_in_row & 0x1) << 30
    word |= (is_extra_userdata & 0x1) << 31
    word |= (mask_flag & 0x1) << 32
    word |= (delta_to_last & 0x1FF) << 33
    return word & ((1 << 42) - 1)


def write_rdec_sched(rate: RateMatrix, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / f"rdec_sched_{rate.m}x{rate.n}ex{rate.z}_w4.txt"
    row_entries, col_entries = build_qc_entries(rate)
    extra_cols = list(range(64, rate.payload_cols))

    lines: list[str] = []
    sched_index = 0
    for row in range(rate.m):
        row_weight = len(row_entries[row])
        extra_entry_by_col = {
            col: next((entry for entry in row_entries[row] if entry.col == col), None)
            for col in extra_cols
        }

        non_extra_entries: list[Entry] = []
        for j in range(1, rate.m):
            target_pre_row = (row + j) % rate.m
            for entry in row_entries[row]:
                pre_entry = prev_in_col(entry, col_entries)
                if pre_entry.row != target_pre_row:
                    continue
                if 64 <= entry.col < rate.payload_cols:
                    continue
                non_extra_entries.append(entry)

        if extra_cols:
            if len(non_extra_entries) >= 3:
                extra_insert_pos = len(non_extra_entries) - 2
            elif len(non_extra_entries) >= 1:
                extra_insert_pos = 1
            else:
                extra_insert_pos = 0
        else:
            extra_insert_pos = 0

        emitted_real_count = 0
        extra_block_emitted = False

        def emit_real(entry: Entry) -> None:
            nonlocal sched_index
            nonlocal emitted_real_count
            pre_entry = prev_in_col(entry, col_entries)
            emitted_real_count += 1
            word = pack_sched_word(
                rate=rate,
                row=row,
                entry=entry,
                pre_entry=pre_entry,
                emitted_real_count=emitted_real_count,
                row_weight=row_weight,
            )
            lines.append(f"9'd{sched_index:<6}:mmem_rdt=42'h{word:011X};")
            sched_index += 1

        def emit_extra_block() -> None:
            nonlocal sched_index
            for extra_col in extra_cols:
                entry = extra_entry_by_col[extra_col]
                if entry is not None:
                    emit_real(entry)
                else:
                    lines.append(
                        f"9'd{sched_index:<6}:mmem_rdt=42'h{0xFFFFFFFFFFF:011X};"
                    )
                    sched_index += 1

        for index in range(len(non_extra_entries) + 1):
            if extra_cols and not extra_block_emitted and index == extra_insert_pos:
                emit_extra_block()
                extra_block_emitted = True
            if index < len(non_extra_entries):
                emit_real(non_extra_entries[index])

        if extra_cols and not extra_block_emitted:
            emit_extra_block()

    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def gcd_extended(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return a, 1, 0
    gcd, x1, y1 = gcd_extended(b, a % b)
    return gcd, y1, x1 - (a // b) * y1


def mod_inverse(value: int, modulus: int) -> int:
    gcd, x, _ = gcd_extended(value, modulus)
    if gcd != 1:
        raise ValueError(f"No modular inverse for {value} mod {modulus}")
    return x % modulus


def row_values_for_payload(rate: RateMatrix, row: int, payload_cols: int) -> list[int]:
    if payload_cols > rate.payload_cols:
        raise ValueError(
            f"Requested payload_cols={payload_cols} exceeds source payload_cols={rate.payload_cols}"
        )
    payload_part = rate.shifts[row][:payload_cols]
    parity_part = rate.shifts[row][rate.payload_cols : rate.payload_cols + rate.m]
    return payload_part + parity_part


def active_shifts_for_row(rate: RateMatrix, row: int, payload_cols: int | None = None) -> list[int]:
    actual_payload_cols = rate.payload_cols if payload_cols is None else payload_cols
    return [shift for shift in row_values_for_payload(rate, row, actual_payload_cols) if shift >= 0]


def derive_first_element_from_last(active_shifts: list[int], delta: int, z: int) -> int:
    if not active_shifts:
        raise ValueError("Cannot derive first_element from an empty active-shift list")
    if delta == 0 or len(active_shifts) == 1:
        return active_shifts[0]
    return (active_shifts[-1] - ((len(active_shifts) - 1) * delta)) % z


def build_lut_configs(
    rates: list[RateMatrix], payload_cols_min: int, payload_cols_max: int
) -> list[LutConfig]:
    configs: list[LutConfig] = []
    for rate in sorted(rates, key=lambda item: item.m):
        for payload_cols in range(payload_cols_min, payload_cols_max + 1):
            if payload_cols > rate.payload_cols:
                continue
            configs.append(LutConfig(m=rate.m, payload_cols=payload_cols, rate=rate))
    if not configs:
        raise ValueError("No LUT configs built")
    return configs


def derive_row_delta_from_rates(rates: list[RateMatrix], row: int, z: int) -> int:
    deltas: set[int] = set()
    for rate in rates:
        if row >= rate.m:
            continue
        shifts = active_shifts_for_row(rate, row, rate.payload_cols)
        for left, right in zip(shifts, shifts[1:]):
            delta = (right - left) % z
            if delta != 0:
                deltas.add(delta)
    if not deltas:
        return 0
    ordered = sorted(deltas)
    base_delta = ordered[0]
    for value in ordered[1:]:
        base_delta = math.gcd(base_delta, value)
    if base_delta == 0:
        return 0
    for value in ordered:
        if value % base_delta != 0:
            raise ValueError(f"Inconsistent delta for row slot {row}: {ordered}")
    return base_delta


def derive_rtl_wrap_tuple(first_element: int, last_element: int, delta: int, z: int, k_limit: int = 50) -> tuple[int, int, int, int]:
    init_base = (first_element - delta + z) % z
    for shift_dist in range(z):
        shift_base = (init_base + shift_dist) % z
        for delta_num in range(k_limit):
            if (last_element + shift_base + delta_num * delta + shift_dist) % z == shift_base:
                return init_base, shift_base, shift_dist, delta_num
    raise ValueError(
        f"Unable to derive RTL wrap tuple for first={first_element}, last={last_element}, delta={delta}, z={z}"
    )


def solve_delta_num_with_fixed_shift_base(
    first_element: int,
    last_element: int,
    delta: int,
    shift_base: int,
    z: int,
    k_limit: int = 100,
) -> tuple[int, int] | None:
    init_base = (first_element - delta + z) % z
    shift_dist = (shift_base - init_base + z) % z
    for delta_num in range(k_limit):
        if (last_element + shift_base + delta_num * delta + shift_dist) % z == shift_base:
            return shift_dist, delta_num
    return None


def shift_offsets_search_order(z: int) -> list[int]:
    order = [0]
    for offset in range(1, z):
        order.append(offset)
        order.append(-offset)
    return order


def derive_wrap_luts(
    rates: list[RateMatrix], z: int, payload_cols_min: int, payload_cols_max: int
) -> dict[str, object]:
    configs = build_lut_configs(rates, payload_cols_min, payload_cols_max)
    config_by_m_k = {(config.m, config.payload_cols): config for config in configs}
    max_m = max(rate.m for rate in rates)
    min_m = min(rate.m for rate in rates)
    deltas = [derive_row_delta_from_rates(rates, row, z) for row in range(max_m)]
    wrap_base_tmp: list[int] = []
    wrap_base_delta_rows: list[list[int]] = []
    rtl_init_base_rows: list[int] = []
    rtl_shift_dist_rows: list[list[int]] = []
    rtl_delta_num_rows: list[list[int]] = []
    lut_labels = [config.label for config in configs]
    lut_index = [
        {
            "index": idx,
            "m": config.m,
            "payload_cols": config.payload_cols,
            "extra_rows": config.m - min_m,
            "extra_userdata_cols": config.payload_cols - payload_cols_min,
            "label": config.label,
        }
        for idx, config in enumerate(configs)
    ]

    for row in range(max_m):
        delta = deltas[row]
        row_configs = [config for config in configs if row < config.m]
        if not row_configs:
            raise ValueError(f"No rates contain row slot {row}")

        ref_config = next(
            (config for config in row_configs if config.m == max_m and config.payload_cols == payload_cols_max),
            None,
        )
        if ref_config is None:
            raise ValueError(f"No reference config M={max_m},K={payload_cols_max} for row slot {row}")

        ref_active = active_shifts_for_row(ref_config.rate, row, ref_config.payload_cols)
        if not ref_active:
            raise ValueError(
                f"Reference config M={ref_config.m},K={ref_config.payload_cols} row slot {row} has no active shifts"
            )

        rtl_init_base, rtl_shift_base, _, _ = derive_rtl_wrap_tuple(
            ref_active[0], ref_active[-1], delta, z, k_limit=100
        )
        chosen_shift_base: int | None = None
        chosen_init_base: int | None = None
        chosen_shift_dists: list[int] | None = None
        chosen_delta_nums: list[int] | None = None

        for offset in shift_offsets_search_order(z):
            candidate_shift_base = (rtl_shift_base + offset) % z
            candidate_shift_dists: list[int] = []
            candidate_delta_nums: list[int] = []
            ok = True
            for config in configs:
                if row >= config.m:
                    candidate_shift_dists.append(0)
                    candidate_delta_nums.append(0)
                    continue

                row_shifts = active_shifts_for_row(config.rate, row, config.payload_cols)
                if not row_shifts:
                    candidate_shift_dists.append(0)
                    candidate_delta_nums.append(0)
                    continue

                solved = solve_delta_num_with_fixed_shift_base(
                    row_shifts[0],
                    row_shifts[-1],
                    delta,
                    candidate_shift_base,
                    z,
                    k_limit=100,
                )
                if solved is None:
                    ok = False
                    break
                shift_dist, delta_num = solved
                if delta_num >= 64:
                    ok = False
                    break
                candidate_shift_dists.append(shift_dist)
                candidate_delta_nums.append(delta_num)

            if ok:
                chosen_shift_base = candidate_shift_base
                chosen_init_base = (ref_active[0] - delta + z) % z
                chosen_shift_dists = candidate_shift_dists
                chosen_delta_nums = candidate_delta_nums
                break

        if chosen_shift_base is None or chosen_shift_dists is None or chosen_delta_nums is None or chosen_init_base is None:
            raise ValueError(f"Unable to find reusable shift_base for row slot {row}")

        rtl_init_base_rows.append(chosen_init_base)
        wrap_base_tmp.append(chosen_shift_base)

        wrap_base_delta_rows.append(chosen_delta_nums)
        rtl_shift_dist_rows.append(chosen_shift_dists)
        rtl_delta_num_rows.append(chosen_delta_nums)

    first_mask_shift_by_m: dict[int, list[int]] = {}
    first_mask_shift_debug: dict[int, list[dict[str, int]]] = {}
    for m in range(min_m, max_m + 1):
        row = m - 1
        delta = deltas[row]
        values: list[int] = []
        debug_items: list[dict[str, int]] = []
        for payload_cols in range(payload_cols_min, payload_cols_max + 1):
            config = config_by_m_k.get((m, payload_cols))
            if config is None:
                raise ValueError(f"Missing LUT config for M={m},K={payload_cols}")
            row_values = row_values_for_payload(config.rate, row, payload_cols)
            active_shifts = [shift for shift in row_values if shift >= 0]
            if not active_shifts:
                raise ValueError(f"No active CPM for first_mask_shift at M={m},K={payload_cols},row={row}")
            # FE is reconstructed from the rightmost active CPM by walking left with DELTA.
            first_element = derive_first_element_from_last(active_shifts, delta, z)
            parity_values = row_values[payload_cols:]
            parity_first = next((shift for shift in parity_values if shift >= 0), None)
            if parity_first is None:
                raise ValueError(
                    f"No active parity CPM for first_mask_shift at M={m},K={payload_cols},row={row}"
                )
            first_mask_shift = (z + first_element - parity_first) % z
            values.append(first_mask_shift)
            debug_items.append(
                {
                    "m": m,
                    "payload_cols": payload_cols,
                    "row": row,
                    "delta": delta,
                    "first_element": first_element,
                    "parity_first": parity_first,
                    "first_mask_shift": first_mask_shift,
                }
            )
        first_mask_shift_by_m[m] = values
        first_mask_shift_debug[m] = debug_items

    return {
        "min_m": min_m,
        "max_m": max_m,
        "m_values": [rate.m for rate in rates],
        "lut_payload_cols": list(range(payload_cols_min, payload_cols_max + 1)),
        "lut_labels": lut_labels,
        "lut_index": lut_index,
        "delta": deltas,
        "wrap_base_tmp": wrap_base_tmp,
        "rtl_init_base": rtl_init_base_rows,
        "wrap_base_delta": wrap_base_delta_rows,
        "rtl_shift_dist": rtl_shift_dist_rows,
        "rtl_delta_num": rtl_delta_num_rows,
        "first_mask_shift_by_m": first_mask_shift_by_m,
        "first_mask_shift_debug": first_mask_shift_debug,
    }


def format_concat(values: Iterable[int], width: int, reverse: bool = True) -> str:
    ordered = list(values)
    if reverse:
        ordered = list(reversed(ordered))
    return "{" + ", ".join(f"{width}'d{value}" for value in ordered) + "}"


def pack_column_bitmap(grid: list[list[int]], rows: int, col: int) -> int:
    value = 0
    for row in range(rows):
        if grid[row][col]:
            value |= 1 << row
    return value


def write_matrix_assignment_svh(rates: list[RateMatrix], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / "ldpc_matrix_assignments_maxk67.svh"
    max_rows = max(rate.m for rate in rates)
    payload_cols = max(rate.payload_cols for rate in rates)
    parity_slots = max_rows
    lines = [
        "// Auto-generated by scripts/ibex_matrix_family/generate_family_artifacts.py",
        "// Full-width (K=67) occupied/fade matrix assignments.",
        "// Bit packing uses row0 as LSB; all vectors are emitted as 17-bit values.",
        "",
        "// occupied",
    ]

    for rate in sorted(rates, key=lambda item: item.m):
        for col in range(payload_cols):
            bitmap = pack_column_bitmap(rate.occupied, rate.m, col)
            lines.append(
                f"assign ldpc_matrix_payload_occupied[{rate.m:2d}][{col:2d}] = 17'b{bitmap:017b};"
            )
        for parity_idx in range(parity_slots):
            if parity_idx < rate.m:
                col = rate.payload_cols + parity_idx
                bitmap = pack_column_bitmap(rate.occupied, rate.m, col)
            else:
                bitmap = 0
            lines.append(
                f"assign ldpc_matrix_parity_occupied[{rate.m:2d}][{parity_idx:2d}] = 17'b{bitmap:017b};"
            )

    lines.append("")
    lines.append("// fade")

    for rate in sorted(rates, key=lambda item: item.m):
        for col in range(payload_cols):
            bitmap = pack_column_bitmap(rate.fade, rate.m, col)
            lines.append(
                f"assign ldpc_matrix_payload_fade[{rate.m:2d}][{col:2d}] = 17'b{bitmap:017b};"
            )
        for parity_idx in range(parity_slots):
            if parity_idx < rate.m:
                col = rate.payload_cols + parity_idx
                bitmap = pack_column_bitmap(rate.fade, rate.m, col)
            else:
                bitmap = 0
            lines.append(
                f"assign ldpc_matrix_parity_fade[{rate.m:2d}][{parity_idx:2d}] = 17'b{bitmap:017b};"
            )

    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def write_lut_outputs(lut: dict[str, object], out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "matrix_lut_summary.json"
    svh_path = out_dir / "ldpc_matrix_lut.svh"
    csv_path = out_dir / "wrap_base_delta_table.csv"
    adjusted_json_path = out_dir / "matrix_lut_adjusted_summary.json"
    adjusted_svh_path = out_dir / "ldpc_matrix_lut_adjusted_extra_rows.svh"
    adjusted_csv_path = out_dir / "wrap_base_delta_adjusted_table.csv"

    json_path.write_text(json.dumps(lut, indent=2) + "\n")

    m_values = lut["m_values"]
    lut_labels = lut["lut_labels"]
    lut_payload_cols = lut["lut_payload_cols"]
    delta_values = lut["delta"]
    wrap_base_tmp = lut["wrap_base_tmp"]
    wrap_base_delta = lut["wrap_base_delta"]
    first_mask_shift_by_m = {
        int(m): values for m, values in lut.get("first_mask_shift_by_m", {}).items()
    }
    max_m = lut["max_m"]
    min_m = lut["min_m"]
    max_extra_rows = max_m - min_m
    k_count = len(lut_payload_cols)

    svh_lines = [
        "// Auto-generated by scripts/ibex_matrix_family/generate_family_artifacts.py",
        f"// Rate order for per-slot delta tables: M={m_values[0]}..{m_values[-1]}, K={lut_payload_cols[0]}..{lut_payload_cols[-1]}",
        f"// Index order is idx = (M-{min_m})*{len(lut_payload_cols)} + (K-{lut_payload_cols[0]}).",
        "// `ldpc_matrix_delta_tmp_adj_ord` and `ldpc_matrix_wrap_base_tmp_adj_ord` are emitted in natural row-slot order.",
        "// `ldpc_matrix_wrap_base_tmp_adj_ord` uses RTL shift_base semantics.",
        "// `ldpc_matrix_first_mask_shiftMM` uses: (Z + FE_last_row - first_parity_CPM_last_row) mod Z.",
        "// FE_last_row is reconstructed from rightmost active CPM by stepping left with DELTA.",
        "// `ldpc_matrix_wrap_base_deltaX_adj_ord` uses RTL delta_num semantics and keeps reversed concatenation for `[W*idx +: W]` access: highest index on the left, index 0 on the right.",
        f"localparam int LDPC_MATRIX_MIN_M = {min_m};",
        f"localparam int LDPC_MATRIX_MAX_M = {max_m};",
        f"localparam int LDPC_MATRIX_LUT_K_MIN = {lut_payload_cols[0]};",
        f"localparam int LDPC_MATRIX_LUT_K_MAX = {lut_payload_cols[-1]};",
        f"localparam int LDPC_MATRIX_RATE_COUNT = {len(m_values)};",
        f"localparam int LDPC_MATRIX_LUT_COUNT = {len(lut_labels)};",
        f"localparam logic [9*{max_m}-1:0] ldpc_matrix_delta_tmp_adj_ord = {format_concat(delta_values, 9, reverse=False)};",
        f"localparam logic [9*{max_m}-1:0] ldpc_matrix_wrap_base_tmp_adj_ord = {format_concat(wrap_base_tmp, 9, reverse=False)};",
    ]

    svh_lines.append("")
    for m in range(min_m, max_m + 1):
        values = first_mask_shift_by_m.get(m)
        if values is None:
            continue
        svh_lines.append(
            f"wire [{len(values)}*9-1:0] ldpc_matrix_first_mask_shift{m:02d} = {format_concat(values, 9, reverse=False)};"
        )

    svh_lines.append("")
    for row, values in enumerate(wrap_base_delta):
        svh_lines.append(
            f"wire [6*{len(values)}-1:0] ldpc_matrix_wrap_base_delta{row}_adj_ord = {format_concat(values, 6)};"
        )

    svh_path.write_text("\n".join(svh_lines) + "\n")

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row_slot", *lut_labels])
        for row, values in enumerate(wrap_base_delta):
            writer.writerow([row, *values])

    adjusted_wrap_base_delta: list[list[int]] = []
    adjusted_labels_by_row: list[list[str]] = []
    adjusted_index_by_row: list[list[dict[str, int | str]]] = []

    for row, values in enumerate(wrap_base_delta):
        offset = max(0, row - 4)
        row_labels: list[str] = []
        row_index: list[dict[str, int | str]] = []
        row_values: list[int] = []
        for adjusted_extra_rows in range(max_extra_rows - offset + 1):
            actual_m = min_m + offset + adjusted_extra_rows
            for payload_cols in lut_payload_cols:
                absolute_idx = (actual_m - min_m) * k_count + (
                    payload_cols - lut_payload_cols[0]
                )
                row_values.append(values[absolute_idx])
                label = f"M{actual_m}_K{payload_cols}"
                row_labels.append(label)
                row_index.append(
                    {
                        "adjusted_index": len(row_values) - 1,
                        "row_slot": row,
                        "offset": offset,
                        "adjusted_extra_rows": adjusted_extra_rows,
                        "m": actual_m,
                        "payload_cols": payload_cols,
                        "label": label,
                    }
                )
        adjusted_wrap_base_delta.append(row_values)
        adjusted_labels_by_row.append(row_labels)
        adjusted_index_by_row.append(row_index)

    adjusted_summary = {
        "min_m": min_m,
        "max_m": max_m,
        "lut_payload_cols": lut_payload_cols,
        "delta": delta_values,
        "wrap_base_tmp": wrap_base_tmp,
        "adjusted_offset_by_row": [max(0, row - 4) for row in range(max_m)],
        "adjusted_labels_by_row": adjusted_labels_by_row,
        "adjusted_index_by_row": adjusted_index_by_row,
        "wrap_base_delta_adjusted": adjusted_wrap_base_delta,
    }
    adjusted_json_path.write_text(json.dumps(adjusted_summary, indent=2) + "\n")

    adjusted_svh_lines = [
        "// Auto-generated by scripts/ibex_matrix_family/generate_family_artifacts.py",
        f"// Adjusted-extra-rows LUTs for legacy DV selectors.",
        f"// For row-slot r, offset = max(0, r-4).",
        f"// Selector for row-slot r is idx = adjusted_extra_rows*{k_count} + (K-{lut_payload_cols[0]}),",
        f"// where adjusted_extra_rows = max((M-{min_m}) - offset, 0).",
        "// `ldpc_matrix_delta_tmp_adj_ord` and `ldpc_matrix_wrap_base_tmp_adj_ord` are emitted in natural row-slot order.",
        "// `ldpc_matrix_wrap_base_tmp_adj_ord` uses RTL shift_base semantics.",
        "// `ldpc_matrix_first_mask_shiftMM` uses: (Z + FE_last_row - first_parity_CPM_last_row) mod Z.",
        "// Each `ldpc_matrix_wrap_base_deltaX_adj_ord` below uses RTL delta_num values with adjusted-extra-rows packing and keeps reversed concatenation for `[W*idx +: W]` access.",
        f"localparam int LDPC_MATRIX_MIN_M = {min_m};",
        f"localparam int LDPC_MATRIX_MAX_M = {max_m};",
        f"localparam int LDPC_MATRIX_LUT_K_MIN = {lut_payload_cols[0]};",
        f"localparam int LDPC_MATRIX_LUT_K_MAX = {lut_payload_cols[-1]};",
        f"localparam logic [9*{max_m}-1:0] ldpc_matrix_delta_tmp_adj_ord = {format_concat(delta_values, 9, reverse=False)};",
        f"localparam logic [9*{max_m}-1:0] ldpc_matrix_wrap_base_tmp_adj_ord = {format_concat(wrap_base_tmp, 9, reverse=False)};",
    ]

    adjusted_svh_lines.append("")
    for m in range(min_m, max_m + 1):
        values = first_mask_shift_by_m.get(m)
        if values is None:
            continue
        adjusted_svh_lines.append(
            f"wire [{len(values)}*9-1:0] ldpc_matrix_first_mask_shift{m:02d} = {format_concat(values, 9, reverse=False)};"
        )

    adjusted_svh_lines.append("")
    for row, values in enumerate(adjusted_wrap_base_delta):
        adjusted_svh_lines.append(
            f"wire [6*{len(values)}-1:0] ldpc_matrix_wrap_base_delta{row}_adj_ord = {format_concat(values, 6)};"
        )

    adjusted_svh_path.write_text("\n".join(adjusted_svh_lines) + "\n")

    with adjusted_csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row_slot", "offset", "adjusted_index", "adjusted_extra_rows", "m", "payload_cols", "label", "value"])
        for row, (index_items, values) in enumerate(zip(adjusted_index_by_row, adjusted_wrap_base_delta)):
            offset = max(0, row - 4)
            for item, value in zip(index_items, values):
                writer.writerow(
                    [
                        row,
                        offset,
                        item["adjusted_index"],
                        item["adjusted_extra_rows"],
                        item["m"],
                        item["payload_cols"],
                        item["label"],
                        value,
                    ]
                )

    return [
        json_path,
        svh_path,
        csv_path,
        adjusted_json_path,
        adjusted_svh_path,
        adjusted_csv_path,
    ]


def write_manifest(
    output_root: Path,
    rates: list[RateMatrix],
    normalized_paths: dict[str, list[Path]],
    bm_schematic_paths: list[Path],
    ens_paths: list[Path],
    rdec_paths: list[Path],
    lut_paths: list[Path],
) -> Path:
    manifest = {
        "output_root": str(output_root),
        "rates": [
            {
                "m": rate.m,
                "n": rate.n,
                "payload_cols": rate.payload_cols,
                "parity_bytes": rate.parity_bytes,
                "extra_bits_of_parity": rate.extra_bits_of_parity,
            }
            for rate in rates
        ],
        "normalized": {
            key: [str(path) for path in paths] for key, paths in normalized_paths.items()
        },
        "bm_schematic": [str(path) for path in bm_schematic_paths],
        "ens": [str(path) for path in ens_paths],
        "rdec_sched": [str(path) for path in rdec_paths],
        "lut": [str(path) for path in lut_paths],
    }
    manifest_path = output_root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate normalized IBEX 13-rate matrix artifacts, encoder schedules, RDEC schedules, and LUT outputs."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("IBEX/ibex_matrix"),
        help="Root directory containing per-rate matrix folders like 13x80/",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("IBEX/ibex_matrix_flat_13rate"),
        help="Flat output directory for normalized matrices and generated artifacts.",
    )
    parser.add_argument(
        "--payload-cols",
        type=int,
        default=67,
        help="Payload/base userdata column count K.",
    )
    parser.add_argument(
        "--circulant-size",
        type=int,
        default=512,
        help="Circulant size Z.",
    )
    parser.add_argument(
        "--m5-parity-bytes",
        type=int,
        default=320,
        help="Parity bytes override for M=5. Other rates use M*(Z/8)-1.",
    )
    parser.add_argument(
        "--lut-k-min",
        type=int,
        default=64,
        help="Minimum payload column count K used when expanding LUT tables.",
    )
    parser.add_argument(
        "--lut-k-max",
        type=int,
        default=67,
        help="Maximum payload column count K used when expanding LUT tables.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    rates = discover_rates(
        input_root=args.input_root,
        payload_cols=args.payload_cols,
        z=args.circulant_size,
        m5_parity_bytes=args.m5_parity_bytes,
    )

    normalized_paths: dict[str, list[Path]] = {
        "matrix": [],
        "fade_matrix": [],
        "occupied_matrix": [],
    }
    ens_paths: list[Path] = []
    rdec_paths: list[Path] = []

    matrix_dir = args.output_root / "matrix"
    fade_dir = args.output_root / "fade_matrix"
    occupied_dir = args.output_root / "occupied_matrix"
    bm_schematic_dir = args.output_root / "bm_schematic"
    ens_dir = args.output_root / "ens"
    rdec_dir = args.output_root / "rdec_sched"
    lut_dir = args.output_root / "lut"
    bm_schematic_paths: list[Path] = []

    for rate in rates:
        normalized_paths["matrix"].append(
            copy_normalized(
                rate.matrix_path,
                matrix_dir,
                normalize_output_name(rate.prefix, "matrix"),
            )
        )
        normalized_paths["fade_matrix"].append(
            copy_normalized(
                rate.fade_path,
                fade_dir,
                normalize_output_name(rate.prefix, "fade"),
            )
        )
        normalized_paths["occupied_matrix"].append(
            copy_normalized(
                rate.occupied_path,
                occupied_dir,
                normalize_output_name(rate.prefix, "occupied"),
            )
        )
        bm_schematic_paths.append(write_bm_schematic(rate, bm_schematic_dir))
        ens_paths.extend(write_ens_files(rate, ens_dir))
        rdec_paths.append(write_rdec_sched(rate, rdec_dir))

    lut = derive_wrap_luts(
        rates,
        z=args.circulant_size,
        payload_cols_min=args.lut_k_min,
        payload_cols_max=args.lut_k_max,
    )
    lut_paths = write_lut_outputs(lut, lut_dir)
    lut_paths.append(write_matrix_assignment_svh(rates, lut_dir))
    manifest_path = write_manifest(
        output_root=args.output_root,
        rates=rates,
        normalized_paths=normalized_paths,
        bm_schematic_paths=bm_schematic_paths,
        ens_paths=ens_paths,
        rdec_paths=rdec_paths,
        lut_paths=lut_paths,
    )

    print(f"[ibex-matrix-family] Generated {len(rates)} rates under {args.output_root}")
    print(f"[ibex-matrix-family] Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
