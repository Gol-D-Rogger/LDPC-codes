#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


FULL_PAYLOAD_COLS = 67
DELTA_TABLE = [
    0,
    13,
    19,
    29,
    41,
    67,
    73,
    79,
    91,
    97,
    103,
    111,
    119,
    127,
    131,
    137,
    149,
]
MATRIX_RE = re.compile(r"^LDPC_(?P<m>\d+)x(?P<n>\d+)ex(?P<z>\d+)_w4_dense5_QC_H_1\.txt$")
PLACEHOLDER_WORD = 0xFFFFFFFFFFF


@dataclass(frozen=True)
class Entry:
    row: int
    col: int
    shift: int


@dataclass
class RateMatrix:
    m: int
    n: int
    z: int
    parity_bytes: int
    shifts: list[list[int]]
    occupied: list[list[int]]
    fade: list[list[int]]

    @property
    def payload_cols(self) -> int:
        return self.n - self.m

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


def read_int_matrix(path: Path) -> list[list[int]]:
    rows: list[list[int]] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        rows.append([int(token) for token in line.split()])
    return rows


def parity_bytes_for_m(m: int, z: int, m5_parity_bytes: int) -> int:
    if m == 5:
        return m5_parity_bytes
    return m * (z // 8) - 1


def load_full_rates(matrix_root: Path, z: int, m5_parity_bytes: int) -> list[RateMatrix]:
    matrix_dir = matrix_root / "matrix"
    occupied_dir = matrix_root / "occupied_matrix"
    fade_dir = matrix_root / "fade_matrix"

    rates: list[RateMatrix] = []
    for matrix_path in sorted(matrix_dir.glob("LDPC_*_QC_H_1.txt")):
        match = MATRIX_RE.match(matrix_path.name)
        if not match:
            continue

        m = int(match.group("m"))
        n = int(match.group("n"))
        file_z = int(match.group("z"))
        if file_z != z:
            continue
        if n - m != FULL_PAYLOAD_COLS:
            continue

        occupied_path = (
            occupied_dir / f"LDPC_{m}x{n}ex{z}_w4_dense5_occupied_1.txt"
        )
        fade_path = fade_dir / f"LDPC_{m}x{n}ex{z}_w4_dense5_fade_1.txt"
        if not occupied_path.exists() or not fade_path.exists():
            raise FileNotFoundError(f"Missing occupied/fade file for {matrix_path.name}")

        shifts = read_int_matrix(matrix_path)
        occupied = read_int_matrix(occupied_path)
        fade = read_int_matrix(fade_path)

        if len(shifts) != m or any(len(row) != n for row in shifts):
            raise ValueError(f"Invalid shift matrix shape in {matrix_path}")
        if len(occupied) != m or any(len(row) != n for row in occupied):
            raise ValueError(f"Invalid occupied matrix shape in {occupied_path}")
        if len(fade) != m or any(len(row) != n for row in fade):
            raise ValueError(f"Invalid fade matrix shape in {fade_path}")

        rates.append(
            RateMatrix(
                m=m,
                n=n,
                z=z,
                parity_bytes=parity_bytes_for_m(m, z, m5_parity_bytes),
                shifts=shifts,
                occupied=occupied,
                fade=fade,
            )
        )

    if not rates:
        raise FileNotFoundError(f"No full-width K=67 matrices found under {matrix_root}")
    return rates


def crop_and_regenerate(full_rate: RateMatrix, target_k: int) -> RateMatrix:
    if target_k < 64 or target_k > FULL_PAYLOAD_COLS:
        raise ValueError(f"Unsupported target K={target_k}")

    target_n = full_rate.m + target_k
    source_cols = full_rate.n
    m = full_rate.m
    z = full_rate.z

    shifts = [[-1 for _ in range(target_n)] for _ in range(m)]
    occupied = [[0 for _ in range(target_n)] for _ in range(m)]
    fade = [[0 for _ in range(target_n)] for _ in range(m)]
    last_element = [0 for _ in range(m)]
    index = [0 for _ in range(m)]

    for row in range(m):
        for col in range(target_n):
            src_col = col if col < target_k else (source_cols - target_n + col)
            src_shift = full_rate.shifts[row][src_col]
            src_occupied = full_rate.occupied[row][src_col]
            src_fade = full_rate.fade[row][src_col]
            if src_fade != 0:
                fade[row][col] = 1
            elif src_occupied != 0 and src_shift >= 0:
                occupied[row][col] = 1
            if src_shift >= 0 and (occupied[row][col] or fade[row][col]):
                last_element[row] = src_shift

    for row in range(m):
        index[row] = last_element[row]

    for rev in range(target_n):
        col = target_n - 1 - rev
        for row in range(m):
            if occupied[row][col] or fade[row][col]:
                shifts[row][col] = index[row]
                index[row] = (index[row] + z - DELTA_TABLE[row]) % z

    return RateMatrix(
        m=m,
        n=target_n,
        z=z,
        parity_bytes=parity_bytes_for_m(m, z, m5_parity_bytes=320),
        shifts=shifts,
        occupied=occupied,
        fade=fade,
    )


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
    for index, candidate in enumerate(entries):
        if candidate.row == entry.row:
            return entries[index - 1] if index > 0 else entries[-1]
    raise ValueError(f"Entry not found in column list: {entry}")


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


def build_rdec_sched_lines(rate: RateMatrix) -> list[str]:
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
        for distance in range(1, rate.m):
            target_pre_row = (row + distance) % rate.m
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
                        f"9'd{sched_index:<6}:mmem_rdt=42'h{PLACEHOLDER_WORD:011X};"
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

    return lines


def validate_sched(lines: list[str], max_col: int) -> None:
    for line in lines:
        match = re.search(r"42'h([0-9A-Fa-f]+)", line)
        if not match:
            continue
        word = int(match.group(1), 16)
        if word == PLACEHOLDER_WORD:
            continue
        col = word & 0x7F
        if col > max_col:
            raise ValueError(
                f"Found scheduler entry with col={col} > max_col={max_col}: {line}"
            )


def write_sched(rate: RateMatrix, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = build_rdec_sched_lines(rate)
    validate_sched(lines, rate.n - 1)
    output_path = output_dir / f"rdec_sched_{rate.m}x{rate.n}ex{rate.z}_w4.txt"
    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate cropped RDEC_sched files from full-width K=67 IBEX matrices. "
            "The cropped view preserves delta spacing by regenerating shifts after crop."
        )
    )
    parser.add_argument(
        "--matrix-root",
        type=Path,
        default=Path("IBEX/ibex_matrix_flat_13rate"),
        help="Root containing matrix/occupied_matrix/fade_matrix directories.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("IBEX/ibex_matrix_flat_13rate/rdec_sched_cropped_k64_67"),
        help="Directory to write generated RDEC_sched files.",
    )
    parser.add_argument(
        "--k-values",
        type=int,
        nargs="+",
        default=[67, 66, 65, 64],
        help="Target K values to generate.",
    )
    parser.add_argument(
        "--m-values",
        type=int,
        nargs="*",
        default=None,
        help="Optional subset of M values to generate.",
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
        help="Parity bytes used for M=5.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    full_rates = load_full_rates(
        matrix_root=args.matrix_root,
        z=args.circulant_size,
        m5_parity_bytes=args.m5_parity_bytes,
    )

    if args.m_values:
        wanted = set(args.m_values)
        full_rates = [rate for rate in full_rates if rate.m in wanted]
        if not full_rates:
            raise ValueError(f"No full-width rates found for M values {sorted(wanted)}")

    generated: list[Path] = []
    for full_rate in full_rates:
        for target_k in args.k_values:
            cropped = crop_and_regenerate(full_rate, target_k)
            generated.append(write_sched(cropped, args.output_root))

    print(
        f"[cropped-rdec] Generated {len(generated)} files under {args.output_root}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
