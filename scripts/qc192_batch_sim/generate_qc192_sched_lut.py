#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


SCHED_LINE_RE = re.compile(
    r"^\s*(\d+)'d\s*([0-9]+)\s*:\s*mmem_rdt\s*=\s*(\d+)'h([0-9A-Fa-f]+)\s*;\s*$"
)


@dataclass(frozen=True)
class SizeContext:
    m: int
    n: int
    size_name: str
    source_dir: Path
    run_matrix_root: Path
    config_path: Path
    log_path: Path
    raw_dir: Path


@dataclass(frozen=True)
class SchedEntry:
    addr_bits: int
    address: int
    data_bits: int
    data_hex: str


def parse_m_range(expr: str) -> list[int]:
    values: set[int] = set()
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
            values.update(range(lo, hi + 1))
        else:
            values.add(int(part))
    if not values:
        raise ValueError(f"Invalid --m-range: {expr}")
    return sorted(values)


def parse_size_name(size_name: str) -> tuple[int, int] | None:
    match = re.match(r"^([0-9]+)x([0-9]+)$", size_name)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def has_flat_matrix_layout(matrix_root: Path) -> bool:
    return all((matrix_root / sub).is_dir() for sub in ("matrix", "fade_matrix", "occupied_matrix"))


def discover_size_dir(matrix_root: Path, m: int) -> Path:
    candidates: list[Path] = []
    for child in matrix_root.iterdir():
        if not child.is_dir():
            continue
        parsed = parse_size_name(child.name)
        if parsed is None:
            continue
        if parsed[0] == m:
            candidates.append(child)
    if not candidates:
        raise FileNotFoundError(f"No size directory for M={m} under {matrix_root}")
    candidates.sort(key=lambda p: p.name)
    return candidates[0]


def find_flat_matrix_triplet(matrix_root: Path, m: int, qc: int) -> tuple[int, Path, Path, Path]:
    n = m + 178
    matrix_path = matrix_root / "matrix" / f"LDPC_{m}x{n}ex{qc}_w4_dense5_QC_H.txt"
    fade_path = matrix_root / "fade_matrix" / f"LDPC_{m}x{n}ex{qc}_w4_dense5_QC_H_fade.txt"
    occupied_path = matrix_root / "occupied_matrix" / f"LDPC_{m}x{n}ex{qc}_w4_dense5_QC_H_occupied.txt"
    if not matrix_path.exists():
        raise FileNotFoundError(f"Missing flat matrix file: {matrix_path}")
    if not fade_path.exists():
        raise FileNotFoundError(f"Missing flat fade file: {fade_path}")
    if not occupied_path.exists():
        raise FileNotFoundError(f"Missing flat occupied file: {occupied_path}")
    return n, matrix_path, fade_path, occupied_path


def patch_leading_int(line: str, value: int) -> str:
    if "//" in line:
        left, right = line.split("//", 1)
        comment = f"//{right.rstrip()}"
    else:
        left, comment = line, ""
    patched, count = re.subn(r"^(\s*)[0-9]+", rf"\g<1>{int(value)}", left, count=1)
    if count != 1:
        raise ValueError(f"Cannot patch numeric line: {line!r}")
    patched = patched.rstrip()
    if comment:
        return f"{patched} {comment}"
    return patched


def build_config_text(template_text: str, *, m: int) -> str:
    lines = template_text.splitlines()
    if len(lines) < 5:
        raise ValueError("Template cnfg must contain at least 5 lines")

    # Per user requirement, only patch line 5:
    # parity size (B) = m * 24 - 4
    lines[4] = patch_leading_int(lines[4], m * 24 - 4)
    return "\n".join(lines) + "\n"


def stage_flat_matrix_layout(
    *,
    matrix_root: Path,
    out_dir: Path,
    m: int,
    qc: int,
    matrix_id: int,
) -> tuple[int, Path]:
    n, matrix_src, fade_src, occupied_src = find_flat_matrix_triplet(matrix_root, m, qc)

    staged_root = out_dir / f"{m}x{n}" / "staged_matrix_root"
    flat_matrix_dir = staged_root / "matrix"
    flat_fade_dir = staged_root / "fade_matrix"
    flat_occupied_dir = staged_root / "occupied_matrix"
    size_matrix_dir = staged_root / f"{m}x{n}" / "matrix"
    size_fade_dir = staged_root / f"{m}x{n}" / "fade_matrix"
    size_occupied_dir = staged_root / f"{m}x{n}" / "occupied_matrix"

    for path in (
        flat_matrix_dir,
        flat_fade_dir,
        flat_occupied_dir,
        size_matrix_dir,
        size_fade_dir,
        size_occupied_dir,
    ):
        path.mkdir(parents=True, exist_ok=True)

    base_matrix = f"LDPC_{m}x{n}ex{qc}_w4_dense5_QC_H.txt"
    base_fade = f"LDPC_{m}x{n}ex{qc}_w4_dense5_QC_H_fade.txt"
    base_occupied = f"LDPC_{m}x{n}ex{qc}_w4_dense5_QC_H_occupied.txt"

    compat_matrix = f"LDPC_{m}x{n}ex{qc}_w4_dense5_QC_H_1_{matrix_id}.txt"
    compat_fade = f"LDPC_{m}x{n}ex{qc}_w4_dense5_fade_1_{matrix_id}.txt"
    compat_occupied = f"LDPC_{m}x{n}ex{qc}_w4_dense5_occupied_1_{matrix_id}.txt"

    for dst_dir in (flat_matrix_dir, size_matrix_dir):
        shutil.copyfile(matrix_src, dst_dir / base_matrix)
    for dst_dir in (flat_fade_dir, size_fade_dir):
        shutil.copyfile(fade_src, dst_dir / base_fade)
    for dst_dir in (flat_occupied_dir, size_occupied_dir):
        shutil.copyfile(occupied_src, dst_dir / base_occupied)

    # Also stage src2-style compatibility names under the sized layout.
    shutil.copyfile(matrix_src, size_matrix_dir / compat_matrix)
    shutil.copyfile(fade_src, size_fade_dir / compat_fade)
    shutil.copyfile(occupied_src, size_occupied_dir / compat_occupied)

    return n, staged_root


def snapshot_matching_files(output_dir: Path, pattern: str) -> dict[Path, int]:
    snap: dict[Path, int] = {}
    for path in output_dir.glob(pattern):
        try:
            snap[path] = path.stat().st_mtime_ns
        except FileNotFoundError:
            continue
    return snap


def pick_updated_file(output_dir: Path, pattern: str, before: dict[Path, int]) -> Path:
    after = snapshot_matching_files(output_dir, pattern)
    updated = [p for p, mtime in after.items() if before.get(p) != mtime]
    if len(updated) == 1:
        return updated[0]
    if len(updated) > 1:
        updated.sort(key=lambda p: after[p], reverse=True)
        return updated[0]
    if after:
        newest = max(after.items(), key=lambda kv: kv[1])[0]
        return newest
    raise FileNotFoundError(f"No generated file matching {pattern} under {output_dir}")


def parse_sched_entries(raw_lines: list[str], *, source: Path) -> list[SchedEntry]:
    entries: list[SchedEntry] = []
    for line_no, line in enumerate(raw_lines, 1):
        stripped = line.strip()
        if not stripped:
            continue
        match = SCHED_LINE_RE.match(stripped)
        if not match:
            raise ValueError(f"{source}:{line_no}: unsupported sched line: {line}")
        entries.append(
            SchedEntry(
                addr_bits=int(match.group(1)),
                address=int(match.group(2)),
                data_bits=int(match.group(3)),
                data_hex=match.group(4).upper(),
            )
        )
    if not entries:
        raise ValueError(f"{source}: no sched entries found")
    return entries


def normalize_sched_entries(entries: list[SchedEntry]) -> tuple[list[str], int]:
    addr_bits_set = {entry.addr_bits for entry in entries}
    data_bits_set = {entry.data_bits for entry in entries}
    if len(addr_bits_set) != 1:
        raise ValueError(f"Inconsistent address widths in sched lines: {sorted(addr_bits_set)}")
    if len(data_bits_set) != 1:
        raise ValueError(f"Inconsistent data widths in sched lines: {sorted(data_bits_set)}")

    addr_bits = entries[0].addr_bits
    data_bits = entries[0].data_bits
    addr_width = max(len(str(entry.address)) for entry in entries)

    lines = [
        f"{addr_bits}'d{entry.address:<{addr_width}}:    mmem_rdt={data_bits}'h{entry.data_hex};"
        for entry in entries
    ]
    return lines, data_bits


def build_verilog_from_template(
    template_text: str,
    *,
    module_name: str,
    sched_lines: list[str],
    default_data_bits: int,
) -> str:
    text, module_subs = re.subn(
        r"(^\s*module\s+)\w+(\s+\(/\*AUTOARG\*/)",
        rf"\1{module_name}\2",
        template_text,
        count=1,
        flags=re.M,
    )
    if module_subs != 1:
        raise ValueError("Failed to replace module line in template")

    body_lines = [
        "always @(*)",
        "begin",
        "    case(mmem_radr)",
    ]
    body_lines.extend(f"    {line}" for line in sched_lines)
    body_lines.append(f"    default: mmem_rdt={default_data_bits}'h0;")
    body_lines.append("    endcase")
    body_lines.append("end")
    body = "\n".join(body_lines) + "\n"

    text, block_subs = re.subn(
        r"always @\(\*\)\s*begin.*?^\s*end\s*$",
        body,
        text,
        count=1,
        flags=re.S | re.M,
    )
    if block_subs != 1:
        raise ValueError("Failed to replace always/case block in template")
    return text


def prepare_size_context(
    *,
    m: int,
    matrix_root: Path,
    out_dir: Path,
    qc: int,
    matrix_id: int,
) -> SizeContext:
    if has_flat_matrix_layout(matrix_root):
        n, run_matrix_root = stage_flat_matrix_layout(
            matrix_root=matrix_root,
            out_dir=out_dir,
            m=m,
            qc=qc,
            matrix_id=matrix_id,
        )
        source_dir = matrix_root
        size_name = f"{m}x{n}"
    else:
        source_dir = discover_size_dir(matrix_root, m)
        parsed = parse_size_name(source_dir.name)
        if parsed is None:
            raise ValueError(f"Invalid size directory name: {source_dir.name}")
        _, n = parsed
        run_matrix_root = matrix_root
        size_name = source_dir.name

    size_out = out_dir / size_name
    config_dir = size_out / "configs"
    log_dir = size_out / "logs"
    raw_dir = size_out / "raw_sched"
    for path in (config_dir, log_dir, raw_dir):
        path.mkdir(parents=True, exist_ok=True)

    return SizeContext(
        m=m,
        n=n,
        size_name=size_name,
        source_dir=source_dir,
        run_matrix_root=run_matrix_root,
        config_path=config_dir / f"qc192_sched_M{m}.cnfg",
        log_path=log_dir / f"ssd_fc_src2_qc192_M{m}.log",
        raw_dir=raw_dir,
    )


def run_generator(
    *,
    exe: Path,
    config_path: Path,
    matrix_root: Path,
    matrix_id: int,
    snr: str,
    ibex_output_dir: Path,
    lldec_pattern: str,
    rdec_pattern: str,
    log_path: Path,
) -> tuple[Path, Path]:
    before_lldec = snapshot_matching_files(ibex_output_dir, lldec_pattern)
    before_rdec = snapshot_matching_files(ibex_output_dir, rdec_pattern)

    with log_path.open("w", encoding="utf-8") as log_handle:
        subprocess.run(
            [
                str(exe),
                "LDPC",
                str(config_path.resolve()),
                "AWGN",
                snr,
                str(matrix_id),
                str(matrix_root.resolve()),
            ],
            cwd=str(exe.parent),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )

    lldec_path = pick_updated_file(ibex_output_dir, lldec_pattern, before_lldec)
    rdec_path = pick_updated_file(ibex_output_dir, rdec_pattern, before_rdec)
    return lldec_path, rdec_path


def convert_sched_to_verilog(
    *,
    raw_sched_path: Path,
    template_text: str,
    module_name: str,
    output_path: Path,
) -> None:
    raw_lines = raw_sched_path.read_text(encoding="utf-8").splitlines()
    entries = parse_sched_entries(raw_lines, source=raw_sched_path)
    sched_lines, default_bits = normalize_sched_entries(entries)
    verilog_text = build_verilog_from_template(
        template_text,
        module_name=module_name,
        sched_lines=sched_lines,
        default_data_bits=default_bits,
    )
    output_path.write_text(verilog_text, encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Drive QC192 ssd_fc_src2 runs, capture print_hm_192() generated "
            "lldec/rdec sched dumps, and wrap them into Verilog LUT templates."
        )
    )
    parser.add_argument(
        "--exe",
        type=Path,
        default=Path("IBEX/ssd_fc_src2"),
        help="Path to the QC192 executable, e.g. IBEX/ssd_fc_src2.",
    )
    parser.add_argument(
        "--template-cnfg",
        type=Path,
        required=True,
        help="Template .cnfg file. The script only patches line 5 parity size to m*24-4.",
    )
    parser.add_argument(
        "--matrix-root",
        type=Path,
        required=True,
        help=(
            "QC192 matrix root. Supports either "
            "<root>/<M>x<N>/{matrix,fade_matrix,occupied_matrix}/... "
            "or flat <root>/{matrix,fade_matrix,occupied_matrix}/LDPC_Mx(M+178)... files."
        ),
    )
    parser.add_argument(
        "--lldec-template",
        type=Path,
        required=True,
        help="Verilog template for the LLDEC LUT.",
    )
    parser.add_argument(
        "--rdec-template",
        type=Path,
        required=True,
        help="Verilog template for the RDEC LUT.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/qc192_sched_lut"),
        help="Output directory for generated configs, logs, raw sched dumps, and .v files.",
    )
    parser.add_argument(
        "--m-range",
        default="14-45",
        help="M range, e.g. 14-45 or 14,16,18.",
    )
    parser.add_argument(
        "--matrix-id",
        type=int,
        default=1,
        help="Matrix ID passed to the executable.",
    )
    parser.add_argument(
        "--snr",
        default="1",
        help="AWGN SNR argument passed to the executable. Scheduler contents are independent of SNR.",
    )
    parser.add_argument(
        "--qc",
        type=int,
        default=192,
        help="QC size. This script is intended for QC192 and validates --qc=192 by default.",
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Skip M values that do not have matching matrix files/directories under --matrix-root.",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    if args.qc != 192:
        raise ValueError(f"This generator is for QC192. Received --qc={args.qc}")

    exe = args.exe.resolve()
    if not exe.exists():
        raise FileNotFoundError(f"Executable not found: {exe}")

    matrix_root = args.matrix_root.resolve()
    if not matrix_root.exists():
        raise FileNotFoundError(f"Matrix root not found: {matrix_root}")

    template_cnfg_text = args.template_cnfg.read_text(encoding="utf-8")
    lldec_template_text = args.lldec_template.read_text(encoding="utf-8")
    rdec_template_text = args.rdec_template.read_text(encoding="utf-8")
    ibex_output_dir = exe.parent / "output"
    ibex_output_dir.mkdir(parents=True, exist_ok=True)

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    verilog_dir = out_dir / "verilog"
    verilog_dir.mkdir(parents=True, exist_ok=True)

    ms = parse_m_range(args.m_range)
    generated: list[tuple[int, Path, Path]] = []

    for m in ms:
        try:
            ctx = prepare_size_context(
                m=m,
                matrix_root=matrix_root,
                out_dir=out_dir,
                qc=args.qc,
                matrix_id=args.matrix_id,
            )
        except FileNotFoundError:
            if args.skip_missing:
                print(f"[SKIP] No matrix data for M={m} under {matrix_root}", file=sys.stderr)
                continue
            raise

        config_text = build_config_text(template_cnfg_text, m=ctx.m)
        ctx.config_path.write_text(config_text, encoding="utf-8")

        lldec_pattern = f"lldec_sched_{ctx.m}x*ex{args.qc}_w*.txt"
        rdec_pattern = f"rdec_sched_192v_{ctx.m}x*ex{args.qc}_w*.txt"

        print(f"[RUN] M={ctx.m} N={ctx.n} size={ctx.size_name}", file=sys.stderr)
        lldec_raw, rdec_raw = run_generator(
            exe=exe,
            config_path=ctx.config_path,
            matrix_root=ctx.run_matrix_root,
            matrix_id=args.matrix_id,
            snr=args.snr,
            ibex_output_dir=ibex_output_dir,
            lldec_pattern=lldec_pattern,
            rdec_pattern=rdec_pattern,
            log_path=ctx.log_path,
        )

        lldec_raw_copy = ctx.raw_dir / lldec_raw.name
        rdec_raw_copy = ctx.raw_dir / rdec_raw.name
        lldec_raw_copy.write_text(lldec_raw.read_text(encoding="utf-8"), encoding="utf-8")
        rdec_raw_copy.write_text(rdec_raw.read_text(encoding="utf-8"), encoding="utf-8")

        module_m = ctx.m
        lldec_module_name = f"lldec_matrix_m{module_m}_lut"
        rdec_module_name = f"rdec_matrix_m{module_m}_lut"
        lldec_output_path = verilog_dir / f"{lldec_module_name}.v"
        rdec_output_path = verilog_dir / f"{rdec_module_name}.v"

        convert_sched_to_verilog(
            raw_sched_path=lldec_raw_copy,
            template_text=lldec_template_text,
            module_name=lldec_module_name,
            output_path=lldec_output_path,
        )
        convert_sched_to_verilog(
            raw_sched_path=rdec_raw_copy,
            template_text=rdec_template_text,
            module_name=rdec_module_name,
            output_path=rdec_output_path,
        )

        generated.append((ctx.m, lldec_output_path, rdec_output_path))
        print(
            f"[DONE] M={ctx.m} -> {lldec_output_path.name}, {rdec_output_path.name}",
            file=sys.stderr,
        )

    print("\nGenerated LUT files:")
    for m, lldec_path, rdec_path in generated:
        print(f"M={m}: {lldec_path} | {rdec_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
