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
TEMPLATE_ASSIGN_RE = re.compile(
    r"^(?P<indent>\s*)assign\s+matrix_col\[\s*(?P<idx>[0-9]+)\s*\]\s*=\s*(?P<bits>[0-9]+)'h(?P<hex>[0-9A-Fa-f]+)\s*;\s*$"
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
class LldecEntry:
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

    shutil.copyfile(matrix_src, size_matrix_dir / compat_matrix)
    shutil.copyfile(fade_src, size_fade_dir / compat_fade)
    shutil.copyfile(occupied_src, size_occupied_dir / compat_occupied)

    return n, staged_root


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
        config_path=config_dir / f"qc192_lldec_M{m}.cnfg",
        log_path=log_dir / f"ssd_fc_src2_qc192_lldec_M{m}.log",
        raw_dir=raw_dir,
    )


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
        return max(after.items(), key=lambda kv: kv[1])[0]
    raise FileNotFoundError(f"No generated file matching {pattern} under {output_dir}")


def run_generator(
    *,
    exe: Path,
    config_path: Path,
    matrix_root: Path,
    matrix_id: int,
    snr: str,
    ibex_output_dir: Path,
    lldec_pattern: str,
    log_path: Path,
) -> Path:
    before_lldec = snapshot_matching_files(ibex_output_dir, lldec_pattern)
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
    return pick_updated_file(ibex_output_dir, lldec_pattern, before_lldec)


def parse_lldec_entries(raw_lines: list[str], *, source: Path) -> list[LldecEntry]:
    entries: list[LldecEntry] = []
    for line_no, line in enumerate(raw_lines, 1):
        stripped = line.strip()
        if not stripped:
            continue
        match = SCHED_LINE_RE.match(stripped)
        if not match:
            raise ValueError(f"{source}:{line_no}: unsupported LLDEC sched line: {line}")
        entries.append(
            LldecEntry(
                address=int(match.group(2)),
                data_bits=int(match.group(3)),
                data_hex=match.group(4).upper(),
            )
        )
    if not entries:
        raise ValueError(f"{source}: no LLDEC sched entries found")
    entries.sort(key=lambda e: e.address)
    return entries


def build_lldec_verilog_from_template(template_text: str, *, module_name: str, entries: list[LldecEntry]) -> str:
    text, module_subs = re.subn(
        r"(^\s*module\s+)\w+(\s+\(/\*AUTOARG\*/)",
        rf"\1{module_name}\2",
        template_text,
        count=1,
        flags=re.M,
    )
    if module_subs != 1:
        raise ValueError("Failed to replace module line in LLDEC template")

    lines = text.splitlines(keepends=True)
    out_lines: list[str] = []
    insert_at: int | None = None
    indent = ""
    data_bits_from_template: int | None = None

    for line in lines:
        match = TEMPLATE_ASSIGN_RE.match(line.rstrip("\n"))
        if match is None:
            out_lines.append(line)
            continue
        if insert_at is None:
            insert_at = len(out_lines)
            indent = match.group("indent")
            data_bits_from_template = int(match.group("bits"))
        # Drop only assign matrix_col lines.

    if insert_at is None:
        raise ValueError("No assign matrix_col[...] lines found in LLDEC template")

    if data_bits_from_template is None:
        data_bits_from_template = entries[0].data_bits

    idx_width = max(3, max(len(str(entry.address)) for entry in entries))
    assign_lines = [
        f"{indent}assign matrix_col[{entry.address:<{idx_width}}]={data_bits_from_template}'h{entry.data_hex};\n"
        for entry in entries
    ]
    out_lines[insert_at:insert_at] = assign_lines
    return "".join(out_lines)


def convert_lldec_to_verilog(
    *,
    raw_sched_path: Path,
    template_text: str,
    module_name: str,
    output_path: Path,
) -> None:
    raw_lines = raw_sched_path.read_text(encoding="utf-8").splitlines()
    entries = parse_lldec_entries(raw_lines, source=raw_sched_path)
    verilog_text = build_lldec_verilog_from_template(
        template_text,
        module_name=module_name,
        entries=entries,
    )
    output_path.write_text(verilog_text, encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Drive QC192 ssd_fc_src2 runs, capture LLDEC sched dumps, and wrap "
            "them into assign matrix_col[...] template files."
        )
    )
    parser.add_argument("--exe", type=Path, default=Path("IBEX/ssd_fc_src2"))
    parser.add_argument("--template-cnfg", type=Path, required=True)
    parser.add_argument("--matrix-root", type=Path, required=True)
    parser.add_argument("--lldec-template", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("output/qc192_lldec_matrix_lut"))
    parser.add_argument("--m-range", default="14-45")
    parser.add_argument("--matrix-id", type=int, default=1)
    parser.add_argument("--snr", default="1")
    parser.add_argument("--qc", type=int, default=192)
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
    ibex_output_dir = exe.parent / "output"
    ibex_output_dir.mkdir(parents=True, exist_ok=True)

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    verilog_dir = out_dir / "verilog"
    verilog_dir.mkdir(parents=True, exist_ok=True)

    generated: list[tuple[int, Path]] = []
    for m in parse_m_range(args.m_range):
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

        ctx.config_path.write_text(build_config_text(template_cnfg_text, m=ctx.m), encoding="utf-8")
        lldec_pattern = f"lldec_sched_{ctx.m}x*ex{args.qc}_w*.txt"

        print(f"[RUN] M={ctx.m} N={ctx.n} size={ctx.size_name}", file=sys.stderr)
        lldec_raw = run_generator(
            exe=exe,
            config_path=ctx.config_path,
            matrix_root=ctx.run_matrix_root,
            matrix_id=args.matrix_id,
            snr=args.snr,
            ibex_output_dir=ibex_output_dir,
            lldec_pattern=lldec_pattern,
            log_path=ctx.log_path,
        )

        lldec_raw_copy = ctx.raw_dir / lldec_raw.name
        lldec_raw_copy.write_text(lldec_raw.read_text(encoding="utf-8"), encoding="utf-8")

        module_name = f"lldec_matrix_m{ctx.m}_lut"
        output_path = verilog_dir / f"{module_name}.v"
        convert_lldec_to_verilog(
            raw_sched_path=lldec_raw_copy,
            template_text=lldec_template_text,
            module_name=module_name,
            output_path=output_path,
        )
        generated.append((ctx.m, output_path))
        print(f"[DONE] M={ctx.m} -> {output_path.name}", file=sys.stderr)

    print("\nGenerated LLDEC LUT files:")
    for m, path in generated:
        print(f"M={m}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
