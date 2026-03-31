#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


PLACEHOLDER_WORD = 0xFFFFFFFFFFF


@dataclass(frozen=True)
class Case:
    m: int
    k: int
    n: int
    user_bytes: int
    parity_bytes: int


def parity_bytes_for_m(m: int) -> int:
    if m == 5:
        return 320
    return m * 64 - 1


def build_case(m: int, k: int) -> Case:
    return Case(
        m=m,
        k=k,
        n=m + k,
        user_bytes=k * 64,
        parity_bytes=parity_bytes_for_m(m),
    )


def build_config_text(case: Case) -> str:
    return "\n".join(
        [
            "0                               // meta data size (B)",
            "512                             // LBA size(B)",
            "0                               // LBA # per DSP CW",
            f"{case.user_bytes:<31d}// userdata size (B)",
            f"{case.parity_bytes:<31d}// parity size (B)",
            f"{case.m:<31d}// row number of base matrix",
            f"{case.n:<31d}// column number of base matrix",
            "512                             // circulant number",
            "5                               // Dense matrix size or rows of E matrix in base matrix",
            "4                               // column weight of base matrix",
            "1                               // maximum simulation number",
            "0                               // maximum error number",
            "FC_RDEC2                        // LDPC decoder",
            "24                              // Fast decoder maximum iteration number",
            "0                               // Iteration to turn on FDEC column skip feature",
            "32                              // Retry decoder maximum iteration number",
            "50                              // Ibex maximum normal iteration number",
            "1024                            // Ibex maximum iteration number",
            "48                              // Ibex post iteration qc threshold",
            "16                              // Ibex post iteration post threshold",
            "1                               // Ibex post enable",
            "0                               // Ibex Syndrome calculate only",
            "1                               // Ibex early terminate disable",
            "3                               // Ibex likelihood width",
            "0.625                           // Retry decoder alpha",
            "1                               // Quantization Mode",
            "8                               // APP/Q bit width",
            "6                               // R bit width",
            "3                               // Fraction bits",
            "7                               // LLR bit width",
            "3                               // LLR fraction bits",
            "MANUAL                          // LLR mode",
            "1                               // Number of NAND soft data bits",
            "0 0.15 -0.15 0.3 -0.3 0.5 -0.5  // VREF values",
            "1.875 -1.875                    // Hard decoding LLR for retry decoder",
            "1060                            // init_synd_wt threshold to skip FDEC",
            "104C11DB7                       // CRC 32 polynomial values",
            "104C11DB7                       // 32-bit LSFR for randomizer",
            "",
        ]
    )


def validate_sched(path: Path, max_col: int) -> None:
    bad: list[tuple[int, int, str]] = []
    for line_no, line in enumerate(path.read_text().splitlines(), 1):
        match = re.search(r"42'h([0-9A-Fa-f]+)", line)
        if not match:
            continue
        word = int(match.group(1), 16)
        if word == PLACEHOLDER_WORD:
            continue
        col = word & 0x7F
        if col > max_col:
            bad.append((line_no, col, line))
    if bad:
        sample = "\n".join(
            f"line {line_no}: col={col}: {text}" for line_no, col, text in bad[:5]
        )
        raise ValueError(
            f"{path} contains col > max_col({max_col}). Sample:\n{sample}"
        )


def validate_sched_lines(lines: list[str], max_col: int) -> None:
    bad: list[tuple[int, int, str]] = []
    for line_no, line in enumerate(lines, 1):
        match = re.search(r"42'h([0-9A-Fa-f]+)", line)
        if not match:
            continue
        word = int(match.group(1), 16)
        if word == PLACEHOLDER_WORD:
            continue
        col = word & 0x7F
        if col > max_col:
            bad.append((line_no, col, line))
    if bad:
        sample = "\n".join(
            f"line {line_no}: col={col}: {text}" for line_no, col, text in bad[:5]
        )
        raise ValueError(f"scheduler lines contain col > max_col({max_col}). Sample:\n{sample}")


def build_verilog_from_template(
    template_text: str,
    module_name: str,
    sched_lines: list[str],
) -> str:
    text = re.sub(
        r"module\s+\w+\s+\(/\*AUTOARG\*/",
        f"module {module_name} (/*AUTOARG*/",
        template_text,
        count=1,
    )

    body_lines = ["always @(*)", "begin", "    case(mmem_radr)"]
    body_lines.extend(f"    {line}" for line in sched_lines)
    body_lines.append("    default  :mmem_rdt=42'h0;")
    body_lines.append("    endcase")
    body_lines.append("end")
    body = "\n".join(body_lines) + "\n"

    text = re.sub(
        r"always @\(\*\)\s*begin.*?^end\s*$",
        body,
        text,
        count=1,
        flags=re.S | re.M,
    )
    return text


def run_case(
    case: Case,
    ibex_dir: Path,
    output_root: Path,
    snr: str,
    template_text: str,
) -> Path:
    config_dir = output_root / "configs"
    log_dir = output_root / "logs"
    verilog_dir = output_root / "verilog"
    config_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    verilog_dir.mkdir(parents=True, exist_ok=True)

    config_path = config_dir / f"Ibex_src_{case.m}x{case.n}.cnfg"
    log_path = log_dir / f"ssd_fc_{case.m}x{case.n}.log"
    generated_sched = ibex_dir / "output" / f"rdec_sched_{case.m}x{case.n}ex512_w4.txt"
    module_name = f"rdec_matrix_m{case.m}k{case.k}_lut"
    verilog_path = verilog_dir / f"{module_name}.v"

    config_path.write_text(build_config_text(case))
    if generated_sched.exists():
        generated_sched.unlink()

    with log_path.open("w") as log_handle:
        subprocess.run(
            ["./ssd_fc", "LDPC", str(config_path.resolve()), "AWGN", snr],
            cwd=ibex_dir,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=True,
        )

    if not generated_sched.exists():
        raise FileNotFoundError(
            f"Expected scheduler not generated: {generated_sched} (see {log_path})"
        )

    validate_sched(generated_sched, case.n - 1)
    sched_lines = generated_sched.read_text().splitlines()
    validate_sched_lines(sched_lines, case.n - 1)
    verilog_text = build_verilog_from_template(template_text, module_name, sched_lines)
    verilog_path.write_text(verilog_text)
    generated_sched.unlink()
    return verilog_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate cropped RDEC_sched files by driving the current "
            "IBEX/src/ldpc_codec.cpp path through ./ssd_fc."
        )
    )
    parser.add_argument(
        "--ibex-dir",
        type=Path,
        default=Path("IBEX"),
        help="IBEX directory containing ssd_fc and output/.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("IBEX/ibex_matrix_flat_13rate/rdec_lut_from_src_k64_67"),
        help="Directory to store generated Verilog LUTs, configs, and logs.",
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=Path("IBEX/rdec_matrix_m10k67_lut.v"),
        help="Verilog LUT template. Only module name and always block are replaced.",
    )
    parser.add_argument(
        "--m-values",
        type=int,
        nargs="*",
        default=list(range(5, 14)),
        help="Row counts to generate. Current IBEX/src supports M<=13.",
    )
    parser.add_argument(
        "--k-values",
        type=int,
        nargs="+",
        default=[67, 66, 65, 64],
        help="Payload-column counts to generate.",
    )
    parser.add_argument(
        "--snr",
        default="4.7",
        help="AWGN SNR argument passed to ssd_fc. Does not affect scheduler contents.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    ibex_dir = args.ibex_dir.resolve()
    ssd_fc_path = ibex_dir / "ssd_fc"
    src_ldpc_path = ibex_dir / "src" / "ldpc_codec.cpp"
    template_path = args.template.resolve()
    if not ssd_fc_path.exists():
        raise FileNotFoundError(f"Missing binary: {ssd_fc_path}")
    if not template_path.exists():
        raise FileNotFoundError(f"Missing template: {template_path}")
    if src_ldpc_path.exists() and ssd_fc_path.stat().st_mtime < src_ldpc_path.stat().st_mtime:
        print(
            f"[src-rdec] WARNING: {ssd_fc_path} is older than {src_ldpc_path}. "
            "Rebuild if you need exact source/binary alignment."
        )
    template_text = template_path.read_text()

    generated: list[Path] = []
    for m in args.m_values:
        if m < 5 or m > 13:
            raise ValueError(f"M={m} is outside current IBEX/src support range [5, 13]")
        for k in args.k_values:
            if k < 64 or k > 67:
                raise ValueError(f"K={k} is outside supported range [64, 67]")
            case = build_case(m, k)
            generated.append(
                run_case(
                    case,
                    ibex_dir,
                    args.output_root.resolve(),
                    args.snr,
                    template_text,
                )
            )

    print(f"[src-rdec] Generated {len(generated)} LUT files under {args.output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
