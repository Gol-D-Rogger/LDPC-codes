#!/usr/bin/env python3
"""Restore project files from OCR raw layer with conservative cleanup.

This pass is deliberately restore-first: it keeps raw file coverage and removes
mechanical OCR noise without replacing large damaged files with stubs.
"""

from __future__ import annotations

import re
import shutil
import unicodedata
import argparse
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW_DIR = ROOT / "ocr_restore_layers" / "raw"
NORMALIZED_DIR = ROOT / "ocr_restore_layers" / "normalized"
REPAIRED_DIR = ROOT / "ocr_restore_layers" / "repaired"
BACKUP_DIR = ROOT / "ocr_restore_layers" / (
    "previous_short_restore_" + datetime.now().strftime("%Y%m%d_%H%M%S")
)
SOURCE_SUFFIXES = {".h", ".cpp", ".cu"}


TRANSLATIONS = str.maketrans(
    {
        "\ufeff": "",
        "\u200b": "",
        "\u200c": "",
        "\u200d": "",
        "\u2060": "",
        "\u00a0": " ",
        "\u2018": "'",
        "\u2019": "'",
        "\u201c": '"',
        "\u201d": '"',
        "\u2013": "-",
        "\u2014": "-",
        "\u2212": "-",
        "\u00d7": "*",
        "\u2217": "*",
        "\u00f7": "/",
    }
)

TIMESTAMP_RE = re.compile(
    r"(?<![A-Za-z_])\d{1,2}\s*[:.]\s*\d{2}\s*[:.]\s*\d{1,3}(?:\.\d+)?"
)
BROKEN_TIMESTAMP_RE = re.compile(
    r"(?<![A-Za-z_])[A-Z]?\d{1,2}\s*:\s*\d{2,6}(?:\s*:\s*\d{1,3})?"
)
TWO_PART_TIMESTAMP_TAIL_RE = re.compile(
    r"(?<![A-Za-z_])\d{1,2}\s*[:.]\s*\d{1,3}\s*D[0-9A-Za-z]+"
)
TIMESTAMP_TAIL_RE = re.compile(
    r"(?<![A-Za-z_])\d{1,2}\s*[:.]\s*\d{2}\s*[:.]\s*\d{1,3}\s*D[0-9A-Za-z]+"
)
DEVICE_CODE_RE = re.compile(r"(?<![A-Za-z_])D[0-9B][A-Za-z0-9:._<>/-]{1,}")
ATTACHED_DEVICE_RE = re.compile(r"(?<=[0-9])D[0-9B][A-Za-z0-9:._<>/-]{1,}")
STANDALONE_NOISE_RE = re.compile(
    r"^\s*(?:[0-9A-Z]{5,}|[0-9]{1,4}|[A-Z]?:?\d{1,2}[:.]\d{1,3}|[:./|_\\-]+)\s*$"
)


TOKEN_FIXES = (
    ("std: :", "std::"),
    ("logger: :", "logger::"),
    ("decoder: :", "decoder::"),
    (" : :", "::"),
    (": :", "::"),
    ("_FILE_", "__FILE__"),
    ("_LINE_", "__LINE__"),
    ("_syncwarp", "__syncwarp"),
    ("_half2float", "__half2float"),
    ("_float2int_rn", "__float2int_rn"),
    ("float2half", "__float2half"),
    ("__device_inline", "__device__ inline"),
    ("__device_ inline", "__device__ inline"),
    ("__device_ inine", "__device__ inline"),
    ("decoder global info", "decoder_global_info"),
    ("decoder global_info", "decoder_global_info"),
    ("decoder_global info", "decoder_global_info"),
    ("decoder_iput_cw", "decoder_input_cw"),
    ("ldc_decode_kernel", "ldpc_decode_kernel"),
    ("Idpc_decoder_output", "ldpc_decoder_output"),
    ("OptimizedSharedNemory", "OptimizedSharedMemory"),
    ("OptimizedSharedMenory", "OptimizedSharedMemory"),
    ("cw_ shared_mem", "cw_shared_mem"),
    ("cw shared_mem", "cw_shared_mem"),
    ("cw workspace", "cw_workspace"),
    ("shift value", "shift_value"),
    ("layer pre", "layer_pre"),
    ("lane id", "lane_id"),
    ("intidx", "int idx"),
    ("intcol", "int col"),
    ("inti", "int i"),
)

REGEX_TOKEN_FIXES = (
    (re.compile(r"(?<![A-Za-z0-9_])_+device_+(?![A-Za-z0-9_])"), "__device__"),
    (re.compile(r"(?<![A-Za-z0-9_])_+evice_+(?![A-Za-z0-9_])"), "__device__"),
    (re.compile(r"(?<![A-Za-z0-9_])_+global_+(?![A-Za-z0-9_])"), "__global__"),
    (re.compile(r"(?<![A-Za-z0-9_])__gbal_+(?![A-Za-z0-9_])"), "__global__"),
    (re.compile(r"(?<![A-Za-z0-9_])_+shared_+(?![A-Za-z0-9_])"), "__shared__"),
)


def normalize_unicode(text: str) -> str:
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = unicodedata.normalize("NFKC", text)
    return text.translate(TRANSLATIONS)


def drop_non_ascii_noise(line: str) -> str:
    return "".join(ch for ch in line if ch == "\t" or ch == "\n" or ord(ch) < 128)


def clean_line(line: str) -> str:
    line = drop_non_ascii_noise(line.rstrip())
    line = TIMESTAMP_TAIL_RE.sub("", line)
    line = TWO_PART_TIMESTAMP_TAIL_RE.sub("", line)
    line = TIMESTAMP_RE.sub("", line)
    line = BROKEN_TIMESTAMP_RE.sub("", line)
    line = ATTACHED_DEVICE_RE.sub("", line)
    line = DEVICE_CODE_RE.sub("", line)
    for old, new in TOKEN_FIXES:
        line = line.replace(old, new)
    for pattern, replacement in REGEX_TOKEN_FIXES:
        line = pattern.sub(replacement, line)
    line = line.replace("device__global__info", "device_global_info")
    line = line.replace("decoder__global__info", "decoder_global_info")
    line = line.replace("__FILE___", "__FILE__")
    line = line.replace("___half2float", "__half2float")
    line = line.replace("___syncwarp", "__syncwarp")
    line = line.replace("\\|", "\\")
    line = re.sub(r"\b([A-Za-z_][A-Za-z0-9_]*)\s+::\s+", r"\1::", line)
    if line.lstrip().startswith("/ ") and not line.lstrip().startswith("//"):
        line = line.replace("/ ", "// ", 1)
    line = re.sub(r"\s+([,;)\]])", r"\1", line)
    line = re.sub(r"([(\[])\s+", r"\1", line)
    line = re.sub(r"[ \t]{2,}", " ", line).rstrip()
    if not line.strip():
        return ""
    if STANDALONE_NOISE_RE.fullmatch(line):
        return ""
    if re.fullmatch(r"\s*(?:_|__|___|____|[{}();,|\\/:.-]+)\s*", line):
        return ""
    return line


def restore_text(text: str) -> str:
    text = normalize_unicode(text)
    lines = [clean_line(line) for line in text.splitlines()]
    compact = [line for line in lines if line.strip()]
    return "\n".join(compact) + "\n"


def backup_current(files: list[Path]) -> None:
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    for raw_file in files:
        current = ROOT / raw_file.name
        if current.exists():
            shutil.copy2(current, BACKUP_DIR / current.name)


def write_layers(raw_file: Path, restored: str) -> None:
    for directory in (NORMALIZED_DIR, REPAIRED_DIR):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / raw_file.name).write_text(restored, encoding="utf-8")
    (ROOT / raw_file.name).write_text(restored, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="write restored files into normalized/repaired/root layers",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show line-count report without writing files",
    )
    args = parser.parse_args()
    if not args.apply and not args.dry_run:
        parser.error("choose --dry-run to inspect or --apply to write files")

    raw_files = sorted(p for p in RAW_DIR.iterdir() if p.suffix in SOURCE_SUFFIXES)
    if args.apply:
        backup_current(raw_files)
        print(f"backup={BACKUP_DIR}")
    for raw_file in raw_files:
        original = raw_file.read_text(encoding="utf-8", errors="ignore")
        restored = restore_text(original)
        if args.apply:
            write_layers(raw_file, restored)
        raw_lines = original.count("\n") + (0 if original.endswith("\n") else 1)
        out_lines = restored.count("\n")
        print(f"{raw_file.name}: raw_lines={raw_lines} restored_lines={out_lines}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
