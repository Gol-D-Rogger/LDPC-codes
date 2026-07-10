#!/usr/bin/env python3
"""Apply OCR restoration code style after clang-format.

clang-format handles indentation and expression joining. This script adds the
project's preferred breathing room after loop/control blocks without touching
raw OCR backups.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE_SUFFIXES = {".cpp", ".h", ".cu"}
CONTROL_RE = re.compile(r"^\s*(for|while|if|switch)\b")
BRACE_ONLY_RE = re.compile(r"^\s*}\s*$")
SKIP_AFTER_RE = re.compile(r"^\s*(}|else\b|catch\b|while\b|;|,)")


def classify_open(line: str) -> str:
    return "control" if CONTROL_RE.match(line) else "other"


def add_control_block_spacing(text: str) -> str:
    lines = text.splitlines()
    out: list[str] = []
    stack: list[str] = []

    for idx, line in enumerate(lines):
        stripped = line.strip()
        closing_control = False

        if BRACE_ONLY_RE.match(line) and stack:
            closing_control = stack.pop() == "control"

        out.append(line)

        if closing_control:
            next_line = lines[idx + 1] if idx + 1 < len(lines) else ""
            if next_line.strip() and not SKIP_AFTER_RE.match(next_line):
                out.append("")

        opens = line.count("{")
        closes = line.count("}")
        if opens:
            first_type = classify_open(line)
            for open_idx in range(opens):
                stack.append(first_type if open_idx == 0 else "other")

        extra_closes = closes - (1 if BRACE_ONLY_RE.match(line) else 0)
        for _ in range(max(0, extra_closes)):
            if stack:
                stack.pop()

    formatted = "\n".join(out).rstrip() + "\n"
    formatted = re.sub(r"\n{3,}", "\n\n", formatted)
    return formatted


def collect_files() -> list[Path]:
    files = [p for p in ROOT.iterdir() if p.is_file() and p.suffix in SOURCE_SUFFIXES]
    files = [p for p in files if p.name != "decoder_gpu.cu"]
    repaired = ROOT / "ocr_restore_layers" / "repaired"
    if repaired.is_dir():
        files.extend(p for p in repaired.iterdir() if p.is_file() and p.suffix in SOURCE_SUFFIXES)
    files = [p for p in files if p.name != "decoder_gpu.cu"]
    return sorted(files)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="report files that would change")
    args = parser.parse_args()

    files = collect_files()

    changed: list[Path] = []
    for path in files:
        original = path.read_text(encoding="utf-8", errors="ignore")
        formatted_proc = subprocess.run(
            ["clang-format", str(path)],
            cwd=ROOT,
            check=True,
            text=True,
            capture_output=True,
        )
        formatted = add_control_block_spacing(formatted_proc.stdout)
        if formatted != original:
            changed.append(path)
            if not args.check:
                path.write_text(formatted, encoding="utf-8")

    if args.check and changed:
        for path in changed:
            print(path.relative_to(ROOT))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
