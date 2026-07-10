#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import zipfile
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


NS_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


class VrefPlotError(RuntimeError):
    pass


def parse_float(v: Any) -> float | None:
    if v is None:
        return None
    s = str(v).strip()
    if not s:
        return None
    try:
        x = float(s)
    except ValueError:
        return None
    if not math.isfinite(x):
        return None
    return x


def col_index_from_cell_ref(ref: str) -> int:
    letters = "".join(ch for ch in str(ref) if ch.isalpha()).upper()
    if not letters:
        return 0
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - ord("A") + 1)
    return idx - 1


def read_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    try:
        raw = zf.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(raw)
    strings: list[str] = []
    for si in root.findall(f"{NS_MAIN}si"):
        parts: list[str] = []
        for t in si.iter(f"{NS_MAIN}t"):
            parts.append(t.text or "")
        strings.append("".join(parts))
    return strings


def read_xlsx_rows(path: Path) -> list[list[Any]]:
    with zipfile.ZipFile(path) as zf:
        shared = read_shared_strings(zf)
        sheet = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))

    rows: list[list[Any]] = []
    for row_el in sheet.iter(f"{NS_MAIN}row"):
        vals: dict[int, Any] = {}
        for cell in row_el.findall(f"{NS_MAIN}c"):
            ref = cell.attrib.get("r", "")
            col = col_index_from_cell_ref(ref)
            ctype = cell.attrib.get("t", "")
            value: Any = None
            if ctype == "inlineStr":
                texts = [t.text or "" for t in cell.iter(f"{NS_MAIN}t")]
                value = "".join(texts)
            else:
                v_el = cell.find(f"{NS_MAIN}v")
                raw = v_el.text if v_el is not None else None
                if ctype == "s" and raw is not None:
                    try:
                        value = shared[int(raw)]
                    except Exception:
                        value = raw
                elif raw is not None:
                    num = parse_float(raw)
                    value = num if num is not None else raw
            vals[col] = value
        if vals:
            max_col = max(vals)
            rows.append([vals.get(i) for i in range(max_col + 1)])
    return rows


def normalize_header(s: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(s or "").strip().lower())


def find_column(headers: list[Any], candidates: set[str]) -> int | None:
    norm = [normalize_header(h) for h in headers]
    for idx, h in enumerate(norm):
        if h in candidates:
            return idx
    return None


def read_curve_from_xlsx(path: Path, *, complete_only: bool, min_fail_cw: int) -> list[tuple[float, float]]:
    rows = read_xlsx_rows(path)
    if not rows:
        return []
    headers = rows[0]
    raw_idx = find_column(headers, {"rawber", "rber"})
    fer_idx = find_column(headers, {"ldpcfer", "fer"})
    fail_idx = find_column(headers, {"failcw", "fail", "failcount"})
    complete_idx = find_column(headers, {"iscomplete", "complete"})
    if raw_idx is None or fer_idx is None:
        raise VrefPlotError(f"Cannot find RAW_BER/LDPC_FER columns in {path}")
    if min_fail_cw > 0 and fail_idx is None:
        raise VrefPlotError(f"--min-fail-cw requires a FAIL_CW column in {path}")

    pts: list[tuple[float, float]] = []
    for row in rows[1:]:
        if complete_only and complete_idx is not None:
            if complete_idx >= len(row) or parse_float(row[complete_idx]) != 1.0:
                continue
        if min_fail_cw > 0:
            fail_cw = parse_float(row[fail_idx] if fail_idx is not None and fail_idx < len(row) else None)
            if fail_cw is None or fail_cw < min_fail_cw:
                continue
        raw = parse_float(row[raw_idx] if raw_idx < len(row) else None)
        fer = parse_float(row[fer_idx] if fer_idx < len(row) else None)
        if raw is None or fer is None or raw <= 0.0 or fer <= 0.0:
            continue
        pts.append((raw, fer))
    pts.sort(key=lambda x: x[0])
    return pts


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    out: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            case_name = str(row.get("case_name", "")).strip()
            if case_name:
                out[case_name] = {str(k): "" if v is None else str(v) for k, v in row.items()}
            progress = str(row.get("progress_xlsx", "")).strip()
            if progress:
                out[str(Path(progress).resolve())] = {str(k): "" if v is None else str(v) for k, v in row.items()}
    return out


def case_name_from_xlsx(path: Path) -> str:
    stem = path.stem
    if stem.endswith("_progress"):
        return stem[: -len("_progress")]
    return stem


def make_vref_label(value: str, fallback: str) -> str:
    vref = str(value or "").strip() or fallback
    if vref.lower().startswith("vref="):
        return vref
    return f"VREF={vref}"


def label_for(path: Path, metadata: dict[str, dict[str, str]]) -> str:
    row = metadata.get(str(path.resolve()))
    if row is None:
        row = metadata.get(case_name_from_xlsx(path))
    fallback = case_name_from_xlsx(path)
    if row is None:
        return make_vref_label("", fallback)
    vref = row.get("vref_label") or row.get("vref_values") or ""
    return make_vref_label(vref, fallback)


def discover_xlsx(input_root: Path) -> list[Path]:
    return sorted(p.resolve() for p in input_root.rglob("*_progress.xlsx") if p.is_file())


def plot_single(path: Path, pts: list[tuple[float, float]], label: str, out_path: Path) -> None:
    import matplotlib.pyplot as plt  # type: ignore

    out_path.parent.mkdir(parents=True, exist_ok=True)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    plt.figure(figsize=(8, 5.5))
    plt.semilogy(xs, ys, marker="o", linewidth=1.4, markersize=4, label=label)
    plt.xlabel("RBER")
    plt.ylabel("FER")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close()


def is_baseline_curve(path: Path, label: str) -> bool:
    case = case_name_from_xlsx(path).lower()
    label_lower = label.lower()
    return case == "baseline" or label_lower.startswith("baseline")


def plot_overlay(curves: list[tuple[Path, str, list[tuple[float, float]]]], out_path: Path) -> None:
    import matplotlib.pyplot as plt  # type: ignore

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(11, 7))
    try:
        cmap = plt.get_cmap("turbo", max(len(curves), 1))
    except ValueError:
        cmap = plt.get_cmap("hsv", max(len(curves), 1))

    for idx, (path, label, pts) in enumerate(curves):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        linestyle = "--" if is_baseline_curve(path, label) else "-"
        plt.semilogy(xs, ys, linestyle=linestyle, color=cmap(idx), marker="o", linewidth=1.1, markersize=3, label=label)

    plt.xlabel("RBER")
    plt.ylabel("FER")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend(fontsize=6, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0)
    plt.tight_layout(rect=(0, 0, 0.72, 1))
    plt.savefig(out_path, dpi=180)
    plt.close()


def parse_args(argv: list[str]) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Plot VREF sweep RBER-FER curves from AutoFER progress xlsx")
    ap.add_argument("--input-root", default="output/vref_sweep", help="Root containing AutoFER *_progress.xlsx files")
    ap.add_argument("--metadata", default="", help="cases.csv generated by submit_vref_sweep.py")
    ap.add_argument("--out-dir", default="", help="Plot output directory (default: <input-root>/plots)")
    ap.add_argument("--complete-only", action="store_true", help="Only plot rows with IsComplete=1")
    ap.add_argument("--min-fail-cw", type=int, default=0, help="Only plot rows with FAIL_CW >= this value (default: 0)")
    return ap.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    input_root = Path(args.input_root).expanduser().resolve()
    if not input_root.exists():
        raise VrefPlotError(f"input root not found: {input_root}")
    metadata_path = Path(args.metadata).expanduser().resolve() if str(args.metadata).strip() else input_root / "cases.csv"
    metadata = load_metadata(metadata_path)
    out_dir = Path(args.out_dir).expanduser().resolve() if str(args.out_dir).strip() else input_root / "plots"

    xlsx_files = discover_xlsx(input_root)
    if not xlsx_files:
        raise VrefPlotError(f"No *_progress.xlsx files found under {input_root}")

    curves: list[tuple[Path, str, list[tuple[float, float]]]] = []
    for xlsx in xlsx_files:
        pts = read_curve_from_xlsx(
            xlsx,
            complete_only=bool(args.complete_only),
            min_fail_cw=max(0, int(args.min_fail_cw)),
        )
        if not pts:
            continue
        label = label_for(xlsx, metadata)
        curves.append((xlsx, label, pts))
        plot_single(xlsx, pts, label, out_dir / "single" / f"{case_name_from_xlsx(xlsx)}.png")

    if not curves:
        raise VrefPlotError("No plottable RBER/FER points found")
    plot_overlay(curves, out_dir / "all_vref_overlay.png")
    print(f"plotted {len(curves)} curve(s) to {out_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except VrefPlotError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
