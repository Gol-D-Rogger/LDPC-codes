#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


SIZE_RE = re.compile(r"^([0-9]+)x([0-9]+)$")
RANKING_RE = re.compile(r"^ranking_snr(.+)\.csv$")


@dataclass(frozen=True)
class CsvFileInfo:
    path: Path
    size_name: str
    m: int
    n: int
    rel_from_size: Path
    kind: str
    snr_tag: str
    snr_value: float | None
    realtime: bool


@dataclass(frozen=True)
class BestMatrixCandidate:
    point_idx: int
    point_label: str
    size_name: str
    m: int
    n: int
    snr_tag: str
    snr_value: float | None
    matrix_id: int
    ldpc_fer: float
    raw_ber: float | None
    ranking_kind: str
    source_ranking_csv: Path


def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def log(msg: str) -> None:
    print(f"[{_ts()}] {msg}")


def parse_float(v: str | None) -> float | None:
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


def parse_int(v: str | None) -> int | None:
    if v is None:
        return None
    s = str(v).strip()
    if not s:
        return None
    try:
        if "." in s:
            return int(float(s))
        return int(s)
    except ValueError:
        return None


def decode_snr_tag(tag: str) -> float | None:
    t = str(tag).strip()
    if not t:
        return None
    sign = 1.0
    if t.startswith("m"):
        sign = -1.0
        t = t[1:]
    t = t.replace("p", ".")
    try:
        return sign * float(t)
    except ValueError:
        return None


def parse_size_context(csv_path: Path) -> tuple[str, int, int, Path] | None:
    # Expect: <root>/<size>/summary/*.csv
    # or     <root>/<size>/summary/realtime/*.csv
    summary_dir: Path | None = None
    if csv_path.parent.name == "summary":
        summary_dir = csv_path.parent
    elif csv_path.parent.name == "realtime" and csv_path.parent.parent.name == "summary":
        summary_dir = csv_path.parent.parent
    if summary_dir is None:
        return None

    size_dir = summary_dir.parent
    m = SIZE_RE.match(size_dir.name)
    if m is None:
        return None
    size_name = size_dir.name
    m_val, n_val = int(m.group(1)), int(m.group(2))
    rel_from_size = csv_path.relative_to(size_dir)
    return size_name, m_val, n_val, rel_from_size


def classify_file(rel_from_size: Path, file_name: str) -> tuple[str, str, float | None, bool]:
    realtime = "realtime" in rel_from_size.parts
    m = RANKING_RE.match(file_name)
    if m is not None:
        tag = m.group(1)
        if realtime:
            return "ranking_realtime", tag, decode_snr_tag(tag), True
        return "ranking_final", tag, decode_snr_tag(tag), False
    if file_name == "rber_fer_points.csv":
        return "rber_fer_points", "", None, realtime
    if file_name == "fer_ranking.csv":
        return "fer_ranking", "", None, realtime
    return "other", "", None, realtime


def discover_summary_csv(input_root: Path, include_realtime: bool) -> list[CsvFileInfo]:
    out: list[CsvFileInfo] = []
    for path in sorted(input_root.rglob("*.csv")):
        parsed = parse_size_context(path)
        if parsed is None:
            continue
        size_name, m_val, n_val, rel = parsed
        if (not include_realtime) and ("realtime" in rel.parts):
            continue
        kind, snr_tag, snr_value, realtime = classify_file(rel, path.name)
        out.append(
            CsvFileInfo(
                path=path.resolve(),
                size_name=size_name,
                m=m_val,
                n=n_val,
                rel_from_size=rel,
                kind=kind,
                snr_tag=snr_tag,
                snr_value=snr_value,
                realtime=realtime,
            )
        )
    return out


def copy_by_size(files: list[CsvFileInfo], out_dir: Path) -> int:
    copied = 0
    by_size_root = out_dir / "by_size"
    for info in files:
        dst = (by_size_root / info.size_name / info.rel_from_size).resolve()
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(info.path, dst)
        copied += 1
    return copied


def write_file_index(files: list[CsvFileInfo], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "size",
                "m",
                "n",
                "kind",
                "realtime",
                "snr_tag",
                "snr_value",
                "rel_from_size",
                "source_csv",
            ]
        )
        for x in sorted(files, key=lambda f: (f.m, f.n, f.kind, str(f.rel_from_size))):
            w.writerow(
                [
                    x.size_name,
                    x.m,
                    x.n,
                    x.kind,
                    int(x.realtime),
                    x.snr_tag,
                    "" if x.snr_value is None else x.snr_value,
                    str(x.rel_from_size),
                    str(x.path),
                ]
            )


def merge_csv_rows(files: list[CsvFileInfo], out_csv: Path) -> tuple[int, int]:
    meta_cols = [
        "size",
        "m",
        "n",
        "kind",
        "realtime",
        "snr_tag",
        "snr_value",
        "source_csv",
    ]
    data_cols: list[str] = []
    data_seen: set[str] = set()
    merged_rows: list[dict[str, str]] = []
    file_count = 0

    for info in files:
        try:
            text = info.path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        lines = text.splitlines()
        if not lines:
            continue
        reader = csv.DictReader(lines)
        if not reader.fieldnames:
            continue
        file_count += 1
        for fn in reader.fieldnames:
            key = str(fn)
            if key not in data_seen:
                data_seen.add(key)
                data_cols.append(key)

        for row in reader:
            rec: dict[str, str] = {
                "size": info.size_name,
                "m": str(info.m),
                "n": str(info.n),
                "kind": info.kind,
                "realtime": str(int(info.realtime)),
                "snr_tag": info.snr_tag,
                "snr_value": "" if info.snr_value is None else str(info.snr_value),
                "source_csv": str(info.path),
            }
            for k, v in row.items():
                rec[str(k)] = "" if v is None else str(v)
            merged_rows.append(rec)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=meta_cols + data_cols)
        w.writeheader()
        for r in merged_rows:
            w.writerow(r)
    return file_count, len(merged_rows)


def read_best_row_from_ranking(info: CsvFileInfo) -> BestMatrixCandidate | None:
    try:
        text = info.path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return None
    lines = text.splitlines()
    if not lines:
        return None

    reader = csv.DictReader(lines)
    rows: list[tuple[int, float, int, float | None]] = []
    for r in reader:
        matrix_id = parse_int(r.get("matrix_id"))
        ldpc_fer = parse_float(r.get("ldpc_fer"))
        # Ignore unfinished/invalid rows (e.g., FER=0 placeholders).
        if matrix_id is None or ldpc_fer is None or ldpc_fer <= 0.0:
            continue
        rank = parse_int(r.get("rank"))
        if rank is None:
            rank = 10**9
        raw_ber = parse_float(r.get("raw_ber"))
        rows.append((rank, ldpc_fer, matrix_id, raw_ber))

    if not rows:
        return None
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    _, best_fer, best_mid, best_raw = rows[0]
    return BestMatrixCandidate(
        point_idx=0,
        point_label="",
        size_name=info.size_name,
        m=info.m,
        n=info.n,
        snr_tag=info.snr_tag,
        snr_value=info.snr_value,
        matrix_id=best_mid,
        ldpc_fer=best_fer,
        raw_ber=best_raw,
        ranking_kind=info.kind,
        source_ranking_csv=info.path,
    )


def _ranking_sort_key(x: CsvFileInfo) -> tuple[int, float, str]:
    if x.snr_value is not None:
        return (0, float(x.snr_value), x.snr_tag)
    return (1, 0.0, x.snr_tag)


def select_best_matrix_candidates_by_size(files: list[CsvFileInfo]) -> dict[str, dict[int, BestMatrixCandidate]]:
    # return: size -> {point_idx(1/2/3) -> candidate}
    point_labels = {1: "firstPointMatrix", 2: "secondPointMatrix", 3: "thirdPointMatrix"}
    ranking_files = [x for x in files if x.kind in {"ranking_final", "ranking_realtime"}]
    if not ranking_files:
        return {}

    by_size: dict[str, list[CsvFileInfo]] = defaultdict(list)
    for x in ranking_files:
        by_size[x.size_name].append(x)

    out: dict[str, dict[int, BestMatrixCandidate]] = {}
    for size_name, flist in by_size.items():
        # Same size+snr may have both final/realtime; prefer final.
        by_snr: dict[str, list[CsvFileInfo]] = defaultdict(list)
        for fi in flist:
            snr_key = fi.snr_tag if fi.snr_tag else f"snr_{fi.snr_value}"
            by_snr[snr_key].append(fi)

        chosen_per_snr: list[CsvFileInfo] = []
        for _, sgroup in by_snr.items():
            sgroup_sorted = sorted(
                sgroup,
                key=lambda x: (
                    0 if x.kind == "ranking_final" else 1,
                    str(x.path),
                ),
            )
            chosen_per_snr.append(sgroup_sorted[0])

        ordered = sorted(chosen_per_snr, key=_ranking_sort_key)
        size_points: dict[int, BestMatrixCandidate] = {}
        for idx, info in enumerate(ordered[:3], start=1):
            best = read_best_row_from_ranking(info)
            if best is None:
                continue
            size_points[idx] = BestMatrixCandidate(
                point_idx=idx,
                point_label=point_labels[idx],
                size_name=best.size_name,
                m=best.m,
                n=best.n,
                snr_tag=best.snr_tag,
                snr_value=best.snr_value,
                matrix_id=best.matrix_id,
                ldpc_fer=best.ldpc_fer,
                raw_ber=best.raw_ber,
                ranking_kind=best.ranking_kind,
                source_ranking_csv=best.source_ranking_csv,
            )
        if size_points:
            out[size_name] = size_points
    return out


def select_deepest_point_candidate_by_size(files: list[CsvFileInfo]) -> dict[str, BestMatrixCandidate]:
    ranking_files = [x for x in files if x.kind in {"ranking_final", "ranking_realtime"}]
    if not ranking_files:
        return {}

    by_size: dict[str, list[CsvFileInfo]] = defaultdict(list)
    for x in ranking_files:
        by_size[x.size_name].append(x)

    out: dict[str, BestMatrixCandidate] = {}
    for size_name, flist in by_size.items():
        # Same size+snr may have both final/realtime; prefer final.
        by_snr: dict[str, list[CsvFileInfo]] = defaultdict(list)
        for fi in flist:
            snr_key = fi.snr_tag if fi.snr_tag else f"snr_{fi.snr_value}"
            by_snr[snr_key].append(fi)

        chosen_per_snr: list[CsvFileInfo] = []
        for _, sgroup in by_snr.items():
            sgroup_sorted = sorted(
                sgroup,
                key=lambda x: (
                    0 if x.kind == "ranking_final" else 1,
                    str(x.path),
                ),
            )
            chosen_per_snr.append(sgroup_sorted[0])

        # Deepest first (largest SNR).
        ordered_desc = sorted(chosen_per_snr, key=_ranking_sort_key, reverse=True)
        for info in ordered_desc:
            best = read_best_row_from_ranking(info)
            if best is None:
                continue
            out[size_name] = BestMatrixCandidate(
                point_idx=0,
                point_label="deepestPointMatrix",
                size_name=best.size_name,
                m=best.m,
                n=best.n,
                snr_tag=best.snr_tag,
                snr_value=best.snr_value,
                matrix_id=best.matrix_id,
                ldpc_fer=best.ldpc_fer,
                raw_ber=best.raw_ber,
                ranking_kind=best.ranking_kind,
                source_ranking_csv=best.source_ranking_csv,
            )
            break
    return out


def find_matrix_files(size_dir: Path, subdir: str, matrix_id: int) -> list[Path]:
    src = size_dir / subdir
    if not src.exists():
        return []

    mid = int(matrix_id)
    if subdir == "matrix":
        patterns = [
            f"*_H_{mid}_1.txt",
            f"*_H_1_{mid}.txt",
            f"*_{mid}_1.txt",
            f"*_1_{mid}.txt",
        ]
    elif subdir == "fade_matrix":
        patterns = [
            f"*_fade_{mid}_1.txt",
            f"*_fade_1_{mid}.txt",
            f"*fade*_{mid}_1.txt",
            f"*fade*_1_{mid}.txt",
        ]
    elif subdir == "occupied_matrix":
        patterns = [
            f"*_occupied_{mid}_1.txt",
            f"*_occupied_1_{mid}.txt",
            f"*occupied*_{mid}_1.txt",
            f"*occupied*_1_{mid}.txt",
        ]
    else:
        return []

    seen: set[str] = set()
    out: list[Path] = []
    for ptn in patterns:
        for p in sorted(src.glob(ptn)):
            key = str(p.resolve())
            if key in seen:
                continue
            seen.add(key)
            out.append(p.resolve())
    return out


def collect_best_matrix_files(
    *,
    files: list[CsvFileInfo],
    matrix_root: Path,
    best_matrix_dir: Path,
    summary_csv_out: Path,
) -> None:
    selected_by_size = select_best_matrix_candidates_by_size(files)
    if not selected_by_size:
        log("no ranking csv found, skip best matrix collection")
        return

    best_matrix_dir.mkdir(parents=True, exist_ok=True)
    point_labels = {1: "firstPointMatrix", 2: "secondPointMatrix", 3: "thirdPointMatrix"}
    rows: list[dict[str, str]] = []

    for size_name in sorted(selected_by_size.keys()):
        size_selected = selected_by_size[size_name]
        for idx in (1, 2, 3):
            label = point_labels[idx]
            # Flatten all sizes into point-level folders:
            #   firstPointMatrix/{matrix,fade_matrix,occupied_matrix}
            #   secondPointMatrix/{...}
            #   thirdPointMatrix/{...}
            base_dir = best_matrix_dir / label
            (base_dir / "matrix").mkdir(parents=True, exist_ok=True)
            (base_dir / "fade_matrix").mkdir(parents=True, exist_ok=True)
            (base_dir / "occupied_matrix").mkdir(parents=True, exist_ok=True)

            cand = size_selected.get(idx)
            if cand is None:
                rows.append(
                    {
                        "point_label": label,
                        "size": size_name,
                        "m": "",
                        "n": "",
                        "matrix_id": "",
                        "snr_tag": "",
                        "snr_value": "",
                        "ldpc_fer": "",
                        "raw_ber": "",
                        "ranking_kind": "",
                        "source_ranking_csv": "",
                        "copied_matrix_files": "",
                        "copied_fade_files": "",
                        "copied_occupied_files": "",
                    }
                )
                continue

            size_dir = matrix_root / cand.size_name
            if not size_dir.exists():
                log(f"[WARN] matrix size dir not found: {size_dir}")

            copied: dict[str, list[str]] = {"matrix": [], "fade_matrix": [], "occupied_matrix": []}
            for sub in ("matrix", "fade_matrix", "occupied_matrix"):
                src_files = find_matrix_files(size_dir=size_dir, subdir=sub, matrix_id=cand.matrix_id)
                if not src_files:
                    log(f"[WARN] {size_name}/{label}: no {sub} file found for id={cand.matrix_id}")
                    continue
                for src in src_files:
                    # Add size/id prefix to avoid collisions across different M/N folders.
                    dst_name = f"{size_name}__id{cand.matrix_id}__{src.name}"
                    dst = (base_dir / sub / dst_name).resolve()
                    shutil.copy2(src, dst)
                    copied[sub].append(str(dst))

            rows.append(
                {
                    "point_label": label,
                    "size": cand.size_name,
                    "m": str(cand.m),
                    "n": str(cand.n),
                    "matrix_id": str(cand.matrix_id),
                    "snr_tag": cand.snr_tag,
                    "snr_value": "" if cand.snr_value is None else str(cand.snr_value),
                    "ldpc_fer": str(cand.ldpc_fer),
                    "raw_ber": "" if cand.raw_ber is None else str(cand.raw_ber),
                    "ranking_kind": cand.ranking_kind,
                    "source_ranking_csv": str(cand.source_ranking_csv),
                    "copied_matrix_files": ";".join(copied["matrix"]),
                    "copied_fade_files": ";".join(copied["fade_matrix"]),
                    "copied_occupied_files": ";".join(copied["occupied_matrix"]),
                }
            )
            log(
                f"{size_name}/{label}: matrix_id={cand.matrix_id} "
                f"snr={cand.snr_tag} fer={cand.ldpc_fer:g} source={cand.ranking_kind}"
            )

    summary_csv_out.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv_out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "point_label",
                "size",
                "m",
                "n",
                "matrix_id",
                "snr_tag",
                "snr_value",
                "ldpc_fer",
                "raw_ber",
                "ranking_kind",
                "source_ranking_csv",
                "copied_matrix_files",
                "copied_fade_files",
                "copied_occupied_files",
            ],
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)
    log(f"wrote best-matrix summary: {summary_csv_out}")


def collect_deepest_point_matrix_files(
    *,
    files: list[CsvFileInfo],
    matrix_root: Path,
    deepest_matrix_dir: Path,
    summary_csv_out: Path,
) -> None:
    selected = select_deepest_point_candidate_by_size(files)
    if not selected:
        log("no ranking csv found, skip deepest-point matrix collection")
        return

    base_dir = deepest_matrix_dir / "deepestPointMatrix"
    (base_dir / "matrix").mkdir(parents=True, exist_ok=True)
    (base_dir / "fade_matrix").mkdir(parents=True, exist_ok=True)
    (base_dir / "occupied_matrix").mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for size_name in sorted(selected.keys()):
        cand = selected[size_name]
        size_dir = matrix_root / cand.size_name
        if not size_dir.exists():
            log(f"[WARN] matrix size dir not found: {size_dir}")

        copied: dict[str, list[str]] = {"matrix": [], "fade_matrix": [], "occupied_matrix": []}
        for sub in ("matrix", "fade_matrix", "occupied_matrix"):
            src_files = find_matrix_files(size_dir=size_dir, subdir=sub, matrix_id=cand.matrix_id)
            if not src_files:
                log(f"[WARN] {size_name}/deepestPointMatrix: no {sub} file found for id={cand.matrix_id}")
                continue
            for src in src_files:
                dst_name = f"{size_name}__id{cand.matrix_id}__{src.name}"
                dst = (base_dir / sub / dst_name).resolve()
                shutil.copy2(src, dst)
                copied[sub].append(str(dst))

        rows.append(
            {
                "point_label": "deepestPointMatrix",
                "size": cand.size_name,
                "m": str(cand.m),
                "n": str(cand.n),
                "matrix_id": str(cand.matrix_id),
                "snr_tag": cand.snr_tag,
                "snr_value": "" if cand.snr_value is None else str(cand.snr_value),
                "ldpc_fer": str(cand.ldpc_fer),
                "raw_ber": "" if cand.raw_ber is None else str(cand.raw_ber),
                "ranking_kind": cand.ranking_kind,
                "source_ranking_csv": str(cand.source_ranking_csv),
                "copied_matrix_files": ";".join(copied["matrix"]),
                "copied_fade_files": ";".join(copied["fade_matrix"]),
                "copied_occupied_files": ";".join(copied["occupied_matrix"]),
            }
        )
        log(
            f"{size_name}/deepestPointMatrix: matrix_id={cand.matrix_id} "
            f"snr={cand.snr_tag} fer={cand.ldpc_fer:g} source={cand.ranking_kind}"
        )

    summary_csv_out.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv_out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "point_label",
                "size",
                "m",
                "n",
                "matrix_id",
                "snr_tag",
                "snr_value",
                "ldpc_fer",
                "raw_ber",
                "ranking_kind",
                "source_ranking_csv",
                "copied_matrix_files",
                "copied_fade_files",
                "copied_occupied_files",
            ],
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)
    log(f"wrote deepest-point summary: {summary_csv_out}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="One-click organizer for batch_sim summary CSV files")
    ap.add_argument(
        "--input-root",
        default="output/qc192_batch_simple_out",
        help="Root dir containing per-size folders (e.g. 19xNN/summary/...)",
    )
    ap.add_argument(
        "--out-dir",
        default="",
        help="Output dir for organized data (default: <input-root>/summary_organized)",
    )
    ap.add_argument(
        "--matrix-root",
        default="",
        help="Matrix root for collecting best matrix files (optional, e.g. IBEX/ibex_matrix)",
    )
    ap.add_argument(
        "--best-matrix-dir",
        default="",
        help="Output dir for firstPointMatrix/secondPointMatrix/thirdPointMatrix "
        "(default: <out-dir>/best_point_matrices)",
    )
    ap.add_argument(
        "--collect-deepest-per-size",
        action="store_true",
        help="Also collect one best matrix per size at the deepest currently available ranking point",
    )
    ap.add_argument("--no-realtime", action="store_true", help="Exclude summary/realtime/*.csv")
    ap.add_argument("--no-copy", action="store_true", help="Do not copy raw CSVs to by_size/, only generate merged tables")
    return ap.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    input_root = Path(args.input_root).expanduser().resolve()
    if not input_root.exists():
        raise SystemExit(f"[ERROR] input root not found: {input_root}")

    out_dir = Path(args.out_dir).expanduser().resolve() if str(args.out_dir).strip() else (input_root / "summary_organized")
    out_dir.mkdir(parents=True, exist_ok=True)
    include_realtime = not bool(args.no_realtime)

    files = discover_summary_csv(input_root=input_root, include_realtime=include_realtime)
    if not files:
        raise SystemExit(f"[ERROR] no summary csv found under: {input_root}")

    log(f"discovered {len(files)} summary csv file(s)")
    if not args.no_copy:
        copied = copy_by_size(files, out_dir)
        log(f"copied {copied} file(s) to: {out_dir / 'by_size'}")

    index_csv = out_dir / "merged" / "summary_file_index.csv"
    write_file_index(files, index_csv)
    log(f"wrote file index: {index_csv}")

    groups: dict[str, list[CsvFileInfo]] = {
        "all_summary_rows": files,
        "ranking_final_all": [x for x in files if x.kind == "ranking_final"],
        "ranking_realtime_all": [x for x in files if x.kind == "ranking_realtime"],
        "rber_fer_points_all": [x for x in files if x.kind == "rber_fer_points"],
        "fer_ranking_all": [x for x in files if x.kind == "fer_ranking"],
        "other_summary_all": [x for x in files if x.kind == "other"],
    }

    for name, group in groups.items():
        if not group:
            continue
        out_csv = out_dir / "merged" / f"{name}.csv"
        n_files, n_rows = merge_csv_rows(group, out_csv)
        log(f"merged {name}: files={n_files}, rows={n_rows}, out={out_csv}")

    if str(args.matrix_root).strip():
        matrix_root = Path(args.matrix_root).expanduser().resolve()
        if not matrix_root.exists():
            raise SystemExit(f"[ERROR] matrix root not found: {matrix_root}")
        best_matrix_dir = (
            Path(args.best_matrix_dir).expanduser().resolve()
            if str(args.best_matrix_dir).strip()
            else (out_dir / "best_point_matrices")
        )
        collect_best_matrix_files(
            files=files,
            matrix_root=matrix_root,
            best_matrix_dir=best_matrix_dir,
            summary_csv_out=out_dir / "merged" / "best_matrix_by_point.csv",
        )
        if bool(args.collect_deepest_per_size):
            collect_deepest_point_matrix_files(
                files=files,
                matrix_root=matrix_root,
                deepest_matrix_dir=best_matrix_dir,
                summary_csv_out=out_dir / "merged" / "best_matrix_deepest_by_size.csv",
            )

    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
