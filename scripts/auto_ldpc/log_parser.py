"""
Log parsing helpers for LDPC automation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence


STATISTICS_RE = re.compile(r"^\[STATISTICS\]\s*(?P<key>[^:]+):\s*(?P<value>[-+0-9.eE]+)")


@dataclass
class LogParseResult:
    matrix_id: str
    model: str
    snr: float
    metrics: Dict[str, float]
    fer_key: str
    fer_value: Optional[float]
    log_path: Path


def _extract_statistics(lines: Iterable[str], wanted_keys: Sequence[str]) -> Dict[str, float]:
    stats: Dict[str, float] = {}
    wanted = {key: key for key in wanted_keys}
    for raw_line in lines:
        match = STATISTICS_RE.match(raw_line.strip())
        if not match:
            continue
        key = match.group("key").strip()
        if key not in wanted:
            continue
        try:
            stats[key] = float(match.group("value"))
        except ValueError:
            continue
    return stats


def parse_log_file(
    log_path: Path,
    matrix_id: str,
    model: str,
    snr: float,
    stat_keys: Sequence[str],
    fer_key: str,
) -> Optional[LogParseResult]:
    if not log_path.exists():
        return None
    try:
        with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
    except OSError:
        return None
    metrics = _extract_statistics(lines, stat_keys)
    fer_value = metrics.get(fer_key)
    return LogParseResult(
        matrix_id=matrix_id,
        model=model,
        snr=snr,
        metrics=metrics,
        fer_key=fer_key,
        fer_value=fer_value,
        log_path=log_path,
    )


MATRIX_DIR_RE = re.compile(r"matrix(\d+)")
SNR_FILE_RE = re.compile(r"snr[_-]?(-?\d+(?:\.\d+)?)", re.IGNORECASE)


def _scan_matrix_dirs(base_dir: Path) -> Iterator[Path]:
    if not base_dir.exists():
        return iter(())
    for entry in sorted(base_dir.iterdir()):
        if entry.is_dir() and MATRIX_DIR_RE.match(entry.name):
            yield entry


def _detect_matrix_id(dir_path: Path) -> Optional[str]:
    match = MATRIX_DIR_RE.match(dir_path.name)
    if match:
        return match.group(1)
    return None


def _detect_snr(log_path: Path) -> Optional[float]:
    match = SNR_FILE_RE.search(log_path.name)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def scan_logs(
    log_root: Path,
    model: str,
    stat_keys: Sequence[str],
    fer_key: str,
    snr_override: Optional[float] = None,
) -> List[LogParseResult]:
    results: List[LogParseResult] = []
    for matrix_dir in _scan_matrix_dirs(log_root):
        matrix_id = _detect_matrix_id(matrix_dir)
        if not matrix_id:
            continue
        for pattern in ("*.log", "*.out", "*.txt"):
            for log_path in sorted(matrix_dir.glob(pattern)):
                snr_value = snr_override if snr_override is not None else _detect_snr(log_path)
                if snr_value is None:
                    continue
                parsed = parse_log_file(
                    log_path=log_path,
                    matrix_id=matrix_id,
                    model=model,
                    snr=snr_value,
                    stat_keys=stat_keys,
                    fer_key=fer_key,
                )
                if parsed:
                    results.append(parsed)
    return results
