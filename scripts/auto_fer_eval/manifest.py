from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional


@dataclass
class RunRecord:
    axis_type: str
    axis_value: float
    log_path: str
    cmd: list[str]
    raw_ber: Optional[float] = None
    theo_rber: Optional[float] = None
    ldpc_fer: Optional[float] = None
    fail_cw: Optional[int] = None
    total_packets: Optional[int] = None
    retry_dec_avg_iter: Optional[float] = None


def load_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"runs": []}
    return json.loads(path.read_text(encoding="utf-8"))


def save_manifest(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def upsert_run(manifest: dict[str, Any], rec: RunRecord) -> None:
    runs: list[dict[str, Any]] = manifest.setdefault("runs", [])

    matches = [
        r
        for r in runs
        if r.get("axis_type") == rec.axis_type and _axis_equal(rec.axis_type, r.get("axis_value"), rec.axis_value)
    ]
    if matches:
        _update_run_dict(matches[0], rec)
        for extra in matches[1:]:
            runs.remove(extra)
        return

    runs.append(_record_to_dict(rec))


def find_run(manifest: dict[str, Any], axis_type: str, axis_value: float) -> Optional[dict[str, Any]]:
    for r in manifest.get("runs", []):
        if r.get("axis_type") == axis_type and _axis_equal(axis_type, r.get("axis_value"), axis_value):
            return r
    return None


def _update_run_dict(dst: dict[str, Any], rec: RunRecord) -> None:
    dst.update(_record_to_dict(rec))


def _record_to_dict(rec: RunRecord) -> dict[str, Any]:
    return {
        "axis_type": rec.axis_type,
        "axis_value": float(rec.axis_value),
        "log_path": rec.log_path,
        "cmd": rec.cmd,
        "raw_ber": rec.raw_ber,
        "theo_rber": rec.theo_rber,
        "ldpc_fer": rec.ldpc_fer,
        "fail_cw": rec.fail_cw,
        "total_packets": rec.total_packets,
        "retry_dec_avg_iter": rec.retry_dec_avg_iter,
    }


def _axis_equal(axis_type: str, a: Any, b: Any) -> bool:
    if a is None or b is None:
        return False
    if axis_type == "k":
        return int(round(float(a))) == int(round(float(b)))
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-9)
