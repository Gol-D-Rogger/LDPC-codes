from __future__ import annotations

import json
import math
import os
import tempfile
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
    # Strip in-memory index before serializing
    serializable = {k: v for k, v in data.items() if k != "_idx"}
    content = json.dumps(serializable, indent=2, sort_keys=True)
    # Atomic write: tmp file + rename to avoid corruption on crash
    fd, tmp = tempfile.mkstemp(
        dir=str(path.parent), suffix=".tmp", prefix=".manifest_"
    )
    try:
        os.write(fd, content.encode("utf-8"))
        os.fsync(fd)
        os.close(fd)
        os.replace(tmp, str(path))
    except BaseException:
        try:
            os.close(fd)
        except OSError:
            pass
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _ensure_index(manifest: dict[str, Any]) -> dict[tuple, int]:
    """Return (and lazily build) an O(1) lookup index for runs."""
    idx = manifest.get("_idx")
    if idx is not None:
        return idx  # type: ignore[return-value]
    idx = {}
    for i, r in enumerate(manifest.get("runs", [])):
        at = r.get("axis_type")
        av = r.get("axis_value")
        if at is not None and av is not None:
            idx[_axis_key(at, av)] = i
    manifest["_idx"] = idx
    return idx


def upsert_run(manifest: dict[str, Any], rec: RunRecord) -> None:
    runs: list[dict[str, Any]] = manifest.setdefault("runs", [])
    idx = _ensure_index(manifest)
    key = _axis_key(rec.axis_type, rec.axis_value)

    pos = idx.get(key)
    if pos is not None and pos < len(runs):
        _update_run_dict(runs[pos], rec)
    else:
        idx[key] = len(runs)
        runs.append(_record_to_dict(rec))


def find_run(manifest: dict[str, Any], axis_type: str, axis_value: float) -> Optional[dict[str, Any]]:
    idx = _ensure_index(manifest)
    key = _axis_key(axis_type, axis_value)
    pos = idx.get(key)
    if pos is not None:
        runs = manifest.get("runs", [])
        if pos < len(runs):
            return runs[pos]
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


def _axis_key(axis_type: str, axis_value: Any) -> tuple[str, int | float]:
    """Canonical key for axis-based lookups."""
    if axis_type == "k":
        return (axis_type, int(round(float(axis_value))))
    return (axis_type, round(float(axis_value), 9))
