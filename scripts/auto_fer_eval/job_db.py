from __future__ import annotations

import json
import os
import tempfile
import time
from pathlib import Path
from typing import Any, Optional


def _now() -> float:
    return time.time()


def _safe_load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"jobs": []}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {"jobs": []}
    if not isinstance(data, dict):
        return {"jobs": []}
    jobs = data.get("jobs", [])
    if not isinstance(jobs, list):
        return {"jobs": []}
    return {"jobs": jobs}


class JobDB:
    """
    A lightweight persistent job status DB.

    - Stores one record per submitted job (keyed by backend+job_id).
    - Appends state transitions so we can audit what happened after the fact.
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        self.data = _safe_load_json(path)
        self._by_key: dict[str, dict[str, Any]] = {}
        for rec in self.data.get("jobs", []):
            if isinstance(rec, dict) and isinstance(rec.get("key"), str):
                self._by_key[rec["key"]] = rec
        self._dirty = False

    def upsert(self, rec: dict[str, Any]) -> None:
        key = rec.get("key")
        if not isinstance(key, str) or not key:
            raise ValueError("JobDB.upsert requires rec['key'] (non-empty string).")

        existing = self._by_key.get(key)
        if existing is None:
            rec2 = dict(rec)
            rec2.setdefault("history", [])
            rec2.setdefault("last_state", "SUBMITTED")
            rec2.setdefault("submitted_at", _now())
            rec2.setdefault("last_update", rec2.get("submitted_at"))
            self.data["jobs"].append(rec2)
            self._by_key[key] = rec2
            self._dirty = True
            return

        # Merge (do not delete existing fields).
        for k, v in rec.items():
            if v is None:
                continue
            existing[k] = v
        existing["last_update"] = _now()
        self._dirty = True

    def set_state(self, key: str, state: str, *, t: Optional[float] = None) -> None:
        if not key:
            return
        rec = self._by_key.get(key)
        if rec is None:
            # Allow setting state even if caller forgot to upsert first.
            self.upsert({"key": key})
            rec = self._by_key[key]

        t0 = _now() if t is None else float(t)
        prev = rec.get("last_state")
        if prev != state:
            hist = rec.setdefault("history", [])
            if isinstance(hist, list):
                hist.append({"t": t0, "state": state})
            rec["last_state"] = state
            rec["last_update"] = t0
            self._dirty = True

    def mark_cancel(self, key: str, *, reason: str, t: Optional[float] = None) -> None:
        t0 = _now() if t is None else float(t)
        self.upsert({"key": key, "cancel_reason": reason})
        self.set_state(key, "KILLED", t=t0)
        self.upsert({"key": key, "finished_at": t0})

    def flush(self) -> None:
        if not self._dirty:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        content = json.dumps(self.data, indent=2, sort_keys=True)
        # Atomic write: tmp file + rename to avoid corruption on crash
        fd, tmp = tempfile.mkstemp(
            dir=str(self.path.parent), suffix=".tmp", prefix=".jobs_"
        )
        try:
            os.write(fd, content.encode("utf-8"))
            os.fsync(fd)
            os.close(fd)
            os.replace(tmp, str(self.path))
        except BaseException:
            os.close(fd) if not os.get_inheritable(fd) else None
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise
        self._dirty = False

