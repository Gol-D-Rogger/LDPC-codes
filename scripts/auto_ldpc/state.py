"""
SQLite-backed state tracking for LDPC automation.
"""

from __future__ import annotations

import json
import sqlite3
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


CREATE_TABLES_SQL = """
PRAGMA journal_mode = WAL;
PRAGMA synchronous = NORMAL;

CREATE TABLE IF NOT EXISTS matrices (
    matrix_id TEXT PRIMARY KEY,
    created_at TEXT DEFAULT (datetime('now')),
    updated_at TEXT DEFAULT (datetime('now')),
    anchor_fer REAL,
    anchor_fer_updated_at TEXT,
    err_inj_status TEXT DEFAULT 'pending',
    err_inj_job_id TEXT,
    err_inj_fer REAL,
    err_inj_updated_at TEXT,
    passes_threshold INTEGER DEFAULT 0,
    deep_runs_submitted INTEGER DEFAULT 0,
    total_runs_submitted INTEGER DEFAULT 0,
    notes TEXT
);

CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    matrix_id TEXT NOT NULL,
    model TEXT NOT NULL,
    snr REAL NOT NULL,
    fer REAL,
    metrics TEXT,
    log_path TEXT,
    status TEXT DEFAULT 'completed',
    updated_at TEXT DEFAULT (datetime('now')),
    UNIQUE(matrix_id, model, snr)
);

CREATE TABLE IF NOT EXISTS errors (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    context TEXT,
    message TEXT,
    created_at TEXT DEFAULT (datetime('now'))
);
"""


@dataclass
class MatrixRecord:
    matrix_id: str
    anchor_fer: Optional[float]
    err_inj_status: str
    err_inj_fer: Optional[float]
    passes_threshold: bool
    deep_runs_submitted: bool


class StateStore:
    def __init__(self, db_path: Path) -> None:
        self.db_path = db_path
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(str(db_path))
        self._conn.row_factory = sqlite3.Row
        with self._conn:
            self._conn.executescript(CREATE_TABLES_SQL)

    def close(self) -> None:
        with self._lock:
            self._conn.close()

    def _execute(self, sql: str, params: Sequence[Any] = ()) -> sqlite3.Cursor:
        with self._lock:
            cursor = self._conn.execute(sql, params)
            self._conn.commit()
            return cursor

    def record_error(self, context: str, message: str) -> None:
        self._execute(
            "INSERT INTO errors (context, message) VALUES (?, ?)",
            (context, message),
        )

    def upsert_matrix(self, matrix_id: str) -> None:
        self._execute(
            """
            INSERT INTO matrices (matrix_id) VALUES (?)
            ON CONFLICT(matrix_id) DO UPDATE SET updated_at = datetime('now')
            """,
            (matrix_id,),
        )

    def record_awgn_run(
        self,
        matrix_id: str,
        fer: float,
        metrics: Dict[str, Any],
        snr: float,
        log_path: Optional[Path] = None,
        is_anchor: bool = False,
    ) -> None:
        metrics_json = json.dumps(metrics, ensure_ascii=False)
        self._execute(
            """
            INSERT INTO runs (matrix_id, model, snr, fer, metrics, log_path)
            VALUES (?, 'AWGN', ?, ?, ?, ?)
            ON CONFLICT(matrix_id, model, snr) DO UPDATE SET
                fer = excluded.fer,
                metrics = excluded.metrics,
                log_path = excluded.log_path,
                updated_at = datetime('now')
            """,
            (matrix_id, snr, fer, metrics_json, str(log_path) if log_path else None),
        )
        if is_anchor:
            self._execute(
                """
                INSERT INTO matrices (matrix_id, anchor_fer, anchor_fer_updated_at)
                VALUES (?, ?, datetime('now'))
                ON CONFLICT(matrix_id) DO UPDATE SET
                    anchor_fer = excluded.anchor_fer,
                    anchor_fer_updated_at = excluded.anchor_fer_updated_at,
                    updated_at = datetime('now')
                """,
                (matrix_id, fer),
            )

    def update_err_inj_result(
        self,
        matrix_id: str,
        fer: float,
        metrics: Dict[str, Any],
        log_path: Optional[Path] = None,
        passes_threshold: bool = False,
    ) -> None:
        metrics_json = json.dumps(metrics, ensure_ascii=False)
        self._execute(
            """
            INSERT INTO runs (matrix_id, model, snr, fer, metrics, log_path)
            VALUES (?, 'ERR_INJ', 0, ?, ?, ?)
            ON CONFLICT(matrix_id, model, snr) DO UPDATE SET
                fer = excluded.fer,
                metrics = excluded.metrics,
                log_path = excluded.log_path,
                updated_at = datetime('now')
            """,
            (matrix_id, fer, metrics_json, str(log_path) if log_path else None),
        )
        self._execute(
            """
            INSERT INTO matrices (matrix_id, err_inj_status, err_inj_fer,
                                  err_inj_job_id, err_inj_updated_at, passes_threshold)
            VALUES (?, 'completed', ?, NULL, datetime('now'), ?)
            ON CONFLICT(matrix_id) DO UPDATE SET
                err_inj_status = 'completed',
                err_inj_fer = excluded.err_inj_fer,
                err_inj_job_id = NULL,
                err_inj_updated_at = excluded.err_inj_updated_at,
                passes_threshold = excluded.passes_threshold,
                updated_at = datetime('now')
            """,
            (matrix_id, fer, 1 if passes_threshold else 0),
        )

    def mark_err_inj_submitted(self, matrix_id: str, job_id: str) -> None:
        self._execute(
            """
            INSERT INTO matrices (matrix_id, err_inj_status, err_inj_job_id)
            VALUES (?, 'queued', ?)
            ON CONFLICT(matrix_id) DO UPDATE SET
                err_inj_status = 'queued',
                err_inj_job_id = excluded.err_inj_job_id,
                updated_at = datetime('now')
            """,
            (matrix_id, job_id),
        )

    def mark_err_inj_running(self, matrix_id: str) -> None:
        self._execute(
            """
            UPDATE matrices
            SET err_inj_status = 'running', updated_at = datetime('now')
            WHERE matrix_id = ?
            """,
            (matrix_id,),
        )

    def set_passes_threshold(self, matrix_id: str, passes: bool) -> None:
        self._execute(
            """
            UPDATE matrices
            SET passes_threshold = ?, updated_at = datetime('now')
            WHERE matrix_id = ?
            """,
            (1 if passes else 0, matrix_id),
        )

    def mark_deep_runs_submitted(self, matrix_id: str) -> None:
        self._execute(
            """
            UPDATE matrices
            SET deep_runs_submitted = 1, updated_at = datetime('now')
            WHERE matrix_id = ?
            """,
            (matrix_id,),
        )

    def increment_total_runs(self, matrix_id: str, count: int = 1) -> None:
        self._execute(
            """
            UPDATE matrices
            SET total_runs_submitted = total_runs_submitted + ?
            WHERE matrix_id = ?
            """,
            (count, matrix_id),
        )

    def count_matrices(self) -> int:
        cursor = self._execute("SELECT COUNT(*) FROM matrices")
        return int(cursor.fetchone()[0])

    def fetch_top_matrices(self, top_n: int) -> List[MatrixRecord]:
        cursor = self._execute(
            """
            SELECT matrix_id, anchor_fer, err_inj_status, err_inj_fer,
                   passes_threshold, deep_runs_submitted
            FROM matrices
            WHERE anchor_fer IS NOT NULL
            ORDER BY anchor_fer ASC
            LIMIT ?
            """,
            (top_n,),
        )
        records: List[MatrixRecord] = []
        for row in cursor.fetchall():
            records.append(
                MatrixRecord(
                    matrix_id=row["matrix_id"],
                    anchor_fer=row["anchor_fer"],
                    err_inj_status=row["err_inj_status"],
                    err_inj_fer=row["err_inj_fer"],
                    passes_threshold=bool(row["passes_threshold"]),
                    deep_runs_submitted=bool(row["deep_runs_submitted"]),
                )
            )
        return records

    def matrices_pending_err_inj(self) -> List[str]:
        cursor = self._execute(
            """
            SELECT matrix_id
            FROM matrices
            WHERE err_inj_status IN ('pending', 'failed')
            """,
        )
        return [row["matrix_id"] for row in cursor.fetchall()]

    def matrices_require_deep_runs(self) -> List[str]:
        cursor = self._execute(
            """
            SELECT matrix_id
            FROM matrices
            WHERE passes_threshold = 1 AND deep_runs_submitted = 0
            """,
        )
        return [row["matrix_id"] for row in cursor.fetchall()]

    def update_err_inj_status(self, matrix_id: str, status: str) -> None:
        self._execute(
            """
            UPDATE matrices
            SET err_inj_status = ?, updated_at = datetime('now')
            WHERE matrix_id = ?
            """,
            (status, matrix_id),
        )

    def ensure_batch_matrices(self, matrix_ids: Iterable[str]) -> None:
        with self._lock:
            cursor = self._conn.cursor()
            cursor.executemany(
                """
                INSERT INTO matrices (matrix_id)
                VALUES (?)
                ON CONFLICT(matrix_id) DO NOTHING
                """,
                ((matrix_id,) for matrix_id in matrix_ids),
            )
            self._conn.commit()

    def fetch_runs_for_xlsx(self) -> List[Dict[str, Any]]:
        cursor = self._execute(
            """
            SELECT matrix_id, model, snr, fer, metrics, log_path, updated_at
            FROM runs
            ORDER BY matrix_id ASC, model ASC, snr ASC
            """
        )
        rows = []
        for row in cursor.fetchall():
            metrics = {}
            if row["metrics"]:
                try:
                    metrics = json.loads(row["metrics"])
                except json.JSONDecodeError:
                    metrics = {}
            rows.append(
                {
                    "matrix_id": row["matrix_id"],
                    "model": row["model"],
                    "snr": row["snr"],
                    "fer": row["fer"],
                    "metrics": metrics,
                    "log_path": row["log_path"],
                    "updated_at": row["updated_at"],
                }
            )
        return rows

    def purge_old_jobs(self, retention_hours: Optional[int] = None) -> None:
        if retention_hours is None:
            return
        cutoff = time.time() - retention_hours * 3600
        cutoff_iso = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(cutoff))
        self._execute(
            """
            DELETE FROM runs
            WHERE updated_at < ?
            """,
            (cutoff_iso,),
        )

    def fetch_awgn_runs_for_matrices(self, matrix_ids: Sequence[str]) -> List[Dict[str, Any]]:
        if not matrix_ids:
            return []
        placeholders = ",".join("?" for _ in matrix_ids)
        sql = f"""
            SELECT matrix_id, snr, fer, updated_at
            FROM runs
            WHERE model = 'AWGN' AND matrix_id IN ({placeholders})
            ORDER BY matrix_id ASC, snr ASC
        """
        cursor = self._execute(sql, tuple(matrix_ids))
        rows: List[Dict[str, Any]] = []
        for row in cursor.fetchall():
            rows.append(
                {
                    "matrix_id": row["matrix_id"],
                    "snr": row["snr"],
                    "fer": row["fer"],
                    "updated_at": row["updated_at"],
                }
            )
        return rows
