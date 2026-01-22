#!/usr/bin/env python3
"""
LSF Shim Registry Module
Manages job state in a JSON database with file locking for concurrency safety.
"""

import json
import os
import time
import fcntl
from pathlib import Path
from typing import Dict, Optional, Any
from contextlib import contextmanager


class LsfRegistry:
    """Thread-safe job registry using file-based locking."""

    def __init__(self, db_path: Optional[str] = None):
        """
        Initialize registry.

        Args:
            db_path: Path to JSON database. If None, uses LSF_SHIM_DB env var
                     or defaults to ~/.lsf_shim/jobs.json
        """
        if db_path is None:
            db_path = os.environ.get('LSF_SHIM_DB')
            if db_path is None:
                db_path = os.path.expanduser('~/.lsf_shim/jobs.json')

        self.db_path = Path(db_path)
        self.lock_path = Path(str(self.db_path) + '.lock')

        # Ensure directory exists
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        # Initialize empty database if it doesn't exist
        if not self.db_path.exists():
            self._write_db({})

    @contextmanager
    def _lock(self):
        """Context manager for file locking."""
        lock_file = open(self.lock_path, 'w')
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
            lock_file.close()

    def _read_db(self) -> Dict[str, Any]:
        """Read database without locking (caller must hold lock)."""
        if not self.db_path.exists():
            return {}

        with open(self.db_path, 'r') as f:
            return json.load(f)

    def _write_db(self, data: Dict[str, Any]):
        """Write database without locking (caller must hold lock)."""
        with open(self.db_path, 'w') as f:
            json.dump(data, f, indent=2)

    def get_next_job_id(self) -> int:
        """Get next available job ID (thread-safe)."""
        with self._lock():
            db = self._read_db()
            if not db:
                return 1
            return max(int(k) for k in db.keys()) + 1

    def create_job(self, job_id: int, cwd: str, cmd: list, log_path: str,
                   queue: Optional[str] = None, job_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Create a new job entry.

        Args:
            job_id: Unique job ID
            cwd: Working directory
            cmd: Command array
            log_path: Output log path
            queue: Queue name (optional)
            job_name: Job name (optional)

        Returns:
            Created job record
        """
        job = {
            'job_id': job_id,
            'state': 'PEND',
            'rc': None,
            'submit_time': time.time(),
            'start_time': None,
            'end_time': None,
            'cwd': cwd,
            'cmd': cmd,
            'log_path': log_path,
            'queue': queue,
            'job_name': job_name,
            'pid': None,
            'pgid': None
        }

        with self._lock():
            db = self._read_db()
            db[str(job_id)] = job
            self._write_db(db)

        return job

    def get_job(self, job_id: int) -> Optional[Dict[str, Any]]:
        """Get job by ID."""
        with self._lock():
            db = self._read_db()
            return db.get(str(job_id))

    def update_job(self, job_id: int, updates: Dict[str, Any]):
        """
        Update job fields.

        Args:
            job_id: Job ID to update
            updates: Dictionary of fields to update
        """
        with self._lock():
            db = self._read_db()
            job_key = str(job_id)
            if job_key in db:
                db[job_key].update(updates)
                self._write_db(db)

    def set_job_state(self, job_id: int, state: str, rc: Optional[int] = None):
        """
        Update job state.

        Args:
            job_id: Job ID
            state: New state (PEND/RUN/DONE/EXIT)
            rc: Return code (for DONE/EXIT states)
        """
        updates = {'state': state}

        if state == 'RUN' and self.get_job(job_id).get('start_time') is None:
            updates['start_time'] = time.time()

        if state in ('DONE', 'EXIT'):
            updates['end_time'] = time.time()
            if rc is not None:
                updates['rc'] = rc

        self.update_job(job_id, updates)

    def set_job_pid(self, job_id: int, pid: int, pgid: int):
        """Set job process ID and process group ID."""
        self.update_job(job_id, {'pid': pid, 'pgid': pgid})

    def list_jobs(self) -> Dict[str, Any]:
        """List all jobs."""
        with self._lock():
            return self._read_db()

    def delete_job(self, job_id: int):
        """Delete a job from registry."""
        with self._lock():
            db = self._read_db()
            job_key = str(job_id)
            if job_key in db:
                del db[job_key]
                self._write_db(db)
