from __future__ import annotations

import re
import subprocess
import time
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class LocalJob:
    cmd: list[str]
    cwd: str
    log_path: Path
    proc: subprocess.Popen[bytes]
    stdout_fh: object


@dataclass(frozen=True)
class LsfJob:
    cmd: list[str]
    cwd: str
    log_path: Path
    job_id: int


@dataclass(frozen=True)
class DryRunJob:
    """A fake job for dry-run mode."""
    cmd: list[str]
    cwd: str
    log_path: Path
    job_id: int
    queue: str


@dataclass(frozen=True)
class JobState:
    """
    A minimal, backend-agnostic job state abstraction.

    state:
      - "PEND": queued
      - "RUN": running
      - "DONE": finished successfully
      - "EXIT": finished with failure or killed
      - "UNKNOWN": cannot be determined (treat as still running, unless log exists)
    """

    state: str
    done: bool
    ok: bool


class DryRunExecutor:
    """
    A fake executor for dry-run mode.
    Prints what would be submitted without actually running anything.
    """

    def __init__(
        self, 
        *, 
        queue_slow: str = "", 
        queue_fast: str = "", 
        log_base_dir: str = "",
        bsub_extra: Optional[list[str]] = None,
        include_cwd: bool = True,
        ) -> None:
        self.queue_slow = queue_slow
        self.queue_fast = queue_fast
        self.log_base_dir = log_base_dir
        self.bsub_extra = list(bsub_extra or [])
        self.include_cwd = bool(include_cwd)
        self._job_counter = 0
        self._submitted_jobs: list[DryRunJob] = []

    def submit(
        self,
        cmd: list[str],
        *,
        cwd: str,
        log_path: Path,
        job_name: str = "",
        queue: str = "",
    ) -> DryRunJob:
        self._job_counter += 1
        log_path = log_path.resolve()

        # Apply log_base_dir transformation (same logic as LsfExecutor)
        actual_log_path = log_path
        if self.log_base_dir:
            parts = log_path.parts
            for i, part in enumerate(parts):
                if part == "adaptive" and i > 0:
                    rel_parts = parts[i - 1 :]
                    actual_log_path = Path(self.log_base_dir).joinpath(*rel_parts)
                    break

        cwd_abs = str(Path(cwd).resolve())

        bsub: list[str] = ["bsub", "-o", str(actual_log_path)]
        if self.include_cwd:
            bsub.extend(["-cwd", cwd_abs])
        use_queue = queue  # 注意：DryRunExecutor 只接受调用方传入的 queue（通常是 queue_slow/queue_fast）
        if use_queue:
            bsub.extend(["-q", use_queue])
        if job_name:
            bsub.extend(["-J", job_name])
        bsub.extend(self.bsub_extra)
        bsub.extend(cmd)

        cmd_str = " ".join(shlex.quote(s) for s in bsub)
        print(f"[dry-run] {cmd_str}")

        job = DryRunJob(cmd=cmd, cwd=cwd, log_path=actual_log_path, job_id=self._job_counter, queue=queue)
        self._submitted_jobs.append(job)
        return job

    def poll(self, job: DryRunJob) -> JobState:
        # In dry-run mode, jobs are always "done" immediately with a special state
        # that signals to skip log parsing
        return JobState(state="DRY_RUN_DONE", done=True, ok=True)

    def cancel(self, job: DryRunJob) -> None:
        print(f"[dry-run] bkill {job.job_id} (snr point)")

    def resume_job(
        self,
        *,
        cmd: list[str],
        cwd: str,
        log_path: Path,
        backend: str,
        job_id: str,
    ) -> Optional[DryRunJob]:
        return None

    def peek_text(self, job: DryRunJob, *, timeout_sec: float = 5.0) -> Optional[str]:
        return None

    def get_submitted_jobs(self) -> list[DryRunJob]:
        return list(self._submitted_jobs)


class LocalExecutor:
    def submit(
        self,
        cmd: list[str],
        *,
        cwd: str,
        log_path: Path,
        job_name: str = "",
        queue: str = "",  # ignored for local executor
    ) -> LocalJob:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = log_path.open("wb")
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
        return LocalJob(cmd=cmd, cwd=cwd, log_path=log_path, proc=proc, stdout_fh=fh)

    def poll(self, job: LocalJob) -> JobState:
        rc = job.proc.poll()
        if rc is None:
            return JobState(state="RUN", done=False, ok=False)
        try:
            job.stdout_fh.close()
        except Exception:
            pass
        if rc == 0:
            return JobState(state="DONE", done=True, ok=True)
        return JobState(state="EXIT", done=True, ok=False)

    def cancel(self, job: LocalJob, *, timeout_sec: float = 5.0) -> None:
        try:
            job.proc.terminate()
        except Exception:
            return
        t0 = time.time()
        while time.time() - t0 < timeout_sec:
            if job.proc.poll() is not None:
                break
            time.sleep(0.1)
        if job.proc.poll() is None:
            try:
                job.proc.kill()
            except Exception:
                pass
        try:
            job.stdout_fh.close()
        except Exception:
            pass

    def resume_job(
        self,
        *,
        cmd: list[str],
        cwd: str,
        log_path: Path,
        backend: str,
        job_id: str,
    ) -> Optional[LocalJob]:
        return None

    def peek_text(self, job: LocalJob, *, timeout_sec: float = 5.0) -> Optional[str]:
        try:
            if not job.log_path.exists():
                return None
            return job.log_path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            return None


class LsfExecutor:
    def __init__(
        self,
        *,
        bsub_cmd: str = "bsub",
        bjobs_cmd: str = "bjobs",
        bkill_cmd: str = "bkill",
        bpeek_cmd: str = "bpeek",
        queue: str = "",
        queue_slow: str = "",
        queue_fast: str = "",
        log_base_dir: str = "",
        bsub_extra: Optional[list[str]] = None,
        bjobs_extra: Optional[list[str]] = None,
        bkill_extra: Optional[list[str]] = None,
        bpeek_extra: Optional[list[str]] = None,
        include_cwd: bool = True,
    ) -> None:
        self.bsub_cmd = bsub_cmd
        self.bjobs_cmd = bjobs_cmd
        self.bkill_cmd = bkill_cmd
        self.bpeek_cmd = bpeek_cmd
        self.queue = queue
        self.queue_slow = queue_slow or queue  # fallback to queue
        self.queue_fast = queue_fast or queue  # fallback to queue
        self.log_base_dir = log_base_dir
        self.bsub_extra = list(bsub_extra or [])
        self.bjobs_extra = list(bjobs_extra or [])
        self.bkill_extra = list(bkill_extra or [])
        self.bpeek_extra = list(bpeek_extra or [])
        self.include_cwd = bool(include_cwd)
        self._poll_warned: set[int] = set()

    @staticmethod
    def _log_has_final_statistics(log_path: Path) -> bool:
        """
        Return True only when the job log reached the final [STATISTICS] block.

        Partial [SIM] blocks are not sufficient to conclude that the LSF job has
        finished; they only prove the process emitted some intermediate metrics.
        """
        if not log_path.exists():
            return False
        try:
            txt = log_path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            return False
        return re.search(r"^\s*\[STATISTICS\]\s+LDPC\s+FER\s*:\s*[0-9eE+\-\.]+\s*$", txt, flags=re.M) is not None

    def _fallback_state_from_log(self, job: LsfJob) -> JobState:
        """
        Fallback when bjobs output is not available.

        Policy:
          - If log already contains final [STATISTICS]: treat as DONE.
          - Else if log exists and non-empty: treat as RUN (best-effort).
          - Else: treat as UNKNOWN/PEND-like.
        """
        try:
            if job.log_path.exists() and job.log_path.stat().st_size > 0:
                if self._log_has_final_statistics(job.log_path):
                    return JobState(state="DONE", done=True, ok=True)
                return JobState(state="RUN", done=False, ok=False)
        except Exception:
            pass
        return JobState(state="UNKNOWN", done=False, ok=False)

    def resume_job(
        self,
        *,
        cmd: list[str],
        cwd: str,
        log_path: Path,
        backend: str,
        job_id: str,
    ) -> Optional[LsfJob]:
        if str(backend).lower() != "lsf":
            return None
        try:
            job_id_i = int(str(job_id))
        except (TypeError, ValueError):
            return None
        return LsfJob(cmd=list(cmd), cwd=str(cwd), log_path=Path(log_path).resolve(), job_id=job_id_i)

    def submit(
        self,
        cmd: list[str],
        *,
        cwd: str,
        log_path: Path,
        job_name: str = "",
        queue: str = "",
    ) -> LsfJob:
        """
        Submit a job to LSF.

        Args:
            cmd: Command to run
            cwd: Working directory
            log_path: Path to write stdout/stderr (relative to out_dir)
            job_name: Optional job name for -J
            queue: Queue to use. If empty, uses self.queue (default queue).
                   Caller can pass self.queue_slow or self.queue_fast.
        """
        log_path = log_path.resolve()

        # If log_base_dir is set, redirect log output to that directory
        # while preserving the relative structure from log_path
        actual_log_path = log_path
        if self.log_base_dir:
            # Extract the relative part (case_name/log_file.log)
            # log_path is typically: out_dir/case_name/adaptive/main/log_file.log
            # We want: log_base_dir/case_name/adaptive/main/log_file.log
            # The caller should pass the full path, we just replace the base
            actual_log_path = Path(self.log_base_dir) / log_path.name
            # Try to preserve more structure if possible
            # Look for "adaptive" in the path to find the case structure
            parts = log_path.parts
            for i, part in enumerate(parts):
                if part == "adaptive" and i > 0:
                    # Found adaptive, take everything from case_name onwards
                    # parts[i-1] is case_name, parts[i] is "adaptive", etc.
                    rel_parts = parts[i - 1 :]
                    actual_log_path = Path(self.log_base_dir).joinpath(*rel_parts)
                    break

        actual_log_path.parent.mkdir(parents=True, exist_ok=True)
        cwd_abs = str(Path(cwd).resolve())

        bsub: list[str] = [self.bsub_cmd, "-o", str(actual_log_path)]
        if self.include_cwd:
            bsub.extend(["-cwd", cwd_abs])
        use_queue = queue or self.queue
        if use_queue:
            bsub.extend(["-q", use_queue])
        if job_name:
            bsub.extend(["-J", job_name])
        bsub.extend(self.bsub_extra)
        bsub.extend(cmd)

        p = subprocess.run(
            bsub,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            check=False,
        )
        out = p.stdout or ""
        m = re.search(r"Job <([0-9]+)>", out)
        if not m:
            raise RuntimeError(f"bsub failed or unrecognized output:\n{out}")
        job_id = int(m.group(1))
        # Return with actual_log_path so poll/parse can find the log
        return LsfJob(cmd=cmd, cwd=cwd, log_path=actual_log_path, job_id=job_id)

    def poll(self, job: LsfJob) -> JobState:
        # Default: query one job id. If bjobs returns non-zero, LSF may have purged
        # DONE jobs; in that case, rely on the presence of the log as the truth.
        bjobs = [self.bjobs_cmd]
        bjobs.extend(self.bjobs_extra)
        bjobs.append(str(job.job_id))
        try:
            p = subprocess.run(
                bjobs,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                errors="replace",
                check=False,
            )
        except Exception as e:
            # Do not crash the whole scheduler on transient LSF CLI issues.
            # Best-effort infer the state from the log file, and warn once per job id.
            if job.job_id not in self._poll_warned:
                ts = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                print(
                    f"[auto_fer_eval][lsf] WARN: bjobs poll failed for job_id={job.job_id}: "
                    f"{type(e).__name__}: {e} [{ts}]"
                )
                print(f"[auto_fer_eval][lsf] WARN: bjobs cmd: {' '.join(bjobs)} [{ts}]")
                self._poll_warned.add(job.job_id)
            return self._fallback_state_from_log(job)
        out = (p.stdout or "").strip()
        if p.returncode != 0 or not out:
            # LSF may have purged DONE jobs, or bjobs may be temporarily unavailable.
            # Prefer a conservative log-based check.
            if self._log_has_final_statistics(job.log_path):
                return JobState(state="DONE", done=True, ok=True)
            if job.log_path.exists() and job.log_path.stat().st_size > 0:
                return JobState(state="RUN", done=False, ok=False)
            return JobState(state="UNKNOWN", done=False, ok=False)

        lines = [ln.strip() for ln in out.splitlines() if ln.strip()]
        # Typical output:
        #   JOBID USER STAT QUEUE ...
        #   12345 foo  RUN  normal ...
        for ln in reversed(lines):
            if ln.startswith("JOBID"):
                continue
            toks = ln.split()
            if not toks:
                continue
            if toks[0] != str(job.job_id):
                continue
            stat = toks[2] if len(toks) >= 3 else ""
            stat = stat.upper()
            if stat == "DONE":
                return JobState(state="DONE", done=True, ok=True)
            if stat in {"EXIT", "ZOMBI"}:
                return JobState(state="EXIT", done=True, ok=False)
            if stat == "UNKWN":
                # UNKWN is often transient (e.g., exec host / sbatchd unreachable).
                # Do NOT map it to EXIT. Prefer log-based truth if available.
                return self._fallback_state_from_log(job)
            if stat in {"PEND", "PSUSP"}:
                return JobState(state="PEND", done=False, ok=False)
            return JobState(state="RUN", done=False, ok=False)

        return JobState(state="UNKNOWN", done=False, ok=False)

    def poll_batch(self, jobs: list[LsfJob]) -> dict[int, JobState]:
        """Query multiple job IDs in a single bjobs call."""
        if not jobs:
            return {}
        bjobs_cmd = [self.bjobs_cmd] + list(self.bjobs_extra) + [str(j.job_id) for j in jobs]
        job_map = {j.job_id: j for j in jobs}
        result: dict[int, JobState] = {}

        try:
            p = subprocess.run(
                bjobs_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                errors="replace",
                check=False,
            )
        except Exception:
            # Fall back to per-job poll on failure
            for j in jobs:
                result[j.job_id] = self.poll(j)
            return result

        out = (p.stdout or "").strip()
        if p.returncode != 0 or not out:
            for j in jobs:
                result[j.job_id] = self._fallback_state_from_log(j)
            return result

        # Parse bjobs tabular output
        parsed_ids: set[int] = set()
        for ln in out.splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("JOBID"):
                continue
            toks = ln.split()
            if not toks:
                continue
            try:
                jid = int(toks[0])
            except (ValueError, IndexError):
                continue
            if jid not in job_map:
                continue
            stat = toks[2].upper() if len(toks) >= 3 else ""
            if stat == "DONE":
                result[jid] = JobState(state="DONE", done=True, ok=True)
            elif stat in {"EXIT", "ZOMBI"}:
                result[jid] = JobState(state="EXIT", done=True, ok=False)
            elif stat == "UNKWN":
                result[jid] = self._fallback_state_from_log(job_map[jid])
            elif stat in {"PEND", "PSUSP"}:
                result[jid] = JobState(state="PEND", done=False, ok=False)
            else:
                result[jid] = JobState(state="RUN", done=False, ok=False)
            parsed_ids.add(jid)

        # Jobs not found in output (purged) — fall back to log
        for j in jobs:
            if j.job_id not in parsed_ids:
                result[j.job_id] = self._fallback_state_from_log(j)

        return result

    def cancel(self, job: LsfJob) -> None:
        bkill = [self.bkill_cmd]
        bkill.extend(self.bkill_extra)
        bkill.append(str(job.job_id))
        subprocess.run(bkill, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)

    def peek_text(self, job: LsfJob, *, timeout_sec: float = 5.0) -> Optional[str]:
        bpeek = [self.bpeek_cmd]
        bpeek.extend(self.bpeek_extra)
        bpeek.append(str(job.job_id))
        try:
            p = subprocess.run(
                bpeek,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                errors="replace",
                check=False,
                timeout=max(0.1, float(timeout_sec)),
            )
        except Exception:
            return None
        out = p.stdout or ""
        if p.returncode != 0 or not out.strip():
            return None
        return out
