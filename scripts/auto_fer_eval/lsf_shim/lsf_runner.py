#!/usr/bin/env python3
"""
LSF Shim Runner Process
Background process that executes jobs and updates their state.
"""

import os
import sys
import time
import signal
import subprocess
from pathlib import Path

# Add script directory to Python path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lsf_registry import LsfRegistry


def run_job(job_id: int):
    """
    Execute a job in the background.

    Args:
        job_id: Job ID to execute
    """
    registry = LsfRegistry()
    job = registry.get_job(job_id)

    if job is None:
        print(f"Error: Job {job_id} not found", file=sys.stderr)
        sys.exit(1)

    # Get PEND duration from environment or use default
    pend_sec = int(os.environ.get('LSF_SHIM_PEND_SEC', '10'))

    try:
        # Wait for PEND period
        time.sleep(pend_sec)

        # Transition to RUN state
        registry.set_job_state(job_id, 'RUN')

        # Prepare log file
        log_path = Path(job['log_path'])
        log_path.parent.mkdir(parents=True, exist_ok=True)

        # Open log file for writing
        with open(log_path, 'w') as log_file:
            # Start the process in a new process group
            # This allows bkill to terminate the entire job tree
            process = subprocess.Popen(
                job['cmd'],
                cwd=job['cwd'],
                stdout=log_file,
                stderr=subprocess.STDOUT,
                preexec_fn=os.setpgrp  # Create new process group
            )

            # Record PID and PGID
            pid = process.pid
            pgid = os.getpgid(pid)
            registry.set_job_pid(job_id, pid, pgid)

            # Wait for process to complete
            rc = process.wait()

            # Update final state based on return code
            if rc == 0:
                registry.set_job_state(job_id, 'DONE', rc)
            else:
                registry.set_job_state(job_id, 'EXIT', rc)

    except Exception as e:
        # On any error, mark job as EXIT
        print(f"Error running job {job_id}: {e}", file=sys.stderr)
        registry.set_job_state(job_id, 'EXIT', rc=-1)
        sys.exit(1)


def main():
    """Main entry point for runner."""
    if len(sys.argv) != 2:
        print("Usage: lsf_runner.py <job_id>", file=sys.stderr)
        sys.exit(1)

    try:
        job_id = int(sys.argv[1])
    except ValueError:
        print(f"Error: Invalid job_id '{sys.argv[1]}'", file=sys.stderr)
        sys.exit(1)

    # Detach from parent process
    # This allows bsub to return immediately while runner continues
    if os.fork() > 0:
        # Parent exits immediately
        sys.exit(0)

    # Child continues as daemon
    os.setsid()  # Create new session

    # Second fork to prevent zombie
    if os.fork() > 0:
        sys.exit(0)

    # Redirect stdin/stdout/stderr to /dev/null
    # (job output goes to log file, not runner's stdout)
    devnull = os.open(os.devnull, os.O_RDWR)
    os.dup2(devnull, sys.stdin.fileno())
    os.dup2(devnull, sys.stdout.fileno())
    os.dup2(devnull, sys.stderr.fileno())
    os.close(devnull)

    # Execute the job
    run_job(job_id)


if __name__ == '__main__':
    main()
