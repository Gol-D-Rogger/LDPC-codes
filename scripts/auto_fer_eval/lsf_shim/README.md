# LSF Shim - Lightweight LSF Simulator

A portable, lightweight simulator for IBM LSF (Load Sharing Facility) commands (`bsub`, `bjobs`, `bkill`) designed for local development and testing without requiring a real LSF cluster.

## Features

- **Minimal dependencies**: Python 3.6+ standard library only
- **Portable**: Works on Linux, macOS, and WSL
- **Self-contained**: Single directory, easy to copy to other projects
- **LSF-compatible output**: Designed to work with tools that parse LSF command output
- **Concurrent-safe**: File-based locking for multi-process job submissions
- **Configurable**: Environment variables for customization

## Quick Start

### Installation

```bash
cd auto_fer_eval/lsf_shim
./install.sh
```

Add to your PATH (add to `~/.bashrc` or `~/.zshrc`):

```bash
export PATH="/path/to/auto_fer_eval/lsf_shim:$PATH"
```

### Basic Usage

```bash
# Submit a job
bsub -o /tmp/output.log -cwd /tmp echo "Hello World"
# Output: Job <1> is submitted to default queue <normal>.

# Check job status
bjobs 1
# Output:
# JOBID   USER    STAT    QUEUE   FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
# 1       user    RUN     normal  localhost   localhost   echo Hello...  Jan 21 14:00

# Cancel a job
bkill 1
# Output: Job <1> is being terminated
```

## Commands

### bsub - Submit Jobs

Submit a job to the LSF shim scheduler.

**Syntax:**
```bash
bsub -o <log_path> [-cwd <working_dir>] [-q <queue>] [-J <job_name>] [other_options...] <command...>
```

**Required Arguments:**
- `-o <log_path>`: Output log file path (stdout/stderr will be redirected here)
- `<command...>`: Command to execute

**Optional Arguments:**
- `-cwd <working_dir>`: Working directory (default: current directory)
- `-q <queue>`: Queue name (recorded but not used for scheduling)
- `-J <job_name>`: Job name (for display purposes)
- Other LSF flags are accepted but ignored

**Output Format:**
```
Job <123> is submitted to default queue <normal>.
```

**Example:**
```bash
bsub -o /tmp/test.log -cwd /tmp -J my_test echo "Test job"
```

### bjobs - Query Job Status

Query the status of a submitted job.

**Syntax:**
```bash
bjobs [options...] <job_id>
```

**Output Format:**
```
JOBID   USER    STAT    QUEUE   FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
123     user    RUN     normal  localhost   localhost   echo Test...  Jan 21 14:00
```

**Job States:**
- `PEND`: Job is pending (waiting for PEND delay, default 10 seconds)
- `RUN`: Job is currently running
- `DONE`: Job completed successfully (exit code 0)
- `EXIT`: Job failed (non-zero exit code) or was killed

**Example:**
```bash
bjobs 123
```

### bkill - Cancel Jobs

Cancel a running or pending job.

**Syntax:**
```bash
bkill [options...] <job_id>
```

**Behavior:**
- Sends SIGTERM to the job's process group
- Waits 2 seconds for graceful termination
- Sends SIGKILL if process hasn't terminated
- Updates job state to EXIT

**Output Format:**
```
Job <123> is being terminated
```

**Example:**
```bash
bkill 123
```

## Configuration

### Environment Variables

#### LSF_SHIM_DB
Path to the job database file.

**Default:** `~/.lsf_shim/jobs.json`

**Example:**
```bash
export LSF_SHIM_DB=/tmp/my_project/lsf_jobs.json
```

#### LSF_SHIM_PEND_SEC
Duration (in seconds) that jobs remain in PEND state before transitioning to RUN.

**Default:** `10`

**Example:**
```bash
export LSF_SHIM_PEND_SEC=5  # Jobs start running after 5 seconds
```

## Architecture

### Components

1. **lsf_registry.py**: Job state management
   - JSON-based persistent storage
   - File locking for concurrent access
   - CRUD operations for job records

2. **lsf_runner.py**: Background job executor
   - Double-fork daemon pattern
   - Process group management
   - Log file redirection
   - State transitions (PEND → RUN → DONE/EXIT)

3. **bsub**: Job submission command
   - Parses arguments flexibly
   - Creates job record
   - Spawns runner process
   - Returns immediately with job ID

4. **bjobs**: Job status query
   - Reads job from registry
   - Validates process state
   - Formats output for LSF compatibility

5. **bkill**: Job cancellation
   - Terminates process group
   - Updates job state to EXIT

### Data Flow

```
┌─────────┐
│  bsub   │ Create job (PEND) → Spawn runner → Return job_id
└─────────┘
     │
     ↓
┌─────────────┐
│ lsf_runner  │ Wait PEND_SEC → RUN → Execute cmd → DONE/EXIT
└─────────────┘
     │
     ↓
┌─────────────┐
│  Registry   │ JSON database with file locking
└─────────────┘
     ↑
     │
┌─────────┬─────────┐
│ bjobs   │ bkill   │ Read/update job state
└─────────┴─────────┘
```

### Job State Machine

```
PEND (10s default)
  │
  ↓
RUN ──────────┬──→ DONE (exit code 0)
  │           │
  │           └──→ EXIT (exit code ≠ 0)
  │
  └──→ EXIT (bkill)
```

## Integration with auto_fer_eval

The LSF shim is designed to work seamlessly with `auto_fer_eval`'s `LsfExecutor`:

1. Add shim to PATH:
   ```bash
   export PATH="/path/to/auto_fer_eval/lsf_shim:$PATH"
   ```

2. Configure `auto_fer_eval` to use LSF executor in your TOML config:
   ```toml
   executor = "lsf"
   ```

3. The shim will handle:
   - Job submission with proper output format
   - Status polling with correct column layout
   - Log file creation at specified paths
   - Job cancellation when needed

## Portability & Reuse

This LSF shim is designed to be easily portable to other projects:

### Copying to Another Project

1. Copy the entire `lsf_shim/` directory to your project
2. Run `./install.sh` in the new location
3. Add to PATH or update command paths in your project config
4. Optionally set `LSF_SHIM_DB` to a project-specific location

### No Code Changes Required

The shim is self-contained and requires no modifications for different projects. All customization is done through:
- Environment variables (`LSF_SHIM_DB`, `LSF_SHIM_PEND_SEC`)
- PATH configuration
- Command-line arguments

### Minimal Dependencies

- Python 3.6+ (standard library only)
- POSIX-compatible OS (Linux, macOS, WSL)
- No external packages required

## Troubleshooting

### Jobs stuck in PEND state

**Cause:** Runner process may have failed to start.

**Solution:**
- Check that `lsf_runner.py` is executable: `chmod +x lsf_runner.py`
- Verify Python 3 is available: `python3 --version`
- Check system logs for errors

### Jobs show RUN but are actually finished

**Cause:** Process died without updating registry.

**Solution:**
- Run `bjobs <job_id>` - it will detect dead processes and update state to EXIT
- Check log file for error messages

### "Job not found" errors

**Cause:** Job ID doesn't exist in registry.

**Solution:**
- Verify job was submitted successfully (check `bsub` output)
- Check that `LSF_SHIM_DB` points to the correct database file
- Ensure database file hasn't been deleted

### Concurrent submission issues

**Cause:** File locking may fail on some filesystems (e.g., NFS).

**Solution:**
- Use a local filesystem for `LSF_SHIM_DB`
- Example: `export LSF_SHIM_DB=/tmp/lsf_jobs.json`

### Log files not created

**Cause:** Parent directory doesn't exist or insufficient permissions.

**Solution:**
- Ensure parent directory of `-o` path exists
- Check write permissions on log directory
- Use absolute paths for `-o` argument

## Testing

A test script is provided to verify the installation:

```bash
./test_shim.sh
```

This will test:
- Job submission and ID parsing
- State transitions (PEND → RUN → DONE)
- Log file creation
- Job cancellation
- Concurrent submissions

## Limitations

This is a **minimal simulator** designed for local testing, not a full LSF replacement:

- **No real scheduling**: Jobs run immediately after PEND delay
- **No resource management**: `-R`, `-W`, `-n` flags are ignored
- **No queue policies**: All jobs run in the same "queue"
- **No job arrays**: `-J "name[1-10]"` syntax not supported
- **No dependencies**: `-w "done(123)"` syntax not supported
- **Limited bjobs output**: Only essential columns are included

For production workloads, use a real LSF cluster.

## License

This LSF shim is provided as-is for development and testing purposes.

## Contributing

To adapt this shim for your project:
1. Copy the directory
2. Modify as needed (all code is in Python, easy to extend)
3. Submit improvements back if they benefit the general use case

## Version

Version: 1.0.0
Last updated: 2026-01-21
