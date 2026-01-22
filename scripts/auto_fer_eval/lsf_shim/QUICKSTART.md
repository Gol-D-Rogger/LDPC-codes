# LSF Shim Quick Start

## Installation (One-time setup)

```bash
cd /mnt/f/MyProject/MaxioProject/LDPC-codes/HREModel/auto_fer_eval/lsf_shim
./install.sh

# Add to PATH (add this to ~/.bashrc or ~/.zshrc)
export PATH="/mnt/f/MyProject/MaxioProject/LDPC-codes/HREModel/auto_fer_eval/lsf_shim:$PATH"
```

## Basic Usage

```bash
# Submit a job
bsub -o /tmp/output.log -cwd /tmp echo "Hello World"
# Output: Job <1> is submitted to default queue <normal>.

# Check job status
bjobs 1
# Shows: PEND → RUN → DONE/EXIT

# Cancel a job
bkill 1
```

## Configuration

```bash
# Optional: Customize database location
export LSF_SHIM_DB=/path/to/your/project/lsf_jobs.json

# Optional: Reduce PEND delay (default: 10 seconds)
export LSF_SHIM_PEND_SEC=5
```

## Integration with auto_fer_eval

1. Add shim to PATH (see above)
2. Configure auto_fer_eval to use LSF executor
3. Run your FER evaluations as normal

## Verification

```bash
# Quick test
bsub -o /tmp/test.log -cwd /tmp echo "Test"
sleep 3
bjobs 1  # Should show DONE
cat /tmp/test.log  # Should contain "Test"
```

## Portability

To use in another project:
1. Copy entire `lsf_shim/` directory
2. Run `./install.sh` in new location
3. Add to PATH or set command paths in project config
4. No code changes needed!

## Files Created

- `lsf_registry.py` - Job database management
- `lsf_runner.py` - Background job executor
- `bsub` - Job submission command
- `bjobs` - Job status query command
- `bkill` - Job cancellation command
- `install.sh` - Installation script
- `test_shim.sh` - Test suite
- `README.md` - Full documentation
- `QUICKSTART.md` - This file

## Troubleshooting

**Jobs stuck in PEND?**
- Check that `lsf_runner.py` is executable
- Verify Python 3.6+ is available

**"Job not found" errors?**
- Verify `LSF_SHIM_DB` points to correct file
- Check that job was submitted successfully

**Log files not created?**
- Ensure parent directory exists
- Use absolute paths for `-o` argument

For more details, see README.md
