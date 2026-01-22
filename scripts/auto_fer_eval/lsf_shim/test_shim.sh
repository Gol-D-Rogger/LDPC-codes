#!/bin/bash
# LSF Shim Test Script
# Tests basic functionality of bsub, bjobs, and bkill

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="/tmp/lsf_shim_test_$$"
TEST_DB="$TEST_DIR/jobs.json"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counter
TESTS_PASSED=0
TESTS_FAILED=0

# Setup
setup() {
    echo "=== LSF Shim Test Suite ==="
    echo ""
    echo "Setting up test environment..."

    # Create test directory
    mkdir -p "$TEST_DIR"

    # Set environment variables
    export LSF_SHIM_DB="$TEST_DB"
    export LSF_SHIM_PEND_SEC=2  # Shorter PEND for faster tests
    export PATH="$SCRIPT_DIR:$PATH"

    echo "✓ Test directory: $TEST_DIR"
    echo "✓ Test database: $TEST_DB"
    echo "✓ PEND duration: 2 seconds"
    echo ""
}

# Cleanup
cleanup() {
    echo ""
    echo "Cleaning up..."
    rm -rf "$TEST_DIR"
    echo "✓ Test directory removed"
}

# Test result helpers
pass() {
    echo -e "${GREEN}✓ PASS${NC}: $1"
    ((TESTS_PASSED++))
}

fail() {
    echo -e "${RED}✗ FAIL${NC}: $1"
    ((TESTS_FAILED++))
}

info() {
    echo -e "${YELLOW}ℹ${NC} $1"
}

# Test 1: bsub returns job ID
test_bsub_output() {
    echo "Test 1: bsub returns job ID in correct format"

    OUTPUT=$(bsub -o "$TEST_DIR/test1.log" -cwd "$TEST_DIR" echo "Test 1")

    if echo "$OUTPUT" | grep -q "Job <[0-9]\+>"; then
        JOB_ID=$(echo "$OUTPUT" | grep -oP 'Job <\K[0-9]+')
        pass "bsub returned job ID: $JOB_ID"
        echo "$JOB_ID" > "$TEST_DIR/job1_id"
    else
        fail "bsub output format incorrect: $OUTPUT"
    fi
    echo ""
}

# Test 2: bjobs shows PEND state initially
test_bjobs_pend() {
    echo "Test 2: bjobs shows PEND state initially"

    OUTPUT=$(bsub -o "$TEST_DIR/test2.log" -cwd "$TEST_DIR" sleep 10)
    JOB_ID=$(echo "$OUTPUT" | grep -oP 'Job <\K[0-9]+')
    echo "$JOB_ID" > "$TEST_DIR/job2_id"

    sleep 0.5  # Give runner time to start

    BJOBS_OUTPUT=$(bjobs "$JOB_ID")

    if echo "$BJOBS_OUTPUT" | grep -q "PEND"; then
        pass "Job $JOB_ID is in PEND state"
    else
        fail "Job $JOB_ID not in PEND state: $BJOBS_OUTPUT"
    fi
    echo ""
}

# Test 3: Job transitions to RUN state
test_bjobs_run() {
    echo "Test 3: Job transitions to RUN state"

    JOB_ID=$(cat "$TEST_DIR/job2_id")

    info "Waiting for PEND period (2 seconds)..."
    sleep 2.5

    BJOBS_OUTPUT=$(bjobs "$JOB_ID")

    if echo "$BJOBS_OUTPUT" | grep -q "RUN"; then
        pass "Job $JOB_ID transitioned to RUN state"
    else
        fail "Job $JOB_ID not in RUN state: $BJOBS_OUTPUT"
    fi
    echo ""
}

# Test 4: bkill terminates job
test_bkill() {
    echo "Test 4: bkill terminates running job"

    JOB_ID=$(cat "$TEST_DIR/job2_id")

    BKILL_OUTPUT=$(bkill "$JOB_ID")

    if echo "$BKILL_OUTPUT" | grep -q "being terminated"; then
        pass "bkill accepted job $JOB_ID"
    else
        fail "bkill output unexpected: $BKILL_OUTPUT"
    fi

    sleep 0.5

    BJOBS_OUTPUT=$(bjobs "$JOB_ID")

    if echo "$BJOBS_OUTPUT" | grep -q "EXIT"; then
        pass "Job $JOB_ID marked as EXIT after bkill"
    else
        fail "Job $JOB_ID not in EXIT state: $BJOBS_OUTPUT"
    fi
    echo ""
}

# Test 5: Job completes successfully (DONE state)
test_job_done() {
    echo "Test 5: Job completes successfully with DONE state"

    OUTPUT=$(bsub -o "$TEST_DIR/test5.log" -cwd "$TEST_DIR" echo "Success")
    JOB_ID=$(echo "$OUTPUT" | grep -oP 'Job <\K[0-9]+')

    info "Waiting for job to complete (PEND + execution)..."
    sleep 3

    BJOBS_OUTPUT=$(bjobs "$JOB_ID")

    if echo "$BJOBS_OUTPUT" | grep -q "DONE"; then
        pass "Job $JOB_ID completed with DONE state"
    else
        fail "Job $JOB_ID not in DONE state: $BJOBS_OUTPUT"
    fi

    # Check log file
    if [ -f "$TEST_DIR/test5.log" ] && grep -q "Success" "$TEST_DIR/test5.log"; then
        pass "Log file created with correct content"
    else
        fail "Log file missing or incorrect"
    fi
    echo ""
}

# Test 6: Job fails with EXIT state
test_job_exit() {
    echo "Test 6: Job fails with EXIT state (non-zero exit code)"

    OUTPUT=$(bsub -o "$TEST_DIR/test6.log" -cwd "$TEST_DIR" sh -c "exit 1")
    JOB_ID=$(echo "$OUTPUT" | grep -oP 'Job <\K[0-9]+')

    info "Waiting for job to complete..."
    sleep 3

    BJOBS_OUTPUT=$(bjobs "$JOB_ID")

    if echo "$BJOBS_OUTPUT" | grep -q "EXIT"; then
        pass "Job $JOB_ID failed with EXIT state"
    else
        fail "Job $JOB_ID not in EXIT state: $BJOBS_OUTPUT"
    fi
    echo ""
}

# Test 7: Unknown flags don't crash bsub
test_unknown_flags() {
    echo "Test 7: Unknown flags are handled gracefully"

    OUTPUT=$(bsub -o "$TEST_DIR/test7.log" -cwd "$TEST_DIR" -unknown_flag value -another_flag echo "Test")

    if echo "$OUTPUT" | grep -q "Job <[0-9]\+>"; then
        JOB_ID=$(echo "$OUTPUT" | grep -oP 'Job <\K[0-9]+')
        pass "bsub handled unknown flags, job ID: $JOB_ID"
    else
        fail "bsub failed with unknown flags"
    fi
    echo ""
}

# Test 8: Concurrent job submissions
test_concurrent_submissions() {
    echo "Test 8: Concurrent job submissions"

    info "Submitting 5 jobs concurrently..."

    for i in {1..5}; do
        bsub -o "$TEST_DIR/concurrent_$i.log" -cwd "$TEST_DIR" echo "Job $i" &
    done

    wait

    # Check that all jobs were created
    sleep 1

    if [ -f "$TEST_DB" ]; then
        JOB_COUNT=$(python3 -c "import json; print(len(json.load(open('$TEST_DB'))))")
        if [ "$JOB_COUNT" -ge 5 ]; then
            pass "All concurrent jobs were registered (total jobs: $JOB_COUNT)"
        else
            fail "Some concurrent jobs were lost (total jobs: $JOB_COUNT)"
        fi
    else
        fail "Database file not created"
    fi
    echo ""
}

# Test 9: Job with custom working directory
test_custom_cwd() {
    echo "Test 9: Job with custom working directory"

    mkdir -p "$TEST_DIR/custom_dir"
    echo "test content" > "$TEST_DIR/custom_dir/test_file.txt"

    OUTPUT=$(bsub -o "$TEST_DIR/test9.log" -cwd "$TEST_DIR/custom_dir" cat test_file.txt)
    JOB_ID=$(echo "$OUTPUT" | grep -oP 'Job <\K[0-9]+')

    info "Waiting for job to complete..."
    sleep 3

    if [ -f "$TEST_DIR/test9.log" ] && grep -q "test content" "$TEST_DIR/test9.log"; then
        pass "Job executed in custom working directory"
    else
        fail "Job did not execute in custom working directory"
    fi
    echo ""
}

# Test 10: bjobs column format
test_bjobs_format() {
    echo "Test 10: bjobs output format (column 3 is STAT)"

    OUTPUT=$(bsub -o "$TEST_DIR/test10.log" -cwd "$TEST_DIR" echo "Format test")
    JOB_ID=$(echo "$OUTPUT" | grep -oP 'Job <\K[0-9]+')

    sleep 0.5

    BJOBS_OUTPUT=$(bjobs "$JOB_ID")

    # Extract third column from second line (first line is header)
    STAT_COL=$(echo "$BJOBS_OUTPUT" | tail -n 1 | awk '{print $3}')

    if [[ "$STAT_COL" =~ ^(PEND|RUN|DONE|EXIT)$ ]]; then
        pass "bjobs column 3 contains valid STAT: $STAT_COL"
    else
        fail "bjobs column 3 is not STAT: $STAT_COL"
    fi
    echo ""
}

# Main test execution
main() {
    setup

    test_bsub_output
    test_bjobs_pend
    test_bjobs_run
    test_bkill
    test_job_done
    test_job_exit
    test_unknown_flags
    test_concurrent_submissions
    test_custom_cwd
    test_bjobs_format

    cleanup

    # Summary
    echo "=== Test Summary ==="
    echo -e "${GREEN}Passed: $TESTS_PASSED${NC}"
    echo -e "${RED}Failed: $TESTS_FAILED${NC}"
    echo ""

    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "${GREEN}All tests passed!${NC}"
        exit 0
    else
        echo -e "${RED}Some tests failed.${NC}"
        exit 1
    fi
}

# Run tests
main
