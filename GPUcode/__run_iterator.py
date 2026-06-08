#!/usr/bin/python3
import os
from typing import Any
import time
import threading

# Constants matching main.cpp defaults / OCR parameter block.
TOTAL_CODEWORDS = 10000000000000
MAX_FAILURE_COUNT = 16
BYTES_USERDATA = 4112
BYTES_PARITY = 472
NAND_STROBES = 7
ITERATION_LIMIT = 32
SAVE_DECODE_FAIL_DATA = False
GPU_MODE = True
FINITE = 1
Q = 9
C = 7
F = 3
LLR = 7
LLR_f = 3

RBER_CONFIG = {
    "start": 0.015,
    "end": 0.023,
    "step": 0.001,
}


def build_project() -> bool:
    return os.system("make clean; make build") == 0


def run_test(rber: float) -> None:
    nand_strobes = 0 if NAND_STROBES == 1 else NAND_STROBES
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    if True:
        output_dir = (
            f"result/{BYTES_USERDATA}_{BYTES_PARITY}_strobe{nand_strobes}_max{ITERATION_LIMIT}"
            f"_F{FINITE}Q{Q}R{C}F{F}LLR{LLR}{LLR_f}_KunyangMatrix"
        )
    elif ERR_INJ_MODE == 1:
        output_dir = (
            f"{BYTES_USERDATA}_{BYTES_PARITY}_strobe{nand_strobes}_max{ITERATION_LIMIT} _post{POST_ITERATION}"
            f"_original_levels"
        )
    else:
        output_dir = (
            f"{BYTES_USERDATA}_ {BYTES_PARITY} strobe{nand_strobes}_max{ITERATION_LIMIT}_post{POST_ITERATION}"
            f" _target_rber_{TARGET_DISTRIBUTION_RBER}_prng_verilog_{PRNG_VERILOG_MODE}"
        )
    os.system(f"mkdir {output_dir}")

    cmd = (
        f"./a.out "
        f"--total_codewords {TOTAL_CODEWORDS} "
        f"--max_failure_count {MAX_FAILURE_COUNT} "
        f"--bytes_of_userdata {BYTES_USERDATA} "
        f"--bytes_of_parity {BYTES_PARITY} "
        f"--nand_strobes {nand_strobes} "
        f"--iteration_limit {ITERATION_LIMIT} "
        f"--rber {rber:.4f} "
        f"--GPU_mode {1 if GPU_MODE else 0} "
        f"--finite_mode {FINITE} "
        f"--finite_q_num {Q} "
        f"--finite_f_num {F} "
        f"--finite_c_num {C} "
        f"--llr_tot_bit {LLR} "
        f"--llr_frac_bit {LLR_f} "
    )

    rber_str = f"{rber:.4f}".replace(".", "p")
    print(f"Running Command:\n{cmd}")
    os.system(f"{cmd} > {output_dir}/rber_{rber_str}.txt")

if __name__ == "__main__":
    if not build_project():
        print("Make rebuild failed. Exiting.")
        exit(1)
    else:
        # build rber list from config and iterate (reversed order)
        start = RBER_CONFIG["start"]
        end = RBER_CONFIG["end"] + 0.000001
        step = RBER_CONFIG["step"]

        # Generate rber list without numpy
        num_steps = int((end - start) / step) + 1
        rbers = [start + i * step for i in range(num_steps)]

        for rber in reversed(rbers):
            run_test(rber)
