#!/usr/bin/env python3
"""
A tiny dummy simulator for local validation of auto_fer_eval scheduling logic.

It mimics the stdout format of the LDPC C++ model:
  - prints [SIM] and [STATISTICS] blocks with RAW BER / TheoRBER / LDPC FER / Total packets simulated.

Usage (matches build_cmd layout):
  dummy_sim.py <sim_mode> <config> <ch_model> <snr> [extra...]
"""

from __future__ import annotations

import math
import sys
import time


def main(argv: list[str]) -> int:
    # build_cmd() passes: <sim_mode> <config> <ch_model> <snr> [extra...]
    # But if exe="python3" and sim_mode="dummy_sim.py", Python consumes the script path,
    # so argv becomes: <config> <ch_model> <snr> [extra...]
    if len(argv) < 3:
        print("usage: dummy_sim.py [<sim_mode>] <config> <ch_model> <snr> [extra...]")
        return 2

    if len(argv) >= 4:
        snr = float(argv[3])
    else:
        snr = float(argv[2])

    # A synthetic "waterfall": log10(FER) is approximately linear in SNR.
    # At snr0, FER≈1; slope decades per dB.
    snr0 = 3.5
    slope = 5.0
    log10_fer = -slope * (snr - snr0)
    fer = 10.0 ** log10_fer
    fer = max(min(fer, 1.0), 1e-12)

    # Couple a fake RAW BER to FER (not physically meaningful; only for testing).
    raw_ber = min(0.5, max(1e-6, fer * 0.05))
    theo_rber = raw_ber * 0.99

    # Simulate runtime: lower FER takes longer (like hitting max packets).
    sleep_s = min(0.8, 0.05 + 0.12 * max(0.0, -math.log10(fer) - 1.0))
    time.sleep(sleep_s)

    total_packets = 1000
    fail_cw = int(round(total_packets * fer))

    print(f"[SIM] Statistical result of {total_packets} packets simulated: ")
    print(f"[SIM] SNR       : {snr:.6f}")
    print(f"[SIM] FAIL CW   : {fail_cw}")
    print(f"[SIM] RAW BER   : {raw_ber:.6e}")
    print(f"[SIM] TheoRBER  : {theo_rber:.6e}")
    print(f"[SIM] LDPC FER  : {fer:.6e}")
    print("[SIM] Finish simulation @ dummy")
    print("--------------------------------------------------------")
    print(f"[STATISTICS] Total packets simulated: {total_packets}")
    print(f"[STATISTICS] RAW BER   : {raw_ber:.6e}")
    print(f"[STATISTICS] TheoRBER  : {theo_rber:.6e}")
    print(f"[STATISTICS] LDPC FER  : {fer:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
