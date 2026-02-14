from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ParsedMetrics:
    raw_ber: Optional[float] = None
    theo_rber: Optional[float] = None
    ldpc_fer: Optional[float] = None
    fail_cw: Optional[int] = None
    total_packets: Optional[int] = None
    retry_dec_avg_iter: Optional[float] = None


_NUM = r"([0-9eE+\-\.]+)"


def parse_log_text(text: str) -> ParsedMetrics:
    """
    Parse metrics from a log string.

    Priority:
      1) [STATISTICS] lines
      2) fallback to the last complete [SIM] block (Statistical result of ...)
    """

    stats = _parse_statistics(text)
    sim = _parse_last_complete_sim_block(text)
    if sim is None:
        return stats

    # Fill missing fields from SIM (stats has higher priority field-by-field).
    return ParsedMetrics(
        raw_ber=stats.raw_ber if stats.raw_ber is not None else sim.raw_ber,
        theo_rber=stats.theo_rber if stats.theo_rber is not None else sim.theo_rber,
        ldpc_fer=stats.ldpc_fer if stats.ldpc_fer is not None else sim.ldpc_fer,
        fail_cw=stats.fail_cw if stats.fail_cw is not None else sim.fail_cw,
        total_packets=stats.total_packets if stats.total_packets is not None else sim.total_packets,
        retry_dec_avg_iter=stats.retry_dec_avg_iter if stats.retry_dec_avg_iter is not None else sim.retry_dec_avg_iter,
    )


def _parse_statistics(text: str) -> ParsedMetrics:
    raw_ber = _last_float(text, rf"^\s*\[STATISTICS\]\s+RAW\s+BER\s*:\s*{_NUM}\s*$")
    theo_rber = _last_float(text, rf"^\s*\[STATISTICS\]\s+TheoRBER\s*:\s*{_NUM}\s*$")
    ldpc_fer = _last_float(text, rf"^\s*\[STATISTICS\]\s+LDPC\s+FER\s*:\s*{_NUM}\s*$")
    fail_cw = _last_int(text, r"^\s*\[STATISTICS\]\s+FAIL\s+CW\s*:\s*([0-9]+)\s*$")
    # Decoder average iterations can have multiple variants depending on the simulator:
    #   - "Retry decoder average iteration" (HREModel)
    #   - "Retry Decoder average iterations" (IBEX)
    #   - "IBEX Decoder average iterations" (IBEX)
    # We parse the last matching "Decoder average iteration(s)" line (case-insensitive).
    retry_it = _last_float(
        text,
        rf"(?i)^\s*\[STATISTICS\]\s+.*\bdecoder\b\s+average\s+iteration(?:s)?\s*:?\s*{_NUM}\s*$",
    )
    total_packets = _last_int(text, r"^\s*\[STATISTICS\]\s+Total packets simulated:\s*([0-9]+)\s*$")
    return ParsedMetrics(
        raw_ber=raw_ber,
        theo_rber=theo_rber,
        ldpc_fer=ldpc_fer,
        fail_cw=fail_cw,
        total_packets=total_packets,
        retry_dec_avg_iter=retry_it,
    )


def _parse_last_complete_sim_block(text: str) -> Optional[ParsedMetrics]:
    # Find all SIM block starts.
    starts = [m.start() for m in re.finditer(r"^\s*\[SIM\]\s+Statistical result of .*?$", text, flags=re.M)]
    if not starts:
        return None

    ends = starts[1:] + [len(text)]

    for s, e in zip(reversed(starts), reversed(ends)):
        blk = text[s:e]
        cand = _parse_sim_block(blk)
        if _is_enough(cand):
            return cand
    return None


def _parse_sim_block(text: str) -> ParsedMetrics:
    total_packets = _last_int(
        text,
        r"^\s*\[SIM\]\s+Statistical result of\s+([0-9]+)\s+packets simulated.*$",
    )
    fail_cw = _last_int(text, r"^\s*\[SIM\]\s+FAIL\s+CW\s*:\s*([0-9]+)\s*$")
    raw_ber = _last_float(text, rf"^\s*\[SIM\]\s+RAW\s+BER\s*:\s*{_NUM}\s*$")
    theo_rber = _last_float(text, rf"^\s*\[SIM\]\s+TheoRBER\s*:\s*{_NUM}\s*$")
    ldpc_fer = _last_float(text, rf"^\s*\[SIM\]\s+LDPC\s+FER\s*:\s*{_NUM}\s*$")
    retry_it = _last_float(
        text,
        rf"(?i)^\s*\[SIM\]\s+.*\bdecoder\b\s+average\s+iteration(?:s)?\s*:?\s*{_NUM}\s*$",
    )
    return ParsedMetrics(
        raw_ber=raw_ber,
        theo_rber=theo_rber,
        ldpc_fer=ldpc_fer,
        fail_cw=fail_cw,
        total_packets=total_packets,
        retry_dec_avg_iter=retry_it,
    )


def _is_enough(m: ParsedMetrics) -> bool:
    # For interval finding, FER is the minimum required signal.
    return m.ldpc_fer is not None


def _last_float(text: str, pattern: str) -> Optional[float]:
    ms = re.findall(pattern, text, flags=re.M)
    if not ms:
        return None
    try:
        return float(ms[-1])
    except ValueError:
        return None


def _last_int(text: str, pattern: str) -> Optional[int]:
    ms = re.findall(pattern, text, flags=re.M)
    if not ms:
        return None
    try:
        return int(ms[-1])
    except ValueError:
        return None
