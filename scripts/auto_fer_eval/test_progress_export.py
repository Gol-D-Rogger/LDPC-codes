from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import auto_fer_eval as afe
from executors import JobState
from planner import StepPolicy


_PARTIAL_LOG = """\
[SIM] Statistical result of 3 packets simulated:
[SIM] FAIL CW   : 3
[SIM] RAW BER   : 1.000000e-01
[SIM] LDPC FER  : 1.000000e+00
"""


_FINAL_LOG = """\
[STATISTICS] Total packets simulated: 10
[STATISTICS] FAIL CW   : 10
[STATISTICS] RAW BER   : 1.000000e-01
[STATISTICS] LDPC FER  : 1.000000e+00
"""


class _FakeJob:
    def __init__(self, *, cmd: list[str], cwd: str, log_path: Path, job_id: int) -> None:
        self.cmd = list(cmd)
        self.cwd = str(cwd)
        self.log_path = Path(log_path).resolve()
        self.job_id = int(job_id)
        self.poll_count = 0


class _FakeExecutor:
    def __init__(self) -> None:
        self._job_id = 0

    def submit(
        self,
        cmd: list[str],
        *,
        cwd: str,
        log_path: Path,
        job_name: str = "",
        queue: str = "",
    ) -> _FakeJob:
        del job_name, queue
        self._job_id += 1
        log_path = Path(log_path).resolve()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("", encoding="utf-8")
        return _FakeJob(cmd=cmd, cwd=cwd, log_path=log_path, job_id=self._job_id)

    def poll(self, job: _FakeJob) -> JobState:
        job.poll_count += 1
        if job.poll_count == 1:
            return JobState(state="RUN", done=False, ok=False)
        job.log_path.write_text(_FINAL_LOG, encoding="utf-8")
        return JobState(state="DONE", done=True, ok=True)

    def resume_job(
        self,
        *,
        cmd: list[str],
        cwd: str,
        log_path: Path,
        backend: str,
        job_id: str,
    ) -> None:
        del cmd, cwd, log_path, backend, job_id
        return None

    def peek_text(self, job: _FakeJob, *, timeout_sec: float = 5.0) -> str | None:
        del timeout_sec
        if job.poll_count < 2:
            return _PARTIAL_LOG
        return _FINAL_LOG


class ProgressExportTests(unittest.TestCase):
    def _make_case(self, root: Path) -> afe.CaseConfig:
        out_dir = root / "out"
        cfg_path = out_dir / "case" / "config_main.cnfg"
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        cfg_path.write_text(
            "\n".join(
                [
                    "0",
                    "0",
                    "0",
                    "0",
                    "0",
                    "0",
                    "0",
                    "0",
                    "0",
                    "10 // maximum error number",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return afe.CaseConfig(
            name="case/adaptive/main",
            workdir=str(root),
            exe="/bin/echo",
            sim_mode="LDPC",
            config=str(cfg_path),
            ch_model="AWGN",
            cmd_extra_args=[],
            out_dir=str(out_dir),
            log_prefix="fer_case",
            axis_type="snr",
            x_start=1.0,
            x_stop=2.0,
            direction=1,
            pilot_step=0.1,
            low_margin=0.0,
            snr_span_cap=1.0,
            fer_hi=1.0,
            fer_lo=1e-6,
            fer_pilot_stop=0.9,
            fer_pilot_too_low=1e-3,
            max_points=10,
            slope_max_decades=5.0,
            step_policy=StepPolicy(
                step_default=0.1,
                step_mid=0.05,
                step_low=0.025,
                fer_mid=1e-2,
                fer_low=1e-4,
            ),
            fit_enable=True,
            fit_fer_hi=1.0,
            fit_fer_lo=1e-4,
            fit_target_fer=1e-6,
            fit_stop_margin=0.0,
            fit_min_points=2,
            fit_max_extend=1.0,
            gate_enable=False,
            gate_zero_k=0,
        )

    def test_completion_stage_exports_progress_when_jobs_are_in_flight(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            case = self._make_case(root)
            log_path = (Path(case.out_dir) / case.name / "fer_case_snr1.log").resolve()
            manifest = {
                "runs": [
                    {
                        "axis_type": "snr",
                        "axis_value": 1.0,
                        "log_path": str(log_path),
                        "cmd": ["/bin/echo", "placeholder"],
                        "ldpc_fer": 1.0,
                        "fail_cw": 3,
                        "total_packets": 3,
                    }
                ]
            }
            executor = _FakeExecutor()
            exports: list[dict[str, object]] = []

            def _capture_export(
                c: afe.CaseConfig,
                manifest_arg: dict[str, object],
                *,
                xlsx_path: Path,
                job_db_path: Path,
                planned_snrs: list[float],
                pending_snrs: list[float],
                in_flight_snrs: list[float],
            ) -> None:
                del c, manifest_arg, job_db_path
                exports.append(
                    {
                        "xlsx_path": str(xlsx_path),
                        "planned": list(planned_snrs),
                        "pending": list(pending_snrs),
                        "in_flight": list(in_flight_snrs),
                    }
                )

            with patch.object(afe, "export_progress_xlsx", side_effect=_capture_export), patch.object(
                afe.time, "sleep", return_value=None
            ):
                completed = afe.adaptive_finalize_incomplete_logs(
                    case,
                    manifest,
                    executor=executor,
                    poll_sec=0.2,
                    fail_fast=True,
                    timeout_log_grace_sec=0.0,
                    required_fail_cw=10,
                    max_in_flight=1,
                    queue="",
                    stages={"main"},
                    export_xlsx_sec=3600.0,
                )

            self.assertEqual(completed, 1)
            self.assertTrue(exports, "expected completion stage to trigger at least one progress export")
            self.assertTrue(
                any(1.0 in snapshot["in_flight"] for snapshot in exports),
                f"expected an export snapshot with in-flight completion jobs, got {exports}",
            )

    def test_completion_stage_writes_progress_xlsx_file(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            case = self._make_case(root)
            log_path = (Path(case.out_dir) / case.name / "fer_case_snr1.log").resolve()
            manifest = {
                "runs": [
                    {
                        "axis_type": "snr",
                        "axis_value": 1.0,
                        "log_path": str(log_path),
                        "cmd": ["/bin/echo", "placeholder"],
                        "ldpc_fer": 1.0,
                        "fail_cw": 3,
                        "total_packets": 3,
                    }
                ]
            }
            executor = _FakeExecutor()

            with patch.object(afe.time, "sleep", return_value=None):
                completed = afe.adaptive_finalize_incomplete_logs(
                    case,
                    manifest,
                    executor=executor,
                    poll_sec=0.2,
                    fail_fast=True,
                    timeout_log_grace_sec=0.0,
                    required_fail_cw=10,
                    max_in_flight=1,
                    queue="",
                    stages={"main"},
                    export_xlsx_sec=3600.0,
                )

            self.assertEqual(completed, 1)
            self.assertTrue(
                afe._progress_xlsx_path_for_case(case).exists(),
                "expected completion stage to materialize the progress workbook",
            )


if __name__ == "__main__":
    unittest.main()
