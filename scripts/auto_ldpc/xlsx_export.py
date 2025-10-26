"""
Export simulation data to XLSX.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path
from typing import Dict, Iterable, List

import json

try:
    from openpyxl import Workbook
except ImportError:  # pragma: no cover
    Workbook = None


HEADERS = [
    "Matrix_ID",
    "Model",
    "SNR",
    "FER",
    "Metrics_JSON",
    "Log_Path",
    "Updated_At",
]


def export_records(records: Iterable[Dict[str, object]], output_path: Path) -> None:
    if Workbook is None:
        raise RuntimeError("未安装 openpyxl，无法导出 XLSX。请运行 `pip install openpyxl`。")

    wb = Workbook()
    ws = wb.active
    ws.title = "results"
    ws.append(HEADERS)

    for record in records:
        metrics = record.get("metrics")
        if isinstance(metrics, (dict, list)):
            metrics_value = json.dumps(metrics, ensure_ascii=False)
        else:
            metrics_value = metrics
        ws.append(
            [
                record.get("matrix_id"),
                record.get("model"),
                record.get("snr"),
                record.get("fer"),
                metrics_value,
                record.get("log_path"),
                record.get("updated_at"),
            ]
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(output_path)
