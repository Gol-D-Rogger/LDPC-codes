"""
Configuration utilities for LDPC automation workflow.

The module loads configuration data from JSON or YAML files and merges it
with sane defaults so downstream code can rely on required keys existing.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import yaml  # type: ignore
except ImportError:  # pragma: no cover
    yaml = None


DEFAULT_CONFIG: Dict[str, Any] = {
    "system": {
        "command_prefix": [],
    },
    "paths": {
        "workspace_root": None,
        "output_root": None,
        "exec": "gen4_ldpc_sim/ssd_fc_dq",
        "config": "gen4_ldpc_sim/config/sdec_dq.cnfg",
        "gen_bin": "GenLDPC/ldpc_gen",
        "matrix_dir": None,
        "sim_out_root": "runs/auto_ldpc",
        "perf_dir": "perf",
        "log_suffix": ".log",
    },
    "generation": {
        "gen_n": None,
        "gen_k": None,
        "gen_count": None,
        "batch_size": 50,
        "queue": None,
        "mode": "lsf",
        "anchor_matrix_dir": None,
        "gen_out_base": None,
    },
    "simulation": {
        "queue": "regr_q",
        "job_prefix": "ldpc_auto",
        "snr_list": [],
        "anchor_snr": None,
        "deep_snr_list": [],
        "dry_run": False,
        "max_jobs": None,
    },
    "selection": {
        "top_n": 5,
        "err_inj_threshold": 1e-3,
    },
    "err_inj": {
        "config": "config/another.cnfg",
        "command": "./ssd_fc LDPC {config} ERR_INJ 340 {matrix_id}",
        "queue": "regr_q",
        "job_prefix": "ldpc_errinj",
    },
    "limits": {
        "max_matrices": 1000,
        "max_runtime_hours": None,
        "sleep_seconds": 120,
        "stop_flag_file": None,
    },
    "log_parser": {
        "awgn": {
            "stat_keys": [
                "Total simulated number",
                "RAW  BER",
                "LDPC BER",
                "LDPC FER",
                "MCRC FER",
                "DATA FER",
                "Retry Decoder average iterations",
            ],
            "fer_key": "LDPC FER",
            "statistics_prefix": "[STATISTICS]",
        },
        "err_inj": {
            "stat_keys": [
                "Total simulated number",
                "RAW  BER",
                "LDPC BER",
                "LDPC FER",
                "MCRC FER",
                "DATA FER",
                "Fast Decoder average iterations",
            ],
            "fer_key": "LDPC FER",
            "statistics_prefix": "[STATISTICS]",
        },
    },
    "xlsx": {
        "enable": True,
        "filename": "simulation_results.xlsx",
    },
    "database": {
        "path": "runs/auto_ldpc/state.db",
    },
    "plot": {
        "enable": True,
        "output": "perf/auto_ldpc/best_plot.png",
        "baseline": {
            "path": None,
            "sheet": "baseline",
            "snr_column": "SNR",
            "fer_column": "FER",
        },
    },
}


class ConfigError(RuntimeError):
    """Raised when the configuration file cannot be loaded or validated."""


def _deep_merge(base: Dict[str, Any], updates: Dict[str, Any]) -> Dict[str, Any]:
    merged: Dict[str, Any] = dict(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merged[key] = _deep_merge(base[key], value)
        else:
            merged[key] = value
    return merged


def _load_raw_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise ConfigError(f"配置文件不存在: {path}")

    if path.suffix.lower() in {".yaml", ".yml"}:
        if yaml is None:
            raise ConfigError("未安装 PyYAML，无法解析 YAML 配置。请运行 `pip install pyyaml`。")
        with path.open("r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
    else:
        with path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)

    if not isinstance(data, dict):
        raise ConfigError("配置文件内容应为对象/字典。")
    return data


def load_config(path: Path) -> Dict[str, Any]:
    """Load configuration file and apply defaults."""

    raw = _load_raw_config(path)
    merged = _deep_merge(DEFAULT_CONFIG, raw)
    return merged


def _make_abs(path_value: Optional[str], base: Path) -> Optional[Path]:
    if path_value in (None, ""):
        return None
    path = Path(path_value)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _expand_placeholders(value: Optional[str], replacements: Dict[str, str]) -> Optional[str]:
    if value is None:
        return None
    result = value
    for key, repl in replacements.items():
        placeholder = f"${{{key}}}"
        result = result.replace(placeholder, repl)
    return result


@dataclass
class ConfigPaths:
    workspace_root: Path
    output_root: Optional[Path]
    exec_path: Path
    config_path: Path
    gen_bin: Path
    matrix_dir: Optional[Path]
    sim_out_root: Path
    perf_dir: Path
    log_suffix: str
    rename_script: Path


def resolve_paths(cfg: Dict[str, Any], base_dir: Optional[Path] = None) -> ConfigPaths:
    paths_cfg = cfg["paths"]
    workspace = paths_cfg.get("workspace_root")
    workspace_root = Path(workspace) if workspace else (base_dir or Path.cwd())
    workspace_root = workspace_root.resolve()

    output_root_cfg = paths_cfg.get("output_root")
    output_root = None
    if output_root_cfg:
        expanded = _expand_placeholders(str(output_root_cfg), {"WORKSPACE_ROOT": str(workspace_root)})
        output_root = _make_abs(expanded, workspace_root)

    replacements = {
        "WORKSPACE_ROOT": str(workspace_root),
        "OUTPUT_ROOT": str(output_root) if output_root else "",
    }

    def to_abs(value: Optional[str], default: Optional[Path] = None) -> Optional[Path]:
        expanded = _expand_placeholders(value, replacements)
        if expanded:
            return _make_abs(expanded, workspace_root)
        return default

    exec_path = to_abs(paths_cfg.get("exec")) or workspace_root / "gen4_ldpc_sim/ssd_fc_dq"
    config_path = to_abs(paths_cfg.get("config")) or workspace_root / "gen4_ldpc_sim/config/sdec_dq.cnfg"
    gen_bin = to_abs(paths_cfg.get("gen_bin")) or workspace_root / "GenLDPC/ldpc_gen"
    matrix_dir_default = output_root / "matrix" if output_root else None
    matrix_dir = to_abs(paths_cfg.get("matrix_dir"), matrix_dir_default)
    sim_out_default = output_root / "runs" if output_root else workspace_root / "runs/auto_ldpc"
    sim_out_root = to_abs(paths_cfg.get("sim_out_root"), sim_out_default) or workspace_root / "runs/auto_ldpc"
    perf_default = output_root / "perf" if output_root else workspace_root / "perf"
    perf_dir = to_abs(paths_cfg.get("perf_dir"), perf_default) or workspace_root / "perf"
    log_suffix = paths_cfg.get("log_suffix") or ".log"

    rename_script = to_abs(paths_cfg.get("rename_script"), workspace_root / "scripts/rename_ldpc_matrices.sh")

    return ConfigPaths(
        workspace_root=workspace_root,
        output_root=output_root,
        exec_path=exec_path.resolve(),
        config_path=config_path.resolve(),
        gen_bin=gen_bin.resolve(),
        matrix_dir=matrix_dir.resolve() if matrix_dir else None,
        sim_out_root=sim_out_root.resolve(),
        perf_dir=perf_dir.resolve(),
        log_suffix=log_suffix,
        rename_script=rename_script.resolve(),
    )


def ensure_lists(cfg: Dict[str, Any]) -> None:
    sim_cfg = cfg.setdefault("simulation", {})
    if not isinstance(sim_cfg.get("snr_list"), list):
        sim_cfg["snr_list"] = list(sim_cfg["snr_list"]) if sim_cfg.get("snr_list") else []
    if not isinstance(sim_cfg.get("deep_snr_list"), list):
        sim_cfg["deep_snr_list"] = list(sim_cfg["deep_snr_list"]) if sim_cfg.get("deep_snr_list") else []

    log_cfg = cfg.setdefault("log_parser", {})
    for key in ("awgn", "err_inj"):
        section = log_cfg.get(key)
        if not isinstance(section, dict):
            continue
        stat_keys = section.get("stat_keys")
        if stat_keys is None:
            continue
        if not isinstance(stat_keys, list):
            section["stat_keys"] = list(stat_keys)

    plot_cfg = cfg.setdefault("plot", {})
    baseline_cfg = plot_cfg.get("baseline")
    if isinstance(baseline_cfg, dict):
        for key in ("path", "sheet", "snr_column", "fer_column"):
            baseline_cfg.setdefault(key, None)


@dataclass
class ConfigBundle:
    raw: Dict[str, Any]
    paths: ConfigPaths
    config_path: Path
    placeholders: Dict[str, str]
    base_dir: Path = field(default_factory=lambda: Path.cwd())

    def expand_string(self, value: Optional[str]) -> Optional[str]:
        return _expand_placeholders(value, self.placeholders)

    def expand_path(self, value: Optional[str], base: Optional[Path] = None) -> Optional[Path]:
        expanded = self.expand_string(value)
        if expanded is None or expanded == "":
            return None
        base_path = base or self.paths.output_root or self.paths.workspace_root
        path = Path(expanded)
        if not path.is_absolute():
            path = base_path / path
        return path.resolve()


def load_bundle(config_path: Path) -> ConfigBundle:
    cfg = load_config(config_path)
    ensure_lists(cfg)
    paths = resolve_paths(cfg, config_path.parent)
    placeholders = {
        "WORKSPACE_ROOT": str(paths.workspace_root),
        "OUTPUT_ROOT": str(paths.output_root) if paths.output_root else "",
    }
    return ConfigBundle(
        raw=cfg,
        paths=paths,
        config_path=config_path,
        placeholders=placeholders,
        base_dir=config_path.parent.resolve(),
    )
