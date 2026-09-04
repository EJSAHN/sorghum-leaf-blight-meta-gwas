from __future__ import annotations

import hashlib
import json
import logging
import platform
import re
import subprocess
import sys
import time
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ID_ALIASES = {
    "BTX6": "BTX623",
    "BTX623": "BTX623",
    "SC748": "SC748_5",
    "SC7485": "SC748_5",
    "SC748_5": "SC748_5",
}


def normalize_id(value: object) -> str:
    text = re.sub(r"\s+", "", str(value)).replace("-", "_").upper()
    return ID_ALIASES.get(text, text)


def slugify(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def rank_int(values: pd.Series | np.ndarray) -> np.ndarray:
    series = pd.Series(values, dtype="float64")
    valid = series.notna()
    out = np.full(len(series), np.nan, dtype=float)
    if valid.sum() == 0:
        return out
    ranks = series[valid].rank(method="average")
    probs = (ranks - 0.5) / valid.sum()
    out[np.flatnonzero(valid.to_numpy())] = stats.norm.ppf(probs.clip(1e-12, 1 - 1e-12))
    return out


def bh_fdr(p_values: np.ndarray | pd.Series) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    safe = np.where(np.isfinite(p), np.clip(p, 0.0, 1.0), 1.0)
    return multipletests(safe, method="fdr_bh")[1]


def lambda_gc_from_p(p_values: np.ndarray | pd.Series) -> float:
    p = np.asarray(p_values, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size == 0:
        return float("nan")
    chi = stats.chi2.isf(np.clip(p, 1e-300, 1.0), df=1)
    return float(np.median(chi) / stats.chi2.ppf(0.5, df=1))


def safe_neglog10(p: np.ndarray | pd.Series) -> np.ndarray:
    arr = np.asarray(p, dtype=float)
    return -np.log10(np.clip(arr, 1e-300, 1.0))


def sha256_file(path: str | Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest().upper()


def bytes_hash(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest().upper()


def dataframe_hash(df: pd.DataFrame) -> str:
    values = pd.util.hash_pandas_object(df, index=True).to_numpy(dtype=np.uint64)
    return bytes_hash(values.tobytes())


def save_json(path: str | Path, obj: Any) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, ensure_ascii=False, default=json_default)


def load_json(path: str | Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, pd.Timestamp):
        return obj.isoformat()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def setup_logger(log_path: str | Path, verbose: bool = True) -> logging.Logger:
    logger = logging.getLogger("leaf_blight_final_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    if verbose:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
    return logger


def load_config(path: str | Path) -> dict[str, Any]:
    return load_json(path)


def ensure_directories(project_root: str | Path) -> dict[str, Path]:
    root = Path(project_root).resolve()
    mapping = {
        "root": root,
        "submitted": root / "00_submitted_snapshot",
        "inputs": root / "01_inputs",
        "pipeline": root / "02_pipeline",
        "config": root / "03_config",
        "intermediate": root / "04_intermediate",
        "results": root / "05_results",
        "tables": root / "06_tables",
        "figures": root / "07_figures",
        "logs": root / "08_logs",
        "manuscript": root / "09_manuscript",
        "response": root / "10_response_letter",
        "release": root / "11_release",
        "cache": root / "04_intermediate" / "genotype_cache",
        "structure": root / "04_intermediate" / "structure_cache",
        "gwas_cache": root / "04_intermediate" / "gwas_cache",
        "full_results": root / "05_results" / "full_scan_results",
    }
    for path in mapping.values():
        path.mkdir(parents=True, exist_ok=True)
    return mapping


def elapsed(start_time: float) -> str:
    seconds = max(0.0, time.time() - start_time)
    if seconds < 60:
        return f"{seconds:.1f} s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)} min {sec:.0f} s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)} h {int(minutes)} min"


def write_text(path: str | Path, text: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8")


def scientific(value: Any, digits: int = 4) -> str:
    try:
        return f"{float(value):.{digits}e}"
    except Exception:
        return ""

def software_versions() -> dict[str, str]:
    packages = [
        "numpy", "pandas", "scipy", "scikit-learn", "statsmodels",
        "openpyxl", "matplotlib",
    ]
    versions: dict[str, str] = {
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
    }
    for package in packages:
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def git_commit(source_root: str | Path) -> str:
    try:
        proc = subprocess.run(
            ["git", "-C", str(Path(source_root)), "rev-parse", "HEAD"],
            check=False, capture_output=True, text=True, encoding="utf-8", errors="replace",
        )
        if proc.returncode == 0:
            return proc.stdout.strip()
    except OSError:
        pass
    return "unavailable"
