from __future__ import annotations

import logging
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from lb_common import elapsed, save_json
from lb_structure import PopulationStructure


def windows_to_wsl(path: str | Path) -> str:
    raw = str(path)
    if len(raw) >= 3 and raw[1:3] == ":\\":
        drive = raw[0].lower()
        rest = raw[3:].replace("\\", "/")
        return f"/mnt/{drive}/{rest}"
    return raw.replace("\\", "/")


def _run_wsl(shell_command: str, log_path: Path, logger: logging.Logger) -> int:
    logger.info("COMMAND: wsl.exe -e sh -lc %s", shell_command)
    start = time.time()
    proc = subprocess.run(
        ["wsl.exe", "-e", "sh", "-lc", shell_command],
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "COMMAND\n" + shell_command + "\n\nSTDOUT\n" + proc.stdout + "\n\nSTDERR\n" + proc.stderr,
        encoding="utf-8",
    )
    logger.info("Exit %d (%s): %s", proc.returncode, elapsed(start), log_path)
    return int(proc.returncode)


def prepare_kinship_text(
    structure: PopulationStructure,
    output_dir: str | Path,
    logger: logging.Logger,
) -> dict[str, Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    rows: list[dict[str, Any]] = []

    def write_matrix(label: str, matrix: np.ndarray, path: Path) -> None:
        k = 0.5 * (np.asarray(matrix, dtype=float) + np.asarray(matrix, dtype=float).T)
        symmetry = float(np.max(np.abs(k - k.T)))
        minimum_eigenvalue = float(np.linalg.eigvalsh(k).min())
        if symmetry > 1e-8 or minimum_eigenvalue < -1e-7:
            raise ValueError(
                f"Kinship matrix {label} failed validation: symmetry={symmetry}, "
                f"min eigenvalue={minimum_eigenvalue}"
            )
        np.savetxt(path, k, fmt="%.15g")
        paths[label] = path
        rows.append(
            {
                "context": label,
                "path": str(path),
                "n_samples": int(k.shape[0]),
                "diag_mean": float(np.mean(np.diag(k))),
                "symmetry_error": symmetry,
                "minimum_eigenvalue": minimum_eigenvalue,
            }
        )

    write_matrix("global", structure.global_kinship, out / "global_kinship.cXX.txt")
    for chrom in range(1, 11):
        u = structure.loco_eigenvectors[chrom - 1]
        d = np.maximum(structure.loco_eigenvalues[chrom - 1], 0.0)
        k = (u * d[None, :]) @ u.T
        write_matrix(f"chr{chrom}", k, out / f"loco_chr{chrom}.cXX.txt")

    pd.DataFrame(rows).to_csv(out / "kinship_manifest.csv", index=False)
    logger.info("Prepared global and 10 chromosome-specific kinship matrices")
    return paths


def write_covariates(
    structure: PopulationStructure,
    output_path: str | Path,
    n_pcs: int,
) -> Path:
    if n_pcs < 0 or n_pcs > structure.pcs.shape[1]:
        raise ValueError(f"Requested {n_pcs} PCs, but only {structure.pcs.shape[1]} are available")
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame({"Intercept": np.ones(len(structure.sample_ids), dtype=int)})
    for idx in range(n_pcs):
        table[f"PC{idx + 1}"] = structure.pcs.reindex(structure.sample_ids).iloc[:, idx].to_numpy(float)
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        table.to_csv(handle, sep="\t", header=False, index=False, float_format="%.15g", lineterminator="\n")
    return path


def read_assoc(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    df = pd.read_csv(p, sep=r"\s+", engine="python", keep_default_na=False)
    rename: dict[str, str] = {}
    for col in df.columns:
        low = str(col).strip().lower()
        if low in {"rs", "snp", "id", "marker"}:
            rename[col] = "gemma_id"
        elif low in {"chr", "chrom", "chromosome"}:
            rename[col] = "CHR_out"
        elif low in {"ps", "pos", "position", "bp"}:
            rename[col] = "pos_out"
        elif low in {"p_wald", "pwald"}:
            rename[col] = "p_wald"
        elif low in {"p_lrt", "plrt"}:
            rename[col] = "p_lrt"
        elif low in {"p_score", "pscore"}:
            rename[col] = "p_score"
        elif low == "beta":
            rename[col] = "beta"
        elif low in {"se", "se_beta", "stderr"}:
            rename[col] = "se"
        elif low in {"l_remle", "lambda", "lambda_reml"}:
            rename[col] = "l_remle"
        elif low in {"logl_h1", "logl1"}:
            rename[col] = "logl_H1"
        elif low in {"af", "allele_freq", "maf"}:
            rename[col] = "af"
        elif low in {"n_miss", "nmiss"}:
            rename[col] = "n_miss"
    out = df.rename(columns=rename)
    if "gemma_id" not in out.columns:
        raise ValueError(f"No marker ID column in {p}; columns={out.columns.tolist()}")
    keep = [
        c
        for c in [
            "gemma_id", "CHR_out", "pos_out", "n_miss", "af", "beta", "se",
            "logl_H1", "l_remle", "p_wald", "p_lrt", "p_score",
        ]
        if c in out.columns
    ]
    out = out[keep].copy()
    out["gemma_id"] = out["gemma_id"].astype(str)
    for col in keep:
        if col != "gemma_id":
            out[col] = pd.to_numeric(out[col], errors="coerce")
    for col in ["p_wald", "p_lrt", "p_score"]:
        if col in out.columns:
            out[col] = out[col].clip(1e-300, 1.0)
    return out


def _find_assoc(work_dir: Path, prefix: str) -> Path:
    candidates = [
        work_dir / "output" / f"{prefix}.assoc.txt",
        work_dir / f"{prefix}.assoc.txt",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = list(work_dir.rglob(f"{prefix}.assoc.txt"))
    if not matches:
        raise FileNotFoundError(f"GEMMA association output not found for {prefix} in {work_dir}")
    return matches[0]


def map_assoc(assoc: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    base = mapping[["variant_index", "gemma_id", "snp", "CHR", "pos", "MAF", "missing_rate"]].copy()
    out = base.merge(assoc, on="gemma_id", how="inner", validate="one_to_one")
    if len(out) != len(assoc):
        raise ValueError(f"Mapped {len(out):,}/{len(assoc):,} GEMMA rows")
    return out.sort_values("variant_index", kind="mergesort").reset_index(drop=True)


def _run_one(
    gemma_wsl: str,
    prefix: Path,
    kinship: Path,
    output_dir: Path,
    output_name: str,
    logger: logging.Logger,
    lmm_mode: int,
    snp_list: Path | None = None,
    covariates: Path | None = None,
    force: bool = False,
) -> Path:
    expected = output_dir / "output" / f"{output_name}.assoc.txt"
    if expected.exists() and not force:
        logger.info("Reusing GEMMA output: %s", expected)
        return expected
    output_dir.mkdir(parents=True, exist_ok=True)
    parts = [
        f"cd {shlex.quote(windows_to_wsl(output_dir))}",
        f"{shlex.quote(gemma_wsl)}",
        f"-bfile {shlex.quote(windows_to_wsl(prefix))}",
        f"-k {shlex.quote(windows_to_wsl(kinship))}",
        f"-lmm {int(lmm_mode)}",
        "-maf 0 -miss 1",
    ]
    if snp_list is not None:
        parts.append(f"-snps {shlex.quote(windows_to_wsl(snp_list))}")
    if covariates is not None:
        parts.append(f"-c {shlex.quote(windows_to_wsl(covariates))}")
    parts.append(f"-o {shlex.quote(output_name)}")
    shell = " && ".join(parts[:2]) + " " + " ".join(parts[2:])
    log_path = output_dir.parent.parent / "logs" / f"{output_name}.log.txt"
    code = _run_wsl(shell, log_path, logger)
    if code != 0:
        raise RuntimeError(f"GEMMA exited with code {code}; see {log_path}")
    return _find_assoc(output_dir, output_name)


def run_gemma_scans(
    gemma_wsl: str,
    prefixes: dict[str, Path],
    mapping_path: str | Path,
    kinship_paths: dict[str, Path],
    chromosome_lists: dict[int, Path],
    output_dir: str | Path,
    logger: logging.Logger,
    lmm_mode: int = 4,
    run_global: bool = True,
    run_loco: bool = True,
    covariates: Path | None = None,
    force: bool = False,
) -> tuple[dict[tuple[str, str], pd.DataFrame], pd.DataFrame]:
    """Run official GEMMA and return mapped global/LOCO association tables."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(mapping_path)
    results: dict[tuple[str, str], pd.DataFrame] = {}
    status_rows: list[dict[str, Any]] = []

    for trait, prefix in prefixes.items():
        if run_global:
            work = out / "global" / trait.lower()
            name = f"gemma_global_{trait.lower()}_lmm{lmm_mode}"
            start = time.time()
            try:
                assoc_path = _run_one(
                    gemma_wsl, prefix, kinship_paths["global"], work, name,
                    logger, lmm_mode, covariates=covariates, force=force,
                )
                table = map_assoc(read_assoc(assoc_path), mapping)
                results[("global", trait)] = table
                status_rows.append(
                    {
                        "context": "global", "trait": trait, "status": "OK",
                        "n_rows": len(table), "elapsed": elapsed(start), "path": str(assoc_path),
                    }
                )
            except Exception as exc:
                status_rows.append(
                    {"context": "global", "trait": trait, "status": "FAILED", "n_rows": 0, "detail": str(exc)}
                )
                raise

        if run_loco:
            parts: list[pd.DataFrame] = []
            start = time.time()
            for chrom in range(1, 11):
                work = out / "loco" / trait.lower() / f"chr{chrom}"
                name = f"gemma_loco_{trait.lower()}_chr{chrom}_lmm{lmm_mode}"
                assoc_path = _run_one(
                    gemma_wsl, prefix, kinship_paths[f"chr{chrom}"], work, name,
                    logger, lmm_mode, snp_list=chromosome_lists[chrom],
                    covariates=covariates, force=force,
                )
                part = map_assoc(read_assoc(assoc_path), mapping)
                part = part.loc[part["CHR"].astype(int) == chrom].copy()
                parts.append(part)
            merged = pd.concat(parts, ignore_index=True).sort_values("variant_index", kind="mergesort")
            if len(merged) != len(mapping):
                raise ValueError(
                    f"LOCO output for {trait} has {len(merged):,} rows; expected {len(mapping):,}"
                )
            results[("loco", trait)] = merged.reset_index(drop=True)
            status_rows.append(
                {
                    "context": "loco", "trait": trait, "status": "OK",
                    "n_rows": len(merged), "elapsed": elapsed(start),
                    "path": str(out / "loco" / trait.lower()),
                }
            )

    status = pd.DataFrame(status_rows)
    status.to_csv(out / "gemma_run_status.csv", index=False)
    return results, status


def save_scan_tables(
    scans: dict[tuple[str, str], pd.DataFrame],
    output_dir: str | Path,
) -> dict[tuple[str, str], Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[tuple[str, str], Path] = {}
    for (context, trait), table in scans.items():
        path = out / f"gemma_{context}_{trait.lower()}_alltests.csv.gz"
        table.to_csv(path, index=False, compression="gzip")
        paths[(context, trait)] = path
    save_json(
        out / "scan_file_manifest.json",
        {f"{context}:{trait}": str(path) for (context, trait), path in paths.items()},
    )
    return paths


def load_scan_tables(output_dir: str | Path) -> dict[tuple[str, str], pd.DataFrame]:
    out = Path(output_dir)
    scans: dict[tuple[str, str], pd.DataFrame] = {}
    for path in out.glob("gemma_*_alltests.csv.gz"):
        stem = path.name.removeprefix("gemma_").removesuffix("_alltests.csv.gz")
        context, trait = stem.split("_", 1)
        scans[(context, trait.upper())] = pd.read_csv(path)
    return scans
