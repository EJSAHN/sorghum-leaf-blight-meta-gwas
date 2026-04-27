#!/usr/bin/env python
"""Excel-only dual-trait ACAT meta-GWAS pipeline for sorghum leaf blight."""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import t as tdist
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

try:
    import allel
except ImportError as exc:
    raise ImportError(
        "scikit-allel is required. Install with: conda env create -f environment.yml"
    ) from exc


@dataclass
class PipelineConfig:
    base_dir: str = "."
    phenotype_file: str = "phenotype.xlsx"
    phenotype_sheet: str = "auto"
    vcf_file: str = "700k.vcf"
    gff_file: str = "Sbicolor_454_v3.1.1.gene.gff3.gz"
    output_file: str = "Supplementary_Data_1.xlsx"
    n_pc: int = 3
    seed: int = 1
    maf_min: float = 0.05
    missing_max: float = 0.20
    fdr_cutoff: float = 0.10
    lead_bin_kb: int = 250
    gene_flank_kb: int = 50
    credible_window_kb: int = 300
    credible_coverage: float = 0.95
    single_trait_top_n: int = 10


def std_id(x: object) -> str:
    s = re.sub(r"\s+", "", str(x))
    return s.replace("-", "_")


def rank_int(y: pd.Series) -> pd.Series:
    y = pd.Series(y, dtype="float")
    r = y.rank(method="average", na_option="keep")
    p = (r - 0.5) / r.notna().sum()
    z = stats.norm.ppf(p.clip(1e-12, 1 - 1e-12))
    return pd.Series(z, index=y.index)


def chr_to_int(values: Iterable[object]) -> pd.Series:
    s = pd.Series(values).astype(str)
    s = s.str.replace(r"^chr", "", regex=True, flags=re.IGNORECASE).str.lstrip("0")
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def acat_p(pvec: Iterable[float]) -> float:
    p = np.asarray(pvec, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p < 1)]
    if p.size == 0:
        return np.nan
    t = np.tan((0.5 - p) * np.pi)
    return float(0.5 - np.arctan(np.mean(t)) / np.pi)


def open_text(path: str):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")


def has_vcf_header(path: str, max_lines: int = 2000) -> bool:
    try:
        with open_text(path) as handle:
            for _ in range(max_lines):
                line = handle.readline()
                if not line:
                    break
                if line.startswith("#CHROM") and "\t" in line:
                    return True
    except Exception:
        return False
    return False


def resolve_vcf(base_dir: str, vcf_file: str) -> str:
    candidate = os.path.join(base_dir, vcf_file)
    if os.path.isfile(candidate) and has_vcf_header(candidate):
        return candidate
    stem = os.path.splitext(candidate)[0]
    for ext in [".vcf", ".vcf.gz", ".vcf.bgz"]:
        p = stem + ext
        if os.path.isfile(p) and has_vcf_header(p):
            return p
    for pattern in ["*.vcf", "*.vcf.gz", "*.vcf.bgz"]:
        for p in Path(base_dir).rglob(pattern):
            if has_vcf_header(str(p)):
                return str(p)
    raise FileNotFoundError(f"No readable VCF found from {candidate}")


def read_phenotypes(cfg: PipelineConfig) -> pd.DataFrame:
    path = os.path.join(cfg.base_dir, cfg.phenotype_file)
    xls = pd.ExcelFile(path)
    sheet = xls.sheet_names[0] if cfg.phenotype_sheet.lower() in {"auto", "first"} else cfg.phenotype_sheet
    df = pd.read_excel(path, sheet_name=sheet)
    df.columns = [str(c).strip().replace("\n", " ") for c in df.columns]

    id_col = next(
        (c for c in df.columns if re.search(r"\b(cultivar|accession|id|genotype|line)\b", c, re.I)),
        None,
    )
    inc_col = next((c for c in df.columns if re.search(r"\binciden", c, re.I)), None)
    sev_col = next((c for c in df.columns if re.search(r"\bsever", c, re.I)), None)
    if id_col is None or inc_col is None or sev_col is None:
        raise ValueError("Could not identify ID, incidence, and severity columns in phenotype file.")

    out = pd.DataFrame(
        {
            "ID_std": df[id_col].map(std_id),
            "INC": pd.to_numeric(df[inc_col], errors="coerce"),
            "SEV": pd.to_numeric(df[sev_col], errors="coerce"),
        }
    )
    out["INC_INT"] = rank_int(out["INC"])
    out["SEV_INT"] = rank_int(out["SEV"])
    out = out.drop_duplicates("ID_std").set_index("ID_std")
    return out


def read_genotypes_and_pcs(cfg: PipelineConfig, pheno: pd.DataFrame):
    vcf_path = resolve_vcf(cfg.base_dir, cfg.vcf_file)
    fields = ["samples", "variants/CHROM", "variants/POS", "variants/ID", "calldata/GT"]
    vcf = allel.read_vcf(vcf_path, fields=fields)
    if vcf is None:
        raise RuntimeError(f"Could not read VCF: {vcf_path}")

    samples = np.array([std_id(s) for s in vcf["samples"]], dtype=object)
    keep_sample_mask = np.array([s in set(pheno.index) for s in samples], dtype=bool)
    kept_samples = [s for s, keep in zip(samples, keep_sample_mask) if keep]
    if len(kept_samples) < 30:
        raise ValueError(f"Too few matched phenotype-genotype samples: {len(kept_samples)}")

    chroms_all = np.array([str(c) for c in vcf["variants/CHROM"]], dtype=object)
    poss_all = np.array(vcf["variants/POS"], dtype=int)
    ids_raw = vcf["variants/ID"]
    chr_int_all = chr_to_int(chroms_all)
    snp_ids_all = np.array(
        [
            i if (i not in (None, ".")) else f"S{str(ch).zfill(2)}_{int(pos)}"
            for ch, pos, i in zip(chr_int_all.astype(object), poss_all, ids_raw)
        ],
        dtype=object,
    )

    gt_all = allel.GenotypeArray(vcf["calldata/GT"])
    gt = gt_all.compress(keep_sample_mask, axis=1)

    ac_all = gt.count_alleles()
    biallelic = np.asarray(ac_all.max_allele() == 1)
    gt_bi = gt[biallelic]

    ac_bi = gt_bi.count_alleles(max_allele=1)
    total_alleles = ac_bi.sum(axis=1).astype(float)
    alt_count = ac_bi[:, 1].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        maf = np.where(total_alleles > 0, np.minimum(alt_count, total_alleles - alt_count) / total_alleles, 0.0)

    called_raw = gt_bi.is_called()
    called = called_raw if called_raw.ndim == 2 else called_raw.all(axis=2)
    missing = 1.0 - called.mean(axis=1)
    pass_mask = (maf >= cfg.maf_min) & (missing <= cfg.missing_max)
    if int(pass_mask.sum()) == 0:
        raise ValueError("No SNPs passed QC. Consider reducing MAF or missingness thresholds.")

    gt_ok = gt_bi[pass_mask]
    chroms = chr_int_all[biallelic].to_numpy()[pass_mask].astype(int)
    poss = poss_all[biallelic][pass_mask].astype(int)
    snp_ids = snp_ids_all[biallelic][pass_mask]

    X = gt_ok.to_n_alt(fill=-1).T.astype("float32")
    X[X < 0] = np.nan
    col_mean = np.nanmean(X, axis=0)
    col_mean = np.where(np.isfinite(col_mean), col_mean, 0.0)
    X = np.where(np.isnan(X), col_mean, X).astype("float32")

    pheno_aligned = pheno.reindex(kept_samples)
    rng = np.random.default_rng(cfg.seed)
    pc = PCA(n_components=min(cfg.n_pc, X.shape[0] - 2), random_state=cfg.seed)
    pcs = pc.fit_transform(X - X.mean(axis=0, keepdims=True))
    pc_df = pd.DataFrame(pcs, index=kept_samples, columns=[f"PC{i+1}" for i in range(pcs.shape[1])])
    return X, chroms, poss, snp_ids, kept_samples, pheno_aligned, pc_df


def gwas_ols(trait: str, X: np.ndarray, pheno: pd.DataFrame, pc_df: pd.DataFrame, chroms, poss, snp_ids) -> pd.DataFrame:
    y = pd.to_numeric(pheno[trait], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(y)
    y = y[keep]
    Xk = X[keep, :]
    C = np.column_stack([np.ones(len(y)), pc_df.to_numpy()[keep, :]])

    beta_c = np.linalg.lstsq(C, y, rcond=None)[0]
    y_res = y - C @ beta_c
    beta_x = np.linalg.lstsq(C, Xk, rcond=None)[0]
    X_res = Xk - C @ beta_x

    numerator = X_res.T @ y_res
    xnorm = np.sqrt(np.sum(X_res ** 2, axis=0))
    ynorm = math.sqrt(float(np.sum(y_res ** 2)))
    corr = numerator / (xnorm * ynorm + 1e-12)
    df = max(len(y) - C.shape[1] - 1, 1)
    tstat = corr * np.sqrt(df / (1 - corr ** 2 + 1e-12))
    pvals = 2 * tdist.sf(np.abs(tstat), df)

    out = pd.DataFrame({"snp": snp_ids, "CHR": chroms, "pos": poss, "p": pvals})
    out = out[out["CHR"].between(1, 10)].sort_values(["CHR", "pos"]).reset_index(drop=True)
    return out


def run_single_trait_and_acat(cfg: PipelineConfig, X, chroms, poss, snp_ids, pheno, pc_df):
    res_inc = gwas_ols("INC_INT", X, pheno, pc_df, chroms, poss, snp_ids)
    res_sev = gwas_ols("SEV_INT", X, pheno, pc_df, chroms, poss, snp_ids)
    acat_df = res_inc[["snp", "CHR", "pos", "p"]].rename(columns={"p": "p_INC_INT"})
    acat_df = acat_df.merge(res_sev[["snp", "p"]].rename(columns={"p": "p_SEV_INT"}), on="snp", how="left")
    acat_df["p_ACAT"] = [acat_p(v) for v in acat_df[["p_INC_INT", "p_SEV_INT"]].to_numpy()]
    acat_df["q_BH"] = multipletests(acat_df["p_ACAT"].fillna(1.0), method="fdr_bh")[1]
    return res_inc, res_sev, acat_df


def select_leads(cfg: PipelineConfig, acat_df: pd.DataFrame) -> pd.DataFrame:
    candidates = acat_df.loc[acat_df["q_BH"] < cfg.fdr_cutoff].sort_values("p_ACAT")
    rows = []
    seen = set()
    bin_size = cfg.lead_bin_kb * 1000
    for _, row in candidates.iterrows():
        key = (int(row["CHR"]), int(row["pos"] // bin_size))
        if key in seen:
            continue
        seen.add(key)
        rows.append(row)
    if not rows:
        return pd.DataFrame(columns=["lead_snp", "CHR", "pos", "p_ACAT", "q_BH"])
    leads = pd.DataFrame(rows).rename(columns={"snp": "lead_snp"})
    return leads[["lead_snp", "CHR", "pos", "p_ACAT", "q_BH"]].reset_index(drop=True)


def credible_set_from_p(region: pd.DataFrame, coverage: float) -> pd.DataFrame:
    if region.empty:
        return pd.DataFrame(columns=["snp", "CHR", "pos", "p_ACAT", "post_w"])
    p = region["p_ACAT"].clip(1e-300, 1.0).astype(float)
    z = stats.norm.isf(p / 2.0) * np.sqrt(2.0)
    weights = np.exp(-0.5 * z ** 2)
    weights = np.asarray(weights, dtype=float)
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        weights = np.repeat(1.0 / len(region), len(region))
    else:
        weights = weights / weights.sum()
    order = np.argsort(-weights)
    cum = np.cumsum(weights[order])
    k = int(np.searchsorted(cum, coverage)) + 1
    out = region.iloc[order[:k]][["snp", "CHR", "pos", "p_ACAT"]].copy()
    out["post_w"] = weights[order[:k]]
    return out


def build_credible_sets(cfg: PipelineConfig, leads: pd.DataFrame, acat_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    window = cfg.credible_window_kb * 1000
    for _, lead in leads.iterrows():
        sub = acat_df.loc[
            (acat_df["CHR"].astype(int) == int(lead["CHR"]))
            & (acat_df["pos"].between(int(lead["pos"]) - window, int(lead["pos"]) + window))
        ].dropna(subset=["p_ACAT"])
        cs = credible_set_from_p(sub, cfg.credible_coverage)
        if not cs.empty:
            cs.insert(0, "lead_snp", lead["lead_snp"])
            rows.append(cs)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["lead_snp", "snp", "CHR", "pos", "p_ACAT", "post_w"])


def parse_gff_genes(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=["CHR", "start", "end", "gene"])
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, start, end, attrs = parts[0], parts[3], parts[4], parts[8]
            chr_int = chr_to_int([chrom]).iloc[0]
            if pd.isna(chr_int) or not (1 <= int(chr_int) <= 10):
                continue
            match = re.search(r"(?:^|;)ID=([^;]+)", attrs)
            if match:
                rows.append((int(chr_int), int(start), int(end), match.group(1)))
    return pd.DataFrame(rows, columns=["CHR", "start", "end", "gene"]).sort_values(["CHR", "start"])


def candidate_genes(cfg: PipelineConfig, leads: pd.DataFrame) -> pd.DataFrame:
    gff_path = os.path.join(cfg.base_dir, cfg.gff_file)
    genes = parse_gff_genes(gff_path)
    out_rows = []
    flank = cfg.gene_flank_kb * 1000
    for _, lead in leads.iterrows():
        text = ""
        if not genes.empty:
            sub = genes.loc[genes["CHR"] == int(lead["CHR"])].copy()
            if not sub.empty:
                pos = int(lead["pos"])
                dist = np.where(
                    (pos >= sub["start"].to_numpy()) & (pos <= sub["end"].to_numpy()),
                    0,
                    np.minimum(np.abs(pos - sub["start"].to_numpy()), np.abs(pos - sub["end"].to_numpy())),
                )
                sub = sub.assign(distance_bp=dist)
                keep = sub.loc[sub["distance_bp"] <= flank].sort_values(["distance_bp", "start"])
                text = "; ".join([f"{r.gene}({int(r.distance_bp)}bp)" for r in keep.itertuples()])
        out_rows.append(
            {
                "lead_snp": lead["lead_snp"],
                "CHR": int(lead["CHR"]),
                "pos": int(lead["pos"]),
                "nearest_genes_≤50kb": text,
            }
        )
    return pd.DataFrame(out_rows)


def single_trait_top(res_inc: pd.DataFrame, res_sev: pd.DataFrame, top_n: int) -> pd.DataFrame:
    inc = res_inc.sort_values("p").head(top_n).copy()
    inc.insert(0, "Trait", "INC_INT")
    sev = res_sev.sort_values("p").head(top_n).copy()
    sev.insert(0, "Trait", "SEV_INT")
    return pd.concat([inc, sev], ignore_index=True)[["Trait", "snp", "CHR", "pos", "p"]]


def miami_metrics(res_inc: pd.DataFrame, res_sev: pd.DataFrame, leads: pd.DataFrame, window_kb: int):
    merged = res_inc.rename(columns={"p": "p_INC"}).merge(
        res_sev[["snp", "p"]].rename(columns={"p": "p_SEV"}), on="snp", how="inner"
    )
    merged["lp_INC"] = -np.log10(merged["p_INC"].clip(1e-300))
    merged["lp_SEV"] = -np.log10(merged["p_SEV"].clip(1e-300))
    corr = pd.DataFrame(
        [
            {
                "scope": "genome",
                "spearman": merged[["lp_INC", "lp_SEV"]].corr(method="spearman").iloc[0, 1],
                "pearson": merged[["lp_INC", "lp_SEV"]].corr(method="pearson").iloc[0, 1],
                "N": len(merged),
            }
        ]
    )

    conc_rows = []
    N = len(merged)
    for frac in [0.01, 0.05, 0.10]:
        k = max(1, int(np.ceil(frac * N)))
        thr_i = merged["lp_INC"].nlargest(k).min()
        thr_s = merged["lp_SEV"].nlargest(k).min()
        A = merged["lp_INC"] >= thr_i
        B = merged["lp_SEV"] >= thr_s
        a = int((A & B).sum())
        b = int(A.sum() - a)
        c = int(B.sum() - a)
        d = int(N - a - b - c)
        table = np.array([[a, b], [c, d]])
        try:
            OR, fisher_p = stats.fisher_exact(table, alternative="two-sided")
        except Exception:
            OR, fisher_p = np.nan, np.nan
        conc_rows.append(
            {
                "frac": frac,
                "N": N,
                "topN": k,
                "overlap": a,
                "jaccard": a / (a + b + c) if (a + b + c) > 0 else np.nan,
                "fold": a / (N * frac * frac) if N > 0 else np.nan,
                "OR": OR,
                "fisher_p": fisher_p,
            }
        )
    concordance = pd.DataFrame(conc_rows)

    lead_rows = []
    window = window_kb * 1000
    for _, lead in leads.iterrows():
        ch = int(lead["CHR"])
        pos = int(lead["pos"])
        inc_sub = res_inc.loc[(res_inc["CHR"].astype(int) == ch) & (res_inc["pos"].between(pos - window, pos + window))]
        sev_sub = res_sev.loc[(res_sev["CHR"].astype(int) == ch) & (res_sev["pos"].between(pos - window, pos + window))]
        lead_rows.append(
            {
                "lead_snp": lead["lead_snp"],
                "CHR": ch,
                "pos": pos,
                "min_p_INC_window": float(inc_sub["p"].min()) if not inc_sub.empty else np.nan,
                "min_p_SEV_window": float(sev_sub["p"].min()) if not sev_sub.empty else np.nan,
            }
        )
    lead_window = pd.DataFrame(lead_rows)
    return concordance, corr, lead_window


def write_workbook(cfg: PipelineConfig, tables: dict[str, pd.DataFrame]) -> str:
    out_path = os.path.join(cfg.base_dir, cfg.output_file)
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, df in tables.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            ws = writer.book[sheet_name]
            ws.freeze_panes = "A2"
            for col in ws.columns:
                max_len = 0
                col_letter = col[0].column_letter
                for cell in col:
                    max_len = max(max_len, len(str(cell.value)) if cell.value is not None else 0)
                ws.column_dimensions[col_letter].width = min(max(max_len + 2, 10), 45)
    return out_path


def run_pipeline(cfg: PipelineConfig) -> str:
    np.random.seed(cfg.seed)
    pheno = read_phenotypes(cfg)
    X, chroms, poss, snp_ids, samples, pheno_aligned, pc_df = read_genotypes_and_pcs(cfg, pheno)
    res_inc, res_sev, acat_df = run_single_trait_and_acat(cfg, X, chroms, poss, snp_ids, pheno_aligned, pc_df)
    leads = select_leads(cfg, acat_df)
    lead_loci = leads[["lead_snp", "CHR", "pos", "p_ACAT"]].copy()
    genes = candidate_genes(cfg, leads)
    credible = build_credible_sets(cfg, leads, acat_df)
    trait_top = single_trait_top(res_inc, res_sev, cfg.single_trait_top_n)
    concordance, correlations, lead_window = miami_metrics(res_inc, res_sev, leads, cfg.credible_window_kb)
    tables = {
        "Lead_Loci_ACAT": lead_loci,
        "Candidate_Genes_50kb": genes,
        "Credible_Set_SNPs_95pct": credible,
        "Single_Trait_GWAS_TopLoci": trait_top,
        "Miami_Concordance": concordance,
        "Miami_Correlations": correlations,
        "Lead_Window_Concordance": lead_window,
    }
    return write_workbook(cfg, tables)


def parse_args() -> PipelineConfig:
    parser = argparse.ArgumentParser(description="Excel-only dual-trait ACAT meta-GWAS pipeline for sorghum leaf blight.")
    parser.add_argument("--base-dir", default=PipelineConfig.base_dir)
    parser.add_argument("--phenotype-file", default=PipelineConfig.phenotype_file)
    parser.add_argument("--phenotype-sheet", default=PipelineConfig.phenotype_sheet)
    parser.add_argument("--vcf-file", default=PipelineConfig.vcf_file)
    parser.add_argument("--gff-file", default=PipelineConfig.gff_file)
    parser.add_argument("--output-file", default=PipelineConfig.output_file)
    parser.add_argument("--n-pc", type=int, default=PipelineConfig.n_pc)
    parser.add_argument("--maf-min", type=float, default=PipelineConfig.maf_min)
    parser.add_argument("--missing-max", type=float, default=PipelineConfig.missing_max)
    parser.add_argument("--fdr-cutoff", type=float, default=PipelineConfig.fdr_cutoff)
    parser.add_argument("--lead-bin-kb", type=int, default=PipelineConfig.lead_bin_kb)
    parser.add_argument("--gene-flank-kb", type=int, default=PipelineConfig.gene_flank_kb)
    parser.add_argument("--credible-window-kb", type=int, default=PipelineConfig.credible_window_kb)
    parser.add_argument("--seed", type=int, default=PipelineConfig.seed)
    args = parser.parse_args()
    return PipelineConfig(**vars(args))


if __name__ == "__main__":
    config = parse_args()
    output = run_pipeline(config)
    print(f"Saved workbook: {output}")
