from __future__ import annotations

import logging
import re
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from lb_common import normalize_id, rank_int, save_json


MISSING_TOKENS = {
    "": np.nan,
    ".": np.nan,
    "MISSING": np.nan,
    "missing": np.nan,
    "Missing": np.nan,
    "missintg": np.nan,
    "MISSINTG": np.nan,
    "NA": np.nan,
    "N/A": np.nan,
}


def _numeric(series: pd.Series) -> pd.Series:
    work = series.copy()
    text = work.astype(str).str.strip()
    missing_mask = work.isna() | text.isin(set(MISSING_TOKENS))
    work = work.mask(missing_mask)
    return pd.to_numeric(work, errors="coerce")


def _find_header_row(raw: pd.DataFrame, required_tokens: tuple[str, ...]) -> int:
    for idx in range(min(20, len(raw))):
        values = " | ".join(str(x).strip().lower() for x in raw.iloc[idx].tolist())
        if all(token.lower() in values for token in required_tokens):
            return idx
    raise ValueError(f"Could not find header row containing {required_tokens}")


def _find_lb_columns(raw: pd.DataFrame, header_row: int) -> tuple[int, int]:
    # In the original field workbooks, the row immediately above the column labels
    # contains the disease abbreviation "LB" over the incidence column, followed by
    # the severity column.
    for label_row in range(max(0, header_row - 3), header_row):
        for col in range(raw.shape[1] - 1):
            label = str(raw.iat[label_row, col]).strip().upper()
            h1 = str(raw.iat[header_row, col]).strip().lower()
            h2 = str(raw.iat[header_row, col + 1]).strip().lower()
            if label == "LB" and h1.startswith("inc") and h2.startswith("sev"):
                return col, col + 1
    raise ValueError("Could not identify leaf-blight incidence/severity columns")


def parse_niger_field(path: str | Path, location: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=0, header=None)
    header_row = _find_header_row(raw, ("plot", "line", "plants"))
    lb_inc_col, lb_sev_col = _find_lb_columns(raw, header_row)

    header = [str(x).strip().lower() for x in raw.iloc[header_row].tolist()]
    plot_col = next(i for i, x in enumerate(header) if "plot" in x)
    line_col = next(i for i, x in enumerate(header) if "line" in x)
    plants_col = next(i for i, x in enumerate(header) if "plant" in x)

    out = raw.iloc[header_row + 1 :, [plot_col, line_col, plants_col, lb_inc_col, lb_sev_col]].copy()
    out.columns = ["Plot", "Cultivar_raw", "Plants", "INC_raw_source", "SEV_raw_source"]
    out["Plot"] = _numeric(out["Plot"])
    out["Plants"] = _numeric(out["Plants"])
    out["INC_raw_source"] = _numeric(out["INC_raw_source"])
    out["SEV_raw_source"] = _numeric(out["SEV_raw_source"])
    out["Cultivar_raw"] = out["Cultivar_raw"].astype(str).str.strip()
    out = out[out["Plot"].notna() & out["Cultivar_raw"].ne("nan")].copy()
    out["Rep"] = (out["Plot"] // 1000).astype("Int64")
    out["Block"] = pd.Series(pd.NA, index=out.index, dtype="Int64")
    out["Location"] = location
    out["Country"] = "Niger"
    out["Year"] = 2022
    return out.reset_index(drop=True)


def parse_senegal_field(path: str | Path) -> pd.DataFrame:
    xls = pd.ExcelFile(path)
    frames: list[pd.DataFrame] = []
    for sheet in xls.sheet_names:
        raw = pd.read_excel(path, sheet_name=sheet, header=None)
        header_row = _find_header_row(raw, ("rep", "block", "plot", "entry"))
        lb_inc_col, lb_sev_col = _find_lb_columns(raw, header_row)

        header = [str(x).strip().lower() for x in raw.iloc[header_row].tolist()]
        rep_col = next(i for i, x in enumerate(header) if x == "rep")
        block_col = next(i for i, x in enumerate(header) if x == "block")
        plot_col = next(i for i, x in enumerate(header) if x == "plot")
        entry_col = next(i for i, x in enumerate(header) if x == "entry")
        plants_col = next(i for i, x in enumerate(header) if "plant" in x)

        location = None
        for value in raw.iloc[:header_row, :5].to_numpy().ravel():
            text = str(value).strip()
            if text.lower() in {"kolda", "ndiaganiao", "kaymor"}:
                location = text
                break
        if location is None:
            raise ValueError(f"Could not determine Senegal location in sheet {sheet!r}")

        out = raw.iloc[
            header_row + 1 :,
            [rep_col, block_col, plot_col, entry_col, plants_col, lb_inc_col, lb_sev_col],
        ].copy()
        out.columns = [
            "Rep",
            "Block",
            "Plot",
            "Cultivar_raw",
            "Plants",
            "INC_raw_source",
            "SEV_raw_source",
        ]
        for col in ["Rep", "Block", "Plot", "Plants", "INC_raw_source", "SEV_raw_source"]:
            out[col] = _numeric(out[col])
        out["Cultivar_raw"] = out["Cultivar_raw"].astype(str).str.strip()
        out = out[out["Cultivar_raw"].ne("nan")].copy()
        out["Location"] = location
        out["Country"] = "Senegal"
        out["Year"] = 2022
        frames.append(out.reset_index(drop=True))
    return pd.concat(frames, ignore_index=True)


def load_corrected_leaf_blight(path: str | Path) -> pd.DataFrame:
    df = pd.read_excel(path)
    required = {"Cultivar", "LB-Incidence", "LB-Severity", "Location", "Country", "Year"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Corrected phenotype workbook is missing columns: {sorted(missing)}")
    out = df[list(required)].copy()
    out = out[["Cultivar", "LB-Incidence", "LB-Severity", "Location", "Country", "Year"]]
    out["Cultivar"] = out["Cultivar"].astype(str).str.strip()
    out["ID_std"] = out["Cultivar"].map(normalize_id)
    out["INC"] = _numeric(out["LB-Incidence"])
    out["SEV"] = _numeric(out["LB-Severity"])
    out["Location"] = out["Location"].astype(str).str.strip()
    out["Country"] = out["Country"].astype(str).str.strip()
    out["Year"] = _numeric(out["Year"]).astype("Int64")
    return out[["Cultivar", "ID_std", "INC", "SEV", "Location", "Country", "Year"]]


def build_canonical_phenotype(inputs_dir: str | Path, intermediate_dir: str | Path, logger: logging.Logger) -> tuple[pd.DataFrame, dict[str, Any]]:
    inputs = Path(inputs_dir)
    corrected_path = inputs / "Phenotype.xlsx"
    raw_dir = inputs / "raw_field"
    maradi = raw_dir / "Maradi_Field_Niger_2022_MAY_2023.xlsx"
    bengou = raw_dir / "Bengou_Field_Niger_2022_MAY_18_2023.xlsx"
    senegal = raw_dir / "Field_data_all_locations_SEN_2022_MAY_18_2023.xlsx"

    corrected = load_corrected_leaf_blight(corrected_path)
    raw_available = all(path.exists() for path in (maradi, bengou, senegal))
    report: dict[str, Any] = {
        "corrected_rows": int(len(corrected)),
        "corrected_unique_ids": int(corrected["ID_std"].nunique()),
        "raw_field_files_available": raw_available,
    }

    if raw_available:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            raw = pd.concat(
                [
                    parse_niger_field(maradi, "Maradi"),
                    parse_niger_field(bengou, "Bengou"),
                    parse_senegal_field(senegal),
                ],
                ignore_index=True,
            )
        raw_complete = raw.dropna(subset=["INC_raw_source", "SEV_raw_source"]).reset_index(drop=True)
        corrected = corrected.dropna(subset=["INC", "SEV"]).reset_index(drop=True)

        if len(raw_complete) != len(corrected):
            raise ValueError(
                f"Raw complete-pair rows ({len(raw_complete)}) do not match corrected rows ({len(corrected)})."
            )
        location_match = raw_complete["Location"].astype(str).to_numpy() == corrected["Location"].astype(str).to_numpy()
        incidence_match = np.isclose(
            raw_complete["INC_raw_source"].to_numpy(float), corrected["INC"].to_numpy(float), equal_nan=True
        )
        raw_ids = raw_complete["Cultivar_raw"].map(normalize_id)
        id_match = raw_ids.to_numpy() == corrected["ID_std"].to_numpy()

        if not bool(location_match.all()):
            bad = np.flatnonzero(~location_match)[:10]
            raise ValueError(f"Raw/corrected location order mismatch at rows {bad.tolist()}")
        if not bool(incidence_match.all()):
            bad = np.flatnonzero(~incidence_match)[:10]
            raise ValueError(f"Raw/corrected incidence mismatch at rows {bad.tolist()}")
        if id_match.mean() < 0.99:
            bad = np.flatnonzero(~id_match)[:10]
            raise ValueError(f"Unexpected raw/corrected ID mismatch at rows {bad.tolist()}")

        canonical = raw_complete.copy()
        canonical["Cultivar"] = corrected["Cultivar"]
        canonical["ID_std"] = corrected["ID_std"]
        canonical["INC"] = corrected["INC"]
        canonical["SEV"] = corrected["SEV"]
        canonical["INC_correction"] = canonical["INC"] - canonical["INC_raw_source"]
        canonical["SEV_correction"] = canonical["SEV"] - canonical["SEV_raw_source"]
        canonical["CorrectionFlag"] = (
            canonical["INC_correction"].abs().fillna(0) > 1e-12
        ) | (canonical["SEV_correction"].abs().fillna(0) > 1e-12)
        report.update(
            {
                "raw_rows": int(len(raw)),
                "raw_complete_pairs": int(len(raw_complete)),
                "raw_corrected_id_match_fraction": float(id_match.mean()),
                "raw_corrected_value_corrections": int(canonical["CorrectionFlag"].sum()),
            }
        )
    else:
        logger.warning("Original field workbooks were not found; replicate/block metadata will be unavailable.")
        canonical = corrected.copy()
        canonical["Cultivar_raw"] = canonical["Cultivar"]
        canonical["INC_raw_source"] = canonical["INC"]
        canonical["SEV_raw_source"] = canonical["SEV"]
        canonical["Rep"] = pd.Series(pd.NA, index=canonical.index, dtype="Int64")
        canonical["Block"] = pd.Series(pd.NA, index=canonical.index, dtype="Int64")
        canonical["Plot"] = np.arange(1, len(canonical) + 1)
        canonical["Plants"] = np.nan
        canonical["INC_correction"] = 0.0
        canonical["SEV_correction"] = 0.0
        canonical["CorrectionFlag"] = False

    canonical["Rep"] = pd.to_numeric(canonical["Rep"], errors="coerce").astype("Int64")
    canonical["Block"] = pd.to_numeric(canonical["Block"], errors="coerce").astype("Int64")
    canonical["Plot"] = pd.to_numeric(canonical["Plot"], errors="coerce")
    canonical["Plants"] = pd.to_numeric(canonical["Plants"], errors="coerce")

    def make_stratum(row: pd.Series) -> str:
        location = str(row["Location"])
        rep = "NA" if pd.isna(row["Rep"]) else str(int(row["Rep"]))
        if pd.isna(row["Block"]):
            return f"{location}_R{rep}"
        return f"{location}_R{rep}_B{int(row['Block'])}"

    canonical["Stratum"] = canonical.apply(make_stratum, axis=1)
    canonical["DiseasePresent"] = canonical["INC"] > 0
    canonical["INC_zero_SEV_zero"] = (canonical["INC"] == 0) & (canonical["SEV"] == 0)
    canonical["INC_zero_SEV_positive"] = (canonical["INC"] == 0) & (canonical["SEV"] > 0)
    canonical["INC_positive_SEV_zero"] = (canonical["INC"] > 0) & (canonical["SEV"] == 0)

    report.update(
        {
            "canonical_rows": int(len(canonical)),
            "canonical_unique_ids": int(canonical["ID_std"].nunique()),
            "locations": {k: int(v) for k, v in canonical.groupby("Location").size().items()},
            "strata": int(canonical["Stratum"].nunique()),
            "inc0_sev0": int(canonical["INC_zero_SEV_zero"].sum()),
            "inc0_sev_positive": int(canonical["INC_zero_SEV_positive"].sum()),
            "inc_positive_sev0": int(canonical["INC_positive_SEV_zero"].sum()),
            "severity_unique_values": sorted(float(v) for v in canonical["SEV"].dropna().unique()),
        }
    )

    intermediate = Path(intermediate_dir)
    intermediate.mkdir(parents=True, exist_ok=True)
    canonical_path = intermediate / "phenotype_plot_level_canonical.csv"
    canonical.to_csv(canonical_path, index=False)
    save_json(intermediate / "phenotype_qc_report.json", report)
    logger.info(
        "Canonical phenotype: %s rows, %s IDs, %s locations, %s strata",
        len(canonical),
        canonical["ID_std"].nunique(),
        canonical["Location"].nunique(),
        canonical["Stratum"].nunique(),
    )
    return canonical, report


def _fit_adjusted_means(
    data: pd.DataFrame,
    value_col: str,
    covariate_cols: list[str] | None = None,
    label: str = "trait",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    covariate_cols = covariate_cols or []
    cols = ["ID_std", "Stratum", value_col, *covariate_cols]
    work = data[cols].copy().dropna(subset=["ID_std", "Stratum", value_col, *covariate_cols])
    work["ID_std"] = work["ID_std"].astype(str)
    work["Stratum"] = work["Stratum"].astype(str)
    accessions = sorted(work["ID_std"].unique())
    strata = sorted(work["Stratum"].unique())
    if len(accessions) < 10 or len(work) <= len(accessions) + len(strata):
        raise ValueError(f"Insufficient observations to fit adjusted means for {label}")

    acc_ref = accessions[0]
    stratum_ref = strata[0]
    acc_cols = accessions[1:]
    stratum_cols = strata[1:]

    n = len(work)
    p = 1 + len(acc_cols) + len(stratum_cols) + len(covariate_cols)
    X = np.zeros((n, p), dtype=float)
    X[:, 0] = 1.0
    col_names = ["Intercept"]
    offset = 1

    acc_map = {name: i for i, name in enumerate(acc_cols)}
    for row_idx, name in enumerate(work["ID_std"]):
        if name != acc_ref:
            X[row_idx, offset + acc_map[name]] = 1.0
    col_names.extend([f"ID[{name}]" for name in acc_cols])
    offset += len(acc_cols)

    stratum_map = {name: i for i, name in enumerate(stratum_cols)}
    for row_idx, name in enumerate(work["Stratum"]):
        if name != stratum_ref:
            X[row_idx, offset + stratum_map[name]] = 1.0
    col_names.extend([f"Stratum[{name}]" for name in stratum_cols])
    offset += len(stratum_cols)

    covariate_centers: dict[str, float] = {}
    for col in covariate_cols:
        values = pd.to_numeric(work[col], errors="coerce").to_numpy(float)
        center = float(np.mean(values))
        covariate_centers[col] = center
        X[:, offset] = values - center
        col_names.append(f"Covariate[{col}]_centered")
        offset += 1

    y = pd.to_numeric(work[value_col], errors="coerce").to_numpy(float)
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    residual = y - X @ beta
    df_resid = max(n - int(rank), 1)
    sigma2 = float(np.dot(residual, residual) / df_resid)
    xtx_inv = np.linalg.pinv(X.T @ X)
    cov_beta = sigma2 * xtx_inv

    # Equal-weighted marginal mean over all observed nuisance strata, with continuous
    # covariates fixed at their observed means (i.e. centered value zero).
    L = np.zeros((len(accessions), p), dtype=float)
    L[:, 0] = 1.0
    for i, name in enumerate(accessions):
        if name != acc_ref:
            L[i, 1 + acc_map[name]] = 1.0
    stratum_start = 1 + len(acc_cols)
    if len(strata) > 0:
        for name in stratum_cols:
            L[:, stratum_start + stratum_map[name]] = 1.0 / len(strata)

    estimate = L @ beta
    variance = np.einsum("ij,jk,ik->i", L, cov_beta, L)
    se = np.sqrt(np.maximum(variance, 0.0))
    n_obs = work.groupby("ID_std").size().reindex(accessions).fillna(0).astype(int)

    result = pd.DataFrame(
        {
            "ID_std": accessions,
            "estimate": estimate,
            "se": se,
            "ci_low": estimate - 1.96 * se,
            "ci_high": estimate + 1.96 * se,
            "n_plots": n_obs.to_numpy(),
        }
    )
    diagnostics = {
        "label": label,
        "value_col": value_col,
        "covariates": covariate_cols,
        "covariate_centers": covariate_centers,
        "n_observations": n,
        "n_accessions": len(accessions),
        "n_strata": len(strata),
        "design_columns": p,
        "design_rank": int(rank),
        "residual_df": int(df_resid),
        "residual_sd": float(np.sqrt(sigma2)),
        "accession_reference": acc_ref,
        "stratum_reference": stratum_ref,
    }
    return result, diagnostics


def build_accession_traits(canonical: pd.DataFrame, intermediate_dir: str | Path, logger: logging.Logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    accession_index = sorted(canonical["ID_std"].unique())
    traits = pd.DataFrame(index=accession_index)
    model_rows: list[dict[str, Any]] = []

    model_specs = [
        ("INC_ADJ_RAW", canonical, "INC", []),
        ("SEV_ADJ_RAW", canonical, "SEV", []),
        ("SEV_POS_ADJ_RAW", canonical.loc[canonical["INC"] > 0].copy(), "SEV", []),
        ("SEV_COND_ADJ_RAW", canonical, "SEV", ["INC"]),
    ]
    for trait_name, data, value_col, covariates in model_specs:
        estimates, diagnostics = _fit_adjusted_means(data, value_col, covariates, label=trait_name)
        traits[trait_name] = estimates.set_index("ID_std")["estimate"].reindex(accession_index)
        traits[f"{trait_name}_LSM_SE"] = estimates.set_index("ID_std")["se"].reindex(accession_index)
        traits[f"{trait_name}_N_PLOTS"] = estimates.set_index("ID_std")["n_plots"].reindex(accession_index)
        model_rows.append(diagnostics)

    simple = canonical.groupby("ID_std").agg(
        INC_SIMPLE_RAW=("INC", "mean"),
        SEV_SIMPLE_RAW=("SEV", "mean"),
        N_PLOTS_TOTAL=("INC", "size"),
        N_LOCATIONS=("Location", "nunique"),
    )
    traits = traits.join(simple, how="left")

    raw_trait_cols = [
        "INC_ADJ_RAW",
        "SEV_ADJ_RAW",
        "SEV_POS_ADJ_RAW",
        "SEV_COND_ADJ_RAW",
        "INC_SIMPLE_RAW",
        "SEV_SIMPLE_RAW",
    ]
    for col in raw_trait_cols:
        traits[col.replace("_RAW", "_INT")] = rank_int(traits[col])

    # Reproduce the submitted V6 phenotype-ingest behavior for forensic comparison:
    # rank-transform all plot rows, then retain only the first row per accession.
    forensic = canonical[["ID_std", "INC", "SEV", "Location", "Stratum"]].copy()
    forensic["SUBMITTED_V6_INC_INT"] = rank_int(forensic["INC"])
    forensic["SUBMITTED_V6_SEV_INT"] = rank_int(forensic["SEV"])
    forensic = forensic.drop_duplicates("ID_std", keep="first").set_index("ID_std")
    traits = traits.join(forensic[["SUBMITTED_V6_INC_INT", "SUBMITTED_V6_SEV_INT", "Location", "Stratum"]], how="left")
    traits = traits.rename(columns={"Location": "SUBMITTED_V6_SOURCE_LOCATION", "Stratum": "SUBMITTED_V6_SOURCE_STRATUM"})

    traits.index.name = "ID_std"
    traits = traits.reset_index()
    model_diagnostics = pd.DataFrame(model_rows)

    intermediate = Path(intermediate_dir)
    intermediate.mkdir(parents=True, exist_ok=True)
    traits.to_csv(intermediate / "phenotype_accession_traits.csv", index=False)
    model_diagnostics.to_csv(intermediate / "phenotype_model_diagnostics.csv", index=False)
    logger.info("Accession traits built for %s phenotype IDs", len(traits))
    return traits, model_diagnostics


def build_population_table(
    traits: pd.DataFrame,
    matched_ids: list[str],
    min_locations: int,
    expected_primary_n: int | None,
    intermediate_dir: str | Path,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    indexed = traits.set_index("ID_std")
    rows: list[dict[str, Any]] = []
    for sid in matched_ids:
        n_locations = int(indexed.loc[sid, "N_LOCATIONS"]) if sid in indexed.index else 0
        n_plots = int(indexed.loc[sid, "N_PLOTS_TOTAL"]) if sid in indexed.index else 0
        rows.append(
            {
                "ID_std": sid,
                "N_LOCATIONS": n_locations,
                "N_PLOTS_TOTAL": n_plots,
                "PRIMARY_ELIGIBLE": bool(n_locations >= min_locations),
                "EXCLUSION_REASON": "" if n_locations >= min_locations else f"fewer than {min_locations} locations",
            }
        )
    table = pd.DataFrame(rows)
    primary_ids = table.loc[table["PRIMARY_ELIGIBLE"], "ID_std"].tolist()
    full_ids = table["ID_std"].tolist()
    if expected_primary_n is not None and len(primary_ids) != int(expected_primary_n):
        raise ValueError(
            f"Primary population size was {len(primary_ids)}, expected {expected_primary_n}. "
            "Inspect population_definition.csv before continuing."
        )
    table.to_csv(Path(intermediate_dir) / "population_definition.csv", index=False)
    logger.info(
        "Analysis populations: primary n=%d (>= %d locations), full sensitivity n=%d",
        len(primary_ids), min_locations, len(full_ids),
    )
    excluded = table.loc[~table["PRIMARY_ELIGIBLE"], ["ID_std", "N_LOCATIONS", "N_PLOTS_TOTAL"]]
    if not excluded.empty:
        logger.info("Primary population exclusions: %s", excluded.to_dict(orient="records"))
    return table, primary_ids, full_ids


def build_leave_one_location_out_traits(
    canonical: pd.DataFrame,
    locations: list[str],
    intermediate_dir: str | Path,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    accession_index = sorted(canonical["ID_std"].unique())
    traits = pd.DataFrame(index=accession_index)
    diagnostics: list[dict[str, Any]] = []

    for omitted in locations:
        subset = canonical.loc[canonical["Location"].astype(str) != str(omitted)].copy()
        label = re.sub(r"[^A-Za-z0-9]+", "_", str(omitted)).strip("_").upper()
        for trait_stub, value_col, covariates in [
            ("INC", "INC", []),
            ("SEV", "SEV", []),
            ("SEV_COND", "SEV", ["INC"]),
        ]:
            name = f"{trait_stub}_LOO_{label}_RAW"
            estimates, diag = _fit_adjusted_means(subset, value_col, covariates, label=name)
            series = estimates.set_index("ID_std")["estimate"].reindex(accession_index)
            traits[name] = series
            traits[name.replace("_RAW", "_INT")] = rank_int(series)
            diag["omitted_location"] = omitted
            diagnostics.append(diag)

    traits.index.name = "ID_std"
    traits = traits.reset_index()
    diag_df = pd.DataFrame(diagnostics)
    intermediate = Path(intermediate_dir)
    traits.to_csv(intermediate / "phenotype_leave_one_location_out_traits.csv", index=False)
    diag_df.to_csv(intermediate / "phenotype_leave_one_location_out_diagnostics.csv", index=False)
    logger.info("Leave-one-location-out traits built for %d locations", len(locations))
    return traits, diag_df


def build_location_specific_traits(
    canonical: pd.DataFrame,
    locations: list[str],
    intermediate_dir: str | Path,
) -> pd.DataFrame:
    accession_index = sorted(canonical["ID_std"].unique())
    out = pd.DataFrame(index=accession_index)
    for location in locations:
        label = re.sub(r"[^A-Za-z0-9]+", "_", str(location)).strip("_").upper()
        sub = canonical.loc[canonical["Location"].astype(str) == str(location)]
        grouped = sub.groupby("ID_std").agg(INC=("INC", "mean"), SEV=("SEV", "mean"), N_PLOTS=("INC", "size"))
        for trait in ["INC", "SEV"]:
            raw_name = f"{trait}_SITE_{label}_RAW"
            series = grouped[trait].reindex(accession_index)
            out[raw_name] = series
            out[raw_name.replace("_RAW", "_INT")] = rank_int(series)
        out[f"N_PLOTS_SITE_{label}"] = grouped["N_PLOTS"].reindex(accession_index)
    out.index.name = "ID_std"
    out = out.reset_index()
    out.to_csv(Path(intermediate_dir) / "phenotype_location_specific_traits.csv", index=False)
    return out
