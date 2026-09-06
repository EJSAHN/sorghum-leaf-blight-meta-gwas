"""OLS diagnostics and matched-estimator candidate effects for sorghum leaf blight."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Any
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import subspace_angles
from lb_power import minimum_detectable_beta
TRAITS = ('INC_ADJ_RAW', 'INC_ADJ_INT', 'SEV_ADJ_RAW', 'SEV_ADJ_INT')
PC_COUNTS = (0, 1, 2, 3, 5)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def first_file(root: Path, relative: tuple[str, ...]) -> Path:
    for name in relative:
        p = root / name
        if p.is_file():
            return p
    raise FileNotFoundError('Required file not found: ' + ' or '.join(str(root / s) for s in relative))


def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError('All p-values must be finite and within [0,1]')
    idx = np.argsort(p, kind='stable')
    ordered = np.minimum.accumulate((p[idx] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    q = np.empty_like(p); q[idx] = np.minimum(ordered, 1)
    return q


def imputed_block(raw: np.ndarray) -> np.ndarray:
    """Impute over every QC sample before selecting the primary sample set."""
    x = np.array(raw, dtype=np.float64, copy=True)
    if not np.isin(x, (-1, 0, 1, 2)).all():
        raise ValueError('Genotype cache must contain only -1,0,1,2')
    called = x >= 0
    counts = called.sum(axis=1)
    if np.any(counts == 0):
        raise ValueError('All-missing genotype row')
    means = np.where(called, x, 0).sum(axis=1) / counts
    a, b = np.where(~called); x[a,b] = means[a]
    return x


def ols_projection(pc: np.ndarray, n_pc: int) -> tuple[np.ndarray, int]:
    n = len(pc)
    cov = np.ones((n,1))
    if n_pc:
        selected = np.asarray(pc[:,:n_pc], dtype=float)
        scale = selected.std(axis=0, ddof=0)
        if np.any(scale <= 1e-12):
            raise ValueError('A requested PC is constant')
        selected = (selected - selected.mean(axis=0)) / scale
        cov = np.column_stack((cov, selected))
    rank = int(np.linalg.matrix_rank(cov))
    if rank != cov.shape[1]:
        raise ValueError('Covariates are rank deficient')
    q, _ = np.linalg.qr(cov, mode='reduced')
    return np.eye(n) - q @ q.T, n-rank-1


def ols_p_block(x: np.ndarray, y: np.ndarray, projection: np.ndarray, df: int) -> np.ndarray:
    """FWL-adjusted two-sided t probabilities for SNP columns and phenotype columns."""
    xr = projection @ x
    yr = projection @ y
    denom = np.sum(xr*xr, axis=0)
    numer = xr.T @ yr
    yss = np.sum(yr*yr, axis=0)
    good = denom > 1e-12
    p = np.ones((x.shape[1], y.shape[1]), dtype=float)
    beta = np.divide(numer, denom[:,None], out=np.zeros_like(numer), where=good[:,None])
    rss = np.maximum(yss[None,:] - numer*beta, 0)
    se = np.sqrt(np.divide(rss/df, denom[:,None], out=np.zeros_like(rss), where=good[:,None]))
    valid = good[:,None] & np.isfinite(se) & (se > 0)
    t = np.divide(beta, se, out=np.zeros_like(beta), where=valid)
    p[valid] = 2*stats.t.sf(np.abs(t[valid]), df)
    return np.clip(p, 1e-300, 1)


def _indexed_csv(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, dtype={"ID_std": str})
    if "ID_std" not in table or table["ID_std"].duplicated().any():
        raise ValueError("Missing or duplicate ID_std in " + str(path))
    return table.set_index("ID_std")


def pc_geometry(reference_pc: np.ndarray, observed_pc: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare columns and nested covariate spaces without changing either basis."""
    a = np.asarray(reference_pc, dtype=float)
    b = np.asarray(observed_pc, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 5:
        raise ValueError("PC matrices must be aligned n-by-5 arrays")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("PC matrices contain missing or non-finite values")
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    sa = np.sqrt(np.sum(a*a, axis=0)); sb = np.sqrt(np.sum(b*b, axis=0))
    if np.any(sa <= 1e-12) or np.any(sb <= 1e-12):
        raise ValueError("A PC is constant")
    cross = a.T @ b / (sa[:,None]*sb[None,:])
    per_pc = []
    for j in range(5):
        best = int(np.argmax(np.abs(cross[:,j])))
        per_pc.append({"PC":f"PC{j+1}", "same_column_correlation":float(cross[j,j]),
                       "same_column_absolute_correlation":float(abs(cross[j,j])),
                       "best_reference_PC":f"PC{best+1}",
                       "best_absolute_correlation":float(abs(cross[best,j])),
                       "column_correlation_check_passed":bool(abs(cross[j,j]) >= 1-1e-8)})
    spaces = []
    for k in (1,2,3,5):
        aa = a[:,:k]/sa[:k]; bb = b[:,:k]/sb[:k]
        if np.linalg.matrix_rank(aa) != k or np.linalg.matrix_rank(bb) != k:
            raise ValueError("PC covariates are rank deficient")
        qa = np.linalg.qr(aa, mode='reduced')[0]
        qb = np.linalg.qr(bb, mode='reduced')[0]
        distance = float(np.linalg.norm(qa@qa.T-qb@qb.T, ord='fro'))
        angle = float(np.max(subspace_angles(aa, bb)))
        spaces.append({"n_pc":k,"projector_frobenius_distance":distance,
                       "maximum_subspace_angle_degrees":float(np.degrees(angle)),
                       "equivalent_within_numerical_tolerance":bool(distance <= 1e-6)})
    return pd.DataFrame(per_pc), pd.DataFrame(spaces)


def reconstruct_pcs(geno: np.ndarray, indices: np.ndarray,
                    marker_indices: np.ndarray, n_pc: int=5) -> np.ndarray:
    """Reproduce the stored structure convention from specified genotype rows."""
    marker_indices = np.asarray(marker_indices, dtype=int)
    if len(marker_indices) < n_pc+1 or len(np.unique(marker_indices)) != len(marker_indices):
        raise ValueError("Invalid or duplicate structure-marker indices")
    if marker_indices.min() < 0 or marker_indices.max() >= len(geno):
        raise ValueError("Structure-marker indices are outside the genotype cache")
    # Impute on the QC sample set: impute on the full QC population,
    # then select the primary population and standardise each marker.
    x = imputed_block(geno[marker_indices])[:,indices]
    means = x.mean(axis=1)
    f = np.clip(means/2.0, 1e-6, 1-1e-6)
    scale = np.sqrt(2*f*(1-f))
    valid = np.isfinite(scale) & (scale > 1e-8) & (np.var(x,axis=1) > 1e-12)
    if valid.sum() < n_pc+1:
        raise ValueError("Too few polymorphic structure markers")
    z = (x[valid]-means[valid,None])/scale[valid,None]
    k = z.T@z/len(z)
    k = (k+k.T)/2
    factor = float(np.mean(np.diag(k)))
    if not np.isfinite(factor) or factor <= 1e-12:
        raise ValueError("Invalid structure covariance scale")
    k /= factor
    values, vectors = np.linalg.eigh((k+k.T)/2)
    order = np.argsort(values)[::-1]
    return vectors[:,order[:n_pc]]*np.sqrt(np.maximum(values[order[:n_pc]],0))[None,:]


def unix_path(p: Path) -> str:
    s=str(p.resolve()).replace('\\','/')
    if os.name=='nt':
        if not re.match(r'^[A-Za-z]:/',s): raise ValueError('WSL requires a drive-letter path')
        return '/mnt/'+s[0].lower()+s[2:]
    return s


def candidate_score_effects(project: Path, reference: Path, output: Path, exe: str|None) -> dict[str, Any]:
    base=project/'04_intermediate'
    plink=base/'plink_primary_n100'
    mapping_path=first_file(plink,('variant_id_mapping.csv.gz',))
    mapping=pd.read_csv(mapping_path)
    candidates=pd.read_csv(reference/'Top_Ranked_Candidates.csv')
    cols=['variant_index','gemma_id','snp','CHR','pos']
    mapping=mapping[cols]
    cand=candidates.merge(mapping,left_on='lead_snp',right_on='snp',suffixes=('_display',''),validate='one_to_one')
    if len(cand)!=len(candidates):raise ValueError('Incomplete candidate mapping')
    for field in ('CHR','pos'):
        if not np.array_equal(cand[field+'_display'].to_numpy(int),cand[field].to_numpy(int)):
            raise ValueError('Candidate genomic position disagrees with mapping: '+field)
    if cand['gemma_id'].duplicated().any(): raise ValueError('Duplicate mapped candidate IDs')
    for trait in ('inc','sev'):
        for ext in ('.bed','.bim','.fam'):
            if not (plink/f'leaf_blight_{trait}{ext}').is_file():raise FileNotFoundError('Missing PLINK input for '+trait)
    binary=Path(exe) if exe else project/'01_tools'/'gemma-0.98.5'
    if not binary.is_file():raise FileNotFoundError('GEMMA executable not found: '+str(binary))
    local_loo=first_file(project/'06_tables',('Final_Leave_One_Location_Out_Candidates.csv','Leave_One_Location_Out_Candidates.csv'))
    loo=pd.read_csv(local_loo)
    # Confirm that the inputs correspond to the reference score-probability summaries.
    expected_population = pd.read_csv(reference/'Population_Definition.csv')
    primary_mask = expected_population['PRIMARY_ELIGIBLE'].astype(str).str.lower().eq('true')
    primary_ids = expected_population.loc[primary_mask,'ID_std'].astype(str).tolist()
    for trait in ('inc', 'sev'):
        fam = pd.read_csv(plink/f'leaf_blight_{trait}.fam', sep=r'\s+', header=None, dtype={0:str,1:str})
        if fam.iloc[:,1].astype(str).tolist() != primary_ids:
            raise ValueError('PLINK individual order differs from the reference: '+trait)
    reference_loo = pd.read_csv(reference/'LOO_Candidates.csv')
    key = ['lead_snp','omitted_location']
    if loo.duplicated(key).any() or reference_loo.duplicated(key).any():
        raise ValueError('Duplicate candidate/omission records')
    aligned = loo.merge(reference_loo, on=key, suffixes=('_local','_reference'), validate='one_to_one')
    if len(aligned) != len(reference_loo):
        raise ValueError('Candidate/omission records do not match the reference')
    for field in ('p_INC_score','p_SEV_score','p_ACAT_score'):
        if not np.allclose(aligned[field+'_local'],aligned[field+'_reference'],rtol=1e-9,atol=1e-12):
            raise ValueError('Omission probabilities differ from the reference table: '+field)
    if 'beta_INC_score' not in loo or 'beta_SEV_score' not in loo:
        raise ValueError('Local omission table has no saved mode-3 coefficients')
    for trait in ('inc','sev'):
        fam = pd.read_csv(plink/f'leaf_blight_{trait}.fam',sep=r'\s+',header=None,dtype={0:str,1:str})
        expected_trait = _indexed_csv(reference/'Accession_Phenotypes.csv').reindex(primary_ids)[trait.upper()+'_ADJ_INT']
        if not np.allclose(fam.iloc[:,5].to_numpy(float),expected_trait.to_numpy(float),rtol=1e-9,atol=1e-10):
            raise ValueError('PLINK phenotype values differ from the reference: '+trait)
    primary={};commands=[]
    for chrom in sorted(cand['CHR'].unique()):
        selected=cand.loc[cand['CHR']==chrom]
        for trait in ('INC','SEV'):
            run=output/'gemma_candidate_score'/f'{trait.lower()}_chr{int(chrom)}';run.mkdir(parents=True)
            snps=run/'selected_snps.txt'
            snps.write_text('\n'.join(selected['gemma_id'])+'\n',encoding='ascii',newline='\n')
            prefix=plink/f'leaf_blight_{trait.lower()}'
            kinship=base/'kinship_primary_n100'/f'loco_chr{int(chrom)}.cXX.txt'
            if not kinship.is_file():raise FileNotFoundError(str(kinship))
            k = np.loadtxt(kinship)
            if k.shape != (len(primary_ids),len(primary_ids)) or not np.isfinite(k).all() or not np.allclose(k,k.T,atol=1e-7):
                raise ValueError('Invalid supplied LOCO kinship matrix: '+str(kinship))
            label=f'candidate_{trait.lower()}_chr{int(chrom)}'
            args=[unix_path(binary),'-bfile',unix_path(prefix),'-k',unix_path(kinship),'-lmm','3','-maf','0','-miss','1','-snps',unix_path(snps),'-o',label]
            command='cd '+shlex.quote(unix_path(run))+' && '+' '.join(map(shlex.quote,args))
            actual=['wsl.exe','-e','sh','-lc',command] if os.name=='nt' else ['sh','-lc',command]
            commands.append(actual)
            (output/'GEMMA_commands.json').write_text(json.dumps(commands,indent=2),encoding='utf-8')
            print(f'GEMMA score-only candidate coefficients: {trait}, chr{int(chrom)}',flush=True)
            with (run/'command.log').open('w',encoding='utf-8') as log:
                proc=subprocess.run(actual,stdout=log,stderr=subprocess.STDOUT,check=False)
            if proc.returncode:raise RuntimeError('GEMMA failed; see '+str(run/'command.log'))
            assoc=run/'output'/f'{label}.assoc.txt'
            if not assoc.is_file():raise FileNotFoundError(str(assoc))
            a=pd.read_csv(assoc,sep=r'\s+')
            if a['rs'].duplicated().any() or set(a['rs'])!=set(selected['gemma_id']):
                raise ValueError('Candidate SNP selection did not match GEMMA output')
            for row in a.to_dict('records'):
                name=selected.loc[selected['gemma_id']==row['rs'],'snp'].iloc[0]
                expected=float(candidates.loc[candidates['lead_snp']==name,f'p_{trait}_score'].iloc[0])
                if not np.isclose(row['p_score'],expected,rtol=3e-6,atol=1e-12):
                    raise ValueError(f'Score p-value changed for {trait} {name}: {row["p_score"]} versus {expected}')
                primary[(name,trait)]=row
    allrows=[]
    for r in loo.to_dict('records'):
        if r['lead_snp'] not in set(candidates['lead_snp']):continue
        for trait in ('INC','SEV'):
            a=primary[(r['lead_snp'],trait)]
            b=float(r[f'beta_{trait}_score'])
            if not np.isfinite([a['beta'],a['se'],b]).all():raise ValueError('Non-finite coefficient')
            allrows.append(dict(lead_snp=r['lead_snp'],trait=trait,omitted_location=r['omitted_location'],primary_beta_null_ML=float(a['beta']),primary_se_null_ML=float(a['se']),omission_beta_null_ML=b,primary_p_score=float(a['p_score']),omission_p_score=float(r[f'p_{trait}_score']),direction_preserved=bool(np.sign(a['beta'])==np.sign(b)),effect_unit='Inverse-normal phenotype units per synthetic allele1 A',estimator='GEMMA mode 3; covariance ratio estimated under the null by ML'))
    if len(allrows)!=len(candidates)*5*2:raise ValueError('Unexpected candidate-by-omission count')
    df=pd.DataFrame(allrows);df.to_csv(output/'Candidate_Effect_Direction_Check.csv',index=False)
    pd.DataFrame([dict(trait=t,n_comparisons=int((df.trait==t).sum()),direction_preserved=int(df.loc[df.trait==t,'direction_preserved'].sum())) for t in ('INC','SEV')]).to_csv(output/'Candidate_Direction_Summary.csv',index=False)
    (output/'GEMMA_commands.json').write_text(json.dumps(commands,indent=2),encoding='utf-8')
    return {'passed':True,'candidate_trait_omission_rows':len(df),'source_mapping_sha256':sha256(mapping_path),'source_omission_sha256':sha256(local_loo),'gemma_binary_sha256':sha256(binary),
            'effect_estimator':'GEMMA mode 3 on both sides',
            'local_input_sha256':{str(q):sha256(q) for q in sorted(set(
                [plink/f'leaf_blight_{t}{ext}' for t in ('inc','sev') for ext in ('.bed','.bim','.fam')]+
                [base/'kinship_primary_n100'/f'loco_chr{int(c)}.cXX.txt' for c in cand.CHR.unique()]
            ))}}


def read_primary_inputs(project: Path, reference: Path, output: Path,
                        expected_dimensions=(100, 567758, 102)) -> dict[str, Any]:
    """Validate sample order, phenotype values, markers and the primary PC space."""
    inter = project / '04_intermediate'
    manifest_path = first_file(inter, ('genotype_cache/cache_manifest.json',))
    binary = first_file(inter, ('genotype_cache/genotypes_int8.bin',))
    traits_path = first_file(inter, ('phenotype_all_accession_traits.csv',
        'phenotype_all_accession_traits_v2.csv', 'phenotype_all_accession_traits_v02.csv'))
    pc_path = first_file(inter, ('structure_cache/PRIMARY_N100/pcs.csv',))
    marker_path = first_file(inter, ('structure_cache/common_ld_pruned_structure_markers.csv',))
    population_path = first_file(inter, ('population_definition.csv',))
    metadata_path = first_file(inter, ('genotype_cache/variant_metadata.csv.gz',))
    meta = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    all_ids = list(map(str, meta['sample_ids']))
    population = pd.read_csv(population_path, dtype={'ID_std': str})
    mask = population['PRIMARY_ELIGIBLE'].astype(str).str.lower().eq('true')
    ids = population.loc[mask, 'ID_std'].astype(str).tolist()
    expected = pd.read_csv(reference/'Population_Definition.csv', dtype={'ID_std': str})
    emask = expected['PRIMARY_ELIGIBLE'].astype(str).str.lower().eq('true')
    if ids != expected.loc[emask,'ID_std'].astype(str).tolist():
        raise ValueError('Primary sample order differs from the reference')
    if len(set(all_ids)) != len(all_ids) or len(set(ids)) != len(ids):
        raise ValueError('Duplicate sample identifiers')
    if (len(ids), int(meta['n_variants']), len(all_ids)) != tuple(expected_dimensions):
        raise ValueError('Unexpected sample or marker dimensions')
    if binary.stat().st_size != int(meta['n_variants'])*len(all_ids):
        raise ValueError('Genotype binary length disagrees with its manifest')
    traits = _indexed_csv(traits_path).reindex(ids)
    expected_traits = _indexed_csv(reference/'Accession_Phenotypes.csv').reindex(ids)
    for name in TRAITS:
        if not np.allclose(traits[name], expected_traits[name], atol=1e-9, rtol=1e-10):
            raise ValueError('Phenotype does not match reference: '+name)
    pc_cols = [f'PC{k}' for k in range(1,6)]
    pc = _indexed_csv(pc_path).reindex(ids)[pc_cols].to_numpy(float)
    refpc = _indexed_csv(reference/'Structure_PCs.csv').reindex(ids)[pc_cols].to_numpy(float)
    cols, spaces = pc_geometry(refpc, pc)
    if not spaces.equivalent_within_numerical_tolerance.all():
        raise ValueError('Stored PCs disagree with the primary-project reference basis')
    markers = pd.read_csv(marker_path)
    refmarkers = pd.read_csv(reference/'Structure_Markers.csv')
    fields = ['variant_index','snp','CHR','pos']
    if not markers[fields].astype(str).reset_index(drop=True).equals(refmarkers[fields].astype(str).reset_index(drop=True)):
        raise ValueError('Structure-marker list differs from the reference')
    metadata = pd.read_csv(metadata_path, usecols=fields)
    mi = markers.variant_index.to_numpy(int)
    if len(metadata) != int(meta['n_variants']) or not np.array_equal(metadata.variant_index, np.arange(len(metadata))):
        raise ValueError('Genotype metadata order is invalid')
    if np.any(mi<0) or np.any(mi>=len(metadata)):
        raise ValueError('Structure-marker index is out of range')
    if not metadata.iloc[mi][fields].astype(str).reset_index(drop=True).equals(markers[fields].astype(str).reset_index(drop=True)):
        raise ValueError('Structure markers disagree with genotype metadata')
    indices = np.array([all_ids.index(s) for s in ids], dtype=int)
    geno = np.memmap(binary, dtype=np.int8, mode='r', shape=(int(meta['n_variants']), len(all_ids)))
    try:
        rebuilt = reconstruct_pcs(geno, indices, mi)
        columns, reconstructed_spaces = pc_geometry(pc, rebuilt)
        columns.to_csv(output/'PC_Reconstruction_Columns.csv', index=False)
        reconstructed_spaces.to_csv(output/'PC_Reconstruction_Subspaces.csv', index=False)
        if not reconstructed_spaces.equivalent_within_numerical_tolerance.all():
            raise ValueError('PC reconstruction from genotype rows failed')
        hashes={str(p.relative_to(project)):sha256(p) for p in
            [manifest_path,binary,traits_path,pc_path,marker_path,population_path,metadata_path]}
        return dict(geno=geno, y=traits[list(TRAITS)].to_numpy(float), indices=indices,
                    pc=pc, ids=ids, hashes=hashes, n_structure_markers=len(mi))
    except BaseException:
        geno._mmap.close()
        raise


def regenerate_ols(inputs: dict[str, Any], reference: Path, output: Path) -> dict[str, Any]:
    """Reproduce raw/INT diagnostic scans with the primary genomic PC basis."""
    geno, y, pc = inputs['geno'], inputs['y'], inputs['pc']
    projections={k:ols_projection(pc,k) for k in PC_COUNTS}
    pv={k:np.empty((len(geno),len(TRAITS)),dtype=float) for k in PC_COUNTS}
    for start in range(0,len(geno),20000):
        end=min(start+20000,len(geno))
        x=imputed_block(geno[start:end])[:,inputs['indices']].T
        for k,(proj,df) in projections.items():
            pv[k][start:end]=ols_p_block(x,y,proj,df)
        print(f'OLS markers: {end:,}/{len(geno):,}',flush=True)
    print('Calculating FDR and distribution summaries...',flush=True)
    rows=[]
    for j,trait in enumerate(TRAITS):
        for k in PC_COUNTS:
            p=pv[k][:,j]
            rows.append(dict(model=f'OLS_PC{k}',trait=trait,
                phenotype_treatment='INT' if trait.endswith('INT') else 'raw',
                n_pc=k,n_samples=len(y),
                lambda_gc=float(np.median(stats.chi2.isf(p,1))/stats.chi2.ppf(.5,1)),
                minimum_p=float(p.min()),fdr_significant_snps=int((bh_adjust(p)<.10).sum())))
    table=pd.DataFrame(rows)
    table.to_csv(output/'Transformation_Checks.csv',index=False)
    if not (reference/'Transformation_Checks.csv').is_file():
        return {'passed':True, 'rows':len(rows),
                'n_structure_markers':inputs.get('n_structure_markers'),
                'reference_comparison_performed':False}
    expected=pd.read_csv(reference/'Transformation_Checks.csv').set_index(['trait','n_pc'])
    comparison=[]
    for row in rows:
        ref=expected.loc[(row['trait'],row['n_pc'])]
        for metric,tolerance,scale in [('lambda_gc',2e-6,'absolute'),('minimum_p',2e-4,'log10'),('fdr_significant_snps',0,'exact')]:
            observed=float(row[metric]);target=float(ref[metric])
            delta=abs(np.log10(observed)-np.log10(target)) if scale=='log10' else abs(observed-target)
            comparison.append(dict(trait=row['trait'],n_pc=row['n_pc'],metric=metric,
                observed=observed,expected=target,difference=delta,tolerance=tolerance,
                comparison_scale=scale,passed=bool(delta<=tolerance)))
    comp=pd.DataFrame(comparison);comp.to_csv(output/'OLS_Reproduction.csv',index=False)
    return {'passed':bool(comp.passed.all()),'rows':len(rows),'n_structure_markers':inputs.get('n_structure_markers'), 'reference_comparison_performed':True}


def reproduce_power(reference: Path, output: Path) -> dict[str, Any]:
    table=pd.read_csv(reference/'Power.csv');errors=[];rows=[]
    for row in table.to_dict('records'):
        value=minimum_detectable_beta(int(row['n_samples']),float(row['maf']),float(row['alpha']),float(row['target_power']))
        errors.append(abs(value-float(row['minimum_detectable_beta_residual_SD_per_allele'])))
        row['minimum_detectable_beta_residual_SD_per_allele']=value;rows.append(row)
    pd.DataFrame(rows).to_csv(output/'Power.csv',index=False)
    return {'passed':bool(errors and max(errors)<1e-8),'rows':len(rows),'maximum_absolute_difference':max(errors) if errors else None}
