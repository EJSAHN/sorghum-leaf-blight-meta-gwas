from pathlib import Path
import sys, json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pandas as pd
from lb_sensitivity import read_primary_inputs, regenerate_ols, reconstruct_pcs, TRAITS, PC_COUNTS, ols_projection, ols_p_block, imputed_block, bh_adjust
from scipy import stats
import pytest

def make_project(tmp_path):
    project=tmp_path/'project'; ref=tmp_path/'reference'; out=tmp_path/'out'
    inter=project/'04_intermediate';cache=inter/'genotype_cache';struct=inter/'structure_cache'
    for p in [cache,struct/'PRIMARY_N100',ref,out]: p.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(9931);g=rng.binomial(2,.35,size=(350,32)).astype(np.int8);g[::8,5]=-1
    all_ids=[f'L{i:03}' for i in range(32)];ids=all_ids[:30];indices=np.arange(30);mi=np.arange(180)
    pc=reconstruct_pcs(g,indices,mi);y=rng.normal(size=(32,4))
    g.tofile(cache/'genotypes_int8.bin')
    (cache/'cache_manifest.json').write_text(json.dumps(dict(n_variants=350,sample_ids=all_ids)))
    pop=pd.DataFrame(dict(ID_std=all_ids,PRIMARY_ELIGIBLE=[True]*30+[False]*2))
    pop.to_csv(inter/'population_definition.csv',index=False);pop.to_csv(ref/'Population_Definition.csv',index=False)
    traits=pd.DataFrame(y,columns=TRAITS);traits.insert(0,'ID_std',all_ids)
    traits.to_csv(inter/'phenotype_all_accession_traits.csv',index=False);traits.to_csv(ref/'Accession_Phenotypes.csv',index=False)
    pcs=pd.DataFrame(pc,columns=[f'PC{i}' for i in range(1,6)]);pcs.insert(0,'ID_std',ids)
    pcs.to_csv(struct/'PRIMARY_N100/pcs.csv',index=False);pcs.to_csv(ref/'Structure_PCs.csv',index=False)
    md=pd.DataFrame(dict(variant_index=np.arange(350),snp=[f'S01_{i+1}' for i in range(350)],CHR=1,pos=np.arange(1,351)))
    md.to_csv(cache/'variant_metadata.csv.gz',index=False);md.iloc[mi].to_csv(struct/'common_ld_pruned_structure_markers.csv',index=False);md.iloc[mi].to_csv(ref/'Structure_Markers.csv',index=False)
    x=imputed_block(g)[:,indices].T;data=[]
    for j,t in enumerate(TRAITS):
        for k in PC_COUNTS:
            m,df=ols_projection(pc,k);p=ols_p_block(x,y[:30],m,df)[:,j]
            data.append(dict(trait=t,n_pc=k,lambda_gc=float(np.median(stats.chi2.isf(p,1))/stats.chi2.ppf(.5,1)),minimum_p=float(p.min()),fdr_significant_snps=int((bh_adjust(p)<.1).sum())))
    pd.DataFrame(data).to_csv(ref/'Transformation_Checks.csv',index=False)
    return project,ref,out

def test_primary_basis_workflow(tmp_path):
    p,r,o=make_project(tmp_path); inp=read_primary_inputs(p,r,o,(30,350,32))
    try:
        result=regenerate_ols(inp,r,o)
        assert result['passed'] and result['rows']==20 and result['n_structure_markers']==180
    finally:inp['geno']._mmap.close()

def test_primary_basis_mismatch_stops(tmp_path):
    p,r,o=make_project(tmp_path);path=p/'04_intermediate/structure_cache/PRIMARY_N100/pcs.csv'
    table=pd.read_csv(path);table.loc[0,'PC1']+=.2;table.to_csv(path,index=False)
    with pytest.raises(ValueError,match='PCs disagree'):read_primary_inputs(p,r,o,(30,350,32))

def test_reference_marker_mismatch_stops(tmp_path):
    p,r,o=make_project(tmp_path);path=r/'Structure_Markers.csv';md=pd.read_csv(path);md.loc[0,'pos']=99;md.to_csv(path,index=False)
    with pytest.raises(ValueError,match='marker list differs'):read_primary_inputs(p,r,o,(30,350,32))
