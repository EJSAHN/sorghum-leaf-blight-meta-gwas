from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from lb_power import build_power_table
from lb_sensitivity import regenerate_ols
from run_sensitivity import prepare_support, run


def test_no_independent_ols_reference_is_not_claimed(tmp_path):
    ref=tmp_path/'ref';out=tmp_path/'out';ref.mkdir();out.mkdir()
    rng=np.random.default_rng(131)
    inp=dict(geno=rng.binomial(2,.3,size=(25,14)),y=rng.normal(size=(12,4)),
             pc=rng.normal(size=(12,5)),indices=np.arange(12),n_structure_markers=14)
    r=regenerate_ols(inp,ref,out)
    assert r['passed'] and r['rows']==20
    assert r['reference_comparison_performed'] is False
    assert not (out/'OLS_Reproduction.csv').exists()


def test_power_output_units():
    result=build_power_table([100],[.05],[.8],{'x':.05/567758})
    assert result.iloc[0].minimum_detectable_beta_residual_SD_per_allele==pytest.approx(2.1663796269322266)
    assert all('variance_explained' not in c for c in result)
    assert 'residual SD per allele' in result.iloc[0].method_note


def test_power_reads_old_project_schema_without_modifying_it(tmp_path):
    p=tmp_path/'project';out=tmp_path/'out';(p/'06_tables').mkdir(parents=True);out.mkdir()
    path=p/'06_tables/Power.csv'
    pd.DataFrame([dict(n_samples=100,maf=.05,alpha=1e-6,target_power=.8,
        minimum_detectable_standardized_per_allele_beta=2,
        approximate_marginal_variance_explained=.4)]).to_csv(path,index=False)
    before=path.read_bytes();support=prepare_support(p,out,'power',None)
    assert path.read_bytes()==before
    new=pd.read_csv(support/'Power.csv')
    assert 'minimum_detectable_beta_residual_SD_per_allele' in new
    assert 'approximate_marginal_variance_explained' not in new


def test_output_cannot_be_inside_input_project(tmp_path):
    p=tmp_path/'project';p.mkdir()
    with pytest.raises(ValueError,match='outside'):
        run(p,p/'newout','power')
    assert not (p/'newout').exists()
