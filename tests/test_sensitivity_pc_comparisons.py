from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pytest
from lb_sensitivity import pc_geometry, reconstruct_pcs, imputed_block


def test_sign_and_scale_invariance():
    r=np.random.default_rng(291);a=r.normal(size=(100,5));b=a*np.array([-2,5,-1,3,.2])+12
    columns,spaces=pc_geometry(a,b)
    assert columns.column_correlation_check_passed.all()
    assert spaces.equivalent_within_numerical_tolerance.all()


def test_swapped_axes_fail_pc1_but_not_joint_pc2():
    r=np.random.default_rng(292);a=np.linalg.qr(r.normal(size=(100,5)))[0];b=a[:,[1,0,2,3,4]]
    columns,spaces=pc_geometry(a,b)
    assert not columns.column_correlation_check_passed.iloc[0]
    assert not spaces.equivalent_within_numerical_tolerance.iloc[0]
    assert spaces.equivalent_within_numerical_tolerance.iloc[1:].all()


def test_real_basis_change_is_not_bypassed():
    r=np.random.default_rng(293);a=r.normal(size=(100,5));b=a.copy();b[:,0]+=r.normal(size=100)*.1
    columns,spaces=pc_geometry(a,b)
    assert not columns.column_correlation_check_passed.iloc[0]
    assert not spaces.equivalent_within_numerical_tolerance.all()


def test_diagnostic_no_finite_pcs():
    a=np.random.default_rng(2).normal(size=(100,5));b=a.copy();b[3,1]=np.nan
    with pytest.raises(ValueError):pc_geometry(a,b)


def test_reconstructed_pcs_match_direct_svd():
    r=np.random.default_rng(294);g=r.binomial(2,.3,size=(650,102)).astype(np.int8);g[::4,1]=-1
    ids=np.arange(100);markers=np.arange(600)
    pc=reconstruct_pcs(g,ids,markers)
    x=imputed_block(g[markers])[:,ids];mu=x.mean(axis=1);p=mu/2
    z=(x-mu[:,None])/np.sqrt(2*p*(1-p))[:,None]
    _,sing,vh=np.linalg.svd(z,full_matrices=False)
    other=vh[:5].T*sing[:5]
    _,spaces=pc_geometry(other,pc)
    assert spaces.equivalent_within_numerical_tolerance.all()


def test_duplicate_structure_markers_rejected():
    g=np.random.default_rng(6).binomial(2,.3,size=(50,100))
    with pytest.raises(ValueError): reconstruct_pcs(g,np.arange(100),np.array([0,1,2,3,4,5,5]))
