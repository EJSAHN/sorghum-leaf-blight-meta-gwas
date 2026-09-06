from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy import stats
from lb_sensitivity import imputed_block, bh_adjust, ols_p_block, ols_projection
from lb_power import minimum_detectable_beta, two_sided_power
import pytest

def test_power_reference():
    beta=minimum_detectable_beta(100,.05,.05/567758,.80)
    assert abs(beta-2.1663796269322266)<1e-8
    assert abs(two_sided_power(beta,100,.05,.05/567758)-.80)<1e-9

def test_bh():
    np.testing.assert_allclose(bh_adjust(np.array([.01,.04,.03,.8])),[.04,.0533333333333,.0533333333333,.8])

def test_impute_before_subset():
    x=imputed_block(np.array([[0,0,2,-1],[0,1,2,-1]],dtype=np.int8))
    np.testing.assert_allclose(x[:,3],[2/3,1])

@pytest.mark.parametrize('k',[0,1,2,3,5])
def test_ols_against_direct_design(k):
    rng=np.random.default_rng(732)
    pc=rng.normal(size=(100,5));x=rng.binomial(2,.25,size=(100,7)).astype(float)
    y=rng.normal(size=(100,4))+pc[:,[0]]*.3
    m,df=ols_projection(pc,k);p=ols_p_block(x,y,m,df)
    c=np.column_stack((np.ones(100),pc[:,:k]))
    for j in range(x.shape[1]):
        design=np.column_stack((c,x[:,j]));coef=np.linalg.lstsq(design,y,rcond=None)[0]
        rss=np.sum((y-design@coef)**2,axis=0)
        se=np.sqrt(rss/df*np.linalg.inv(design.T@design)[-1,-1])
        direct=2*stats.t.sf(np.abs(coef[-1]/se),df)
        np.testing.assert_allclose(p[j],direct,rtol=1e-10,atol=1e-12)

def test_constant_marker_has_p_one():
    rng=np.random.default_rng(733);m,df=ols_projection(rng.normal(size=(100,5)),3)
    np.testing.assert_array_equal(ols_p_block(np.ones((100,2)),rng.normal(size=(100,4)),m,df),np.ones((2,4)))

def test_invalid_power_inputs():
    with pytest.raises(ValueError):two_sided_power(1,2,.5,.05)
