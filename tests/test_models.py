import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from src.morphoscan import models as M, dataio as D


def test_library_normalize_rows_scale():
    c = np.array([[1.0, 3.0], [0.0, 0.0]])
    out = M.library_normalize_log1p(c)
    assert np.allclose(out[0], np.log1p(np.array([0.25, 0.75]) * 1e4))
    assert np.allclose(out[1], 0.0)

def test_library_normalize_sparse_matches_dense():
    from scipy import sparse
    rng = np.random.default_rng(0)
    c = rng.poisson(2, (20, 8)).astype(float)
    d = M.library_normalize_log1p(c)
    s = M.library_normalize_log1p(sparse.csr_matrix(c)).toarray()
    assert np.allclose(d, s, atol=1e-3)

def test_select_top_genes_order():
    x = np.array([[1.0, 5.0, 3.0], [1.0, 5.0, 3.0]])
    assert list(M.select_top_genes(x, 2)) == [1, 2]

def test_standardizer_zero_mean_unit_sd_and_constant_col():
    x = np.array([[1.0, 7.0], [3.0, 7.0], [5.0, 7.0]])
    z = M.Standardizer().fit(x).transform(x)
    assert np.allclose(z.mean(axis=0), 0) and abs(z[:, 0].std() - 1) < 1e-9
    assert np.all(z[:, 1] == 0)

def test_ridge_bank_recovers_multioutput():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(300, 4)); W = rng.normal(size=(4, 3))
    Y = X @ W + 0.5
    b = M.fit_ridge_bank(X, Y, 1e-8)
    assert np.allclose(M.predict_bank(X, b), Y, atol=1e-4)

def test_per_gene_pearson_known():
    T = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 5.0]])
    P = np.array([[2.0, 3.0], [4.0, 2.0], [6.0, 1.0]])
    r = M.per_gene_pearson(P, T)
    assert abs(r[0] - 1.0) < 1e-9 and r[1] < 0

def test_choose_lambda_inner_returns_grid_value():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(120, 3)); Y = X @ rng.normal(size=(3, 2)) + rng.normal(size=(120, 2))
    g = np.repeat(np.arange(4), 30)
    assert M.choose_lambda_inner(X, Y, g) in (0.1, 1.0, 10.0, 100.0)

def test_choose_alpha_smooths_noise():
    rng = np.random.default_rng(3)
    xs, ys = np.meshgrid(np.arange(8), np.arange(8))
    C = np.c_[xs.ravel(), ys.ravel()].astype(float)
    T = (C[:, [0]] + 0.0) @ np.ones((1, 3))
    P = T + rng.normal(0, 3, T.shape)
    a = M.choose_alpha_inner([P], [T], [C], 1.5)
    assert a > 0.0

def test_extract_patch_clamps_edges():
    img = np.arange(100).reshape(10, 10)
    assert D.extract_patch(img, 0, 0, 3).shape == (3, 3)
    assert D.extract_patch(img, 5, 5, 2).shape == (4, 4)
