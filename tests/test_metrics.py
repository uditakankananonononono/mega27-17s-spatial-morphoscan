import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from src.morphoscan import metrics as M


def test_pearson_perfect():
    assert abs(M.pearson_r([1, 2, 3, 4], [2, 4, 6, 8]) - 1.0) < 1e-9

def test_pearson_anticorr():
    assert abs(M.pearson_r([1, 2, 3, 4], [8, 6, 4, 2]) + 1.0) < 1e-9

def test_pearson_constant_zero():
    assert M.pearson_r([1, 1, 1], [1, 2, 3]) == 0.0

def test_spearman_monotone():
    assert M.spearman_r([1, 2, 3, 4], [10, 20, 30, 40]) > 0.999

def test_ridge_recovers_linear():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(200, 5))
    w = np.array([1.0, -2.0, 0.5, 0.0, 3.0])
    y = X @ w + 0.7
    beta = M.ridge_fit(X, y, 1e-8)
    assert np.allclose(beta[:-1], w, atol=1e-3)
    assert abs(beta[-1] - 0.7) < 1e-3

def test_ridge_shrinks():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(100, 3)); y = X @ np.array([1.0, 1.0, 1.0])
    b_small = M.ridge_fit(X, y, 0.01); b_big = M.ridge_fit(X, y, 1e4)
    assert np.linalg.norm(b_big[:-1]) < np.linalg.norm(b_small[:-1])

def test_kernel_rows_sum_to_one():
    d = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.0], [2.0, 1.0, 0.0]])
    W = M.gaussian_kernel_weights(d, 1.0)
    assert np.allclose(W.sum(axis=1), 1.0)
    assert np.all(np.diag(W) == 0)

def test_smooth_alpha_zero_identity():
    y = np.arange(5.0)
    W = np.ones((5, 5)) / 4.0; np.fill_diagonal(W, 0)
    assert np.allclose(M.spatial_smooth(y, W, 0.0), y)

def test_morans_known_value():
    # alternating pattern on a 2-neighbor ring should be strongly negative
    x = np.array([0.0, 1.0, 0.0, 1.0])
    W = np.array([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]], float)
    i = M.morans_i(x, W)
    assert i < -0.9

def test_morans_clustered_positive():
    x = np.array([0.0, 0.0, 1.0, 1.0])
    W = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], float)
    assert M.morans_i(x, W) > 0.9

def test_auc_perfect():
    assert M.auc_mann_whitney([0.9, 0.8, 0.1, 0.2], [1, 1, 0, 0]) == 1.0

def test_auc_chance():
    a = M.auc_mann_whitney([0.5, 0.5, 0.5, 0.5], [1, 1, 0, 0])
    assert a == 0.5

def test_auc_inverse():
    assert M.auc_mann_whitney([0.1, 0.2, 0.8, 0.9], [1, 1, 0, 0]) == 0.0

def test_bh_monotone():
    rej, adj = M.benjamini_hochberg([0.001, 0.01, 0.5, 0.9], 0.05)
    assert adj[0] <= adj[1] <= adj[2] <= adj[3]
    assert rej[0] and not rej[3]

def test_bh_all_significant():
    rej, adj = M.benjamini_hochberg([0.001] * 10, 0.05)
    assert all(rej)

def test_permutation_p_bounds():
    p = M.permutation_pvalue(2.0, np.random.default_rng(0).normal(size=999))
    assert 0.0 < p <= 1.0

def test_wilcoxon_identical():
    assert M.wilcoxon_p([1, 2, 3], [1, 2, 3]) == 1.0

def test_brier_perfect():
    assert M.brier_score([1.0, 0.0], [1, 0]) == 0.0

def test_brier_worst():
    assert M.brier_score([0.0, 1.0], [1, 0]) == 1.0

def test_softmax_sums_one():
    p = M.softmax(np.array([1.0, 2.0, 3.0]))
    assert abs(p.sum() - 1.0) < 1e-12
    assert p[2] > p[1] > p[0]

def test_softmax_temperature():
    p_hot = M.softmax(np.array([1.0, 2.0]), tau=0.1)
    p_cold = M.softmax(np.array([1.0, 2.0]), tau=10.0)
    assert p_hot[1] > p_cold[1]
