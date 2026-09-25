"""Metrics and statistics for MorphoScan. Every function is unit-tested against
known values in tests/ (see tests/test_metrics.py)."""
from __future__ import annotations
import numpy as np
from scipy import stats as _st


def pearson_r(x, y):
    """Pearson product-moment correlation r = cov(x,y)/(sx*sy)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.size < 3 or np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def spearman_r(x, y):
    """Spearman rank correlation (Pearson on ranks)."""
    r, _ = _st.spearmanr(x, y)
    return float(r) if np.isfinite(r) else 0.0


def ridge_fit(X, y, lam):
    """Closed-form ridge: beta = (X^T X + lam I)^-1 X^T y, intercept unpenalized."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n = X.shape[0]
    Xa = np.hstack([X, np.ones((n, 1))])
    P = np.eye(Xa.shape[1]) * lam; P[-1, -1] = 0.0
    beta = np.linalg.solve(Xa.T @ Xa + P, Xa.T @ y)
    return beta


def ridge_predict(X, beta):
    X = np.asarray(X, float)
    Xa = np.hstack([X, np.ones((X.shape[0], 1))])
    return Xa @ beta


def gaussian_kernel_weights(d, length_scale):
    """w_ij = exp(-d_ij^2 / (2 l^2)); rows normalized to sum to 1."""
    d = np.asarray(d, float)
    W = np.exp(-(d ** 2) / (2.0 * length_scale ** 2))
    np.fill_diagonal(W, 0.0)
    s = W.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return W / s


def spatial_smooth(y_hat, W, alpha):
    """Two-stage predictor: y2 = (1-alpha) y_hat + alpha W y_hat."""
    y_hat = np.asarray(y_hat, float)
    return (1.0 - alpha) * y_hat + alpha * (W @ y_hat)


def morans_i(x, W):
    """Moran's I = (n/S0) * sum_ij w_ij (x_i-xbar)(x_j-xbar) / sum_i (x_i-xbar)^2."""
    x = np.asarray(x, float)
    n = x.size
    z = x - x.mean()
    denom = (z ** 2).sum()
    if denom == 0:
        return 0.0
    S0 = W.sum()
    if S0 == 0:
        return 0.0
    num = float((W * np.outer(z, z)).sum())
    return (n / S0) * num / denom


def auc_mann_whitney(scores, labels):
    """ROC AUC = P(score_pos > score_neg) via the Mann-Whitney U statistic."""
    scores = np.asarray(scores, float); labels = np.asarray(labels)
    pos = scores[labels == 1]; neg = scores[labels == 0]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    U, _ = _st.mannwhitneyu(pos, neg, alternative="greater")
    return float(U / (pos.size * neg.size))


def benjamini_hochberg(pvals, alpha=0.05):
    """BH-FDR: return boolean rejected array and adjusted p-values."""
    p = np.asarray(pvals, float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(m); out[order] = adj
    return out <= alpha, out


def permutation_pvalue(stat_obs, stat_null):
    """Two-sided permutation p = (1 + #|null| >= |obs|) / (1 + n_null)."""
    stat_null = np.asarray(stat_null, float)
    return float((1 + np.sum(np.abs(stat_null) >= abs(stat_obs))) / (1 + stat_null.size))


def wilcoxon_p(x, y):
    """Paired Wilcoxon signed-rank p-value (two-sided)."""
    d = np.asarray(x, float) - np.asarray(y, float)
    d = d[d != 0]
    if d.size < 5:
        return 1.0
    return float(_st.wilcoxon(d, alternative="two-sided").pvalue)


def brier_score(probs, labels):
    """Brier score = mean((p - y)^2)."""
    probs = np.asarray(probs, float); labels = np.asarray(labels, float)
    return float(np.mean((probs - labels) ** 2))


def softmax(logits, tau=1.0):
    """alpha_i = exp(s_i/tau) / sum_j exp(s_j/tau)."""
    z = np.asarray(logits, float) / tau
    z = z - z.max()
    e = np.exp(z)
    return e / e.sum()
