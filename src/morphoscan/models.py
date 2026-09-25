"""Model pipeline for MorphoScan Track A (expression) and Track B (tumor scan).
All selection (genes, scalers, lambda, alpha) happens inside training folds only."""
from __future__ import annotations
import numpy as np
from . import metrics


def library_normalize_log1p(counts):
    """Library-size normalize each spot to 1e4 then log1p. Accepts dense or scipy CSR."""
    from scipy import sparse
    if sparse.issparse(counts):
        m = counts.tocsr().astype(np.float32)
        lib = np.asarray(m.sum(axis=1)).ravel()
        lib[lib == 0] = 1.0
        inv = sparse.diags(1e4 / lib)
        m = inv @ m
        m.data = np.log1p(m.data)
        return m
    counts = np.asarray(counts, float)
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(counts / lib * 1e4)


def select_top_genes(normed_train, n_genes):
    """Top-n genes by mean normalized expression on TRAIN data only. CSR-aware."""
    from scipy import sparse
    if sparse.issparse(normed_train):
        means = np.asarray(normed_train.mean(axis=0)).ravel()
    else:
        means = normed_train.mean(axis=0)
    return np.argsort(-means)[:n_genes]


class Standardizer:
    def fit(self, X):
        self.mu = np.asarray(X, float).mean(axis=0)
        self.sd = np.asarray(X, float).std(axis=0)
        self.sd[self.sd == 0] = 1.0
        return self

    def transform(self, X):
        return (np.asarray(X, float) - self.mu) / self.sd


def fit_ridge_bank(X, Y, lam):
    """Vectorized ridge for all output columns at once (shares the Gram solve).
    beta = (Xa^T Xa + lam P)^-1 Xa^T Y with unpenalized intercept."""
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    n = X.shape[0]
    Xa = np.hstack([X, np.ones((n, 1))])
    P = np.eye(Xa.shape[1]) * lam; P[-1, -1] = 0.0
    return np.linalg.solve(Xa.T @ Xa + P, Xa.T @ Y)


def predict_bank(X, bank):
    X = np.asarray(X, float)
    Xa = np.hstack([X, np.ones((X.shape[0], 1))])
    return Xa @ bank


def per_gene_pearson(P, T):
    """Vectorized per-column Pearson between prediction and truth matrices."""
    P = np.asarray(P, float); T = np.asarray(T, float)
    Pc = P - P.mean(axis=0); Tc = T - T.mean(axis=0)
    num = (Pc * Tc).sum(axis=0)
    den = np.sqrt((Pc ** 2).sum(axis=0) * (Tc ** 2).sum(axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    return r


def choose_lambda_inner(X, Y, groups, grid=(0.1, 1.0, 10.0, 100.0)):
    """Pick ridge lambda by group-held-out inner CV on training patients."""
    groups = np.asarray(groups)
    best, best_score = grid[0], -np.inf
    uniq = np.unique(groups)
    for lam in grid:
        rs = []
        for g in uniq:
            tr = groups != g; te = ~tr
            if tr.sum() < 20 or te.sum() < 5:
                continue
            bank = fit_ridge_bank(X[tr], Y[tr], lam)
            P = predict_bank(X[te], bank)
            rs.append(np.nanmedian(per_gene_pearson(P, Y[te])))
        score = np.mean(rs) if rs else -np.inf
        if score > best_score:
            best, best_score = lam, score
    return best


def choose_alpha_inner(sections_pred, sections_true, sections_coords, length_scale, grid=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)):
    """Pick smoothing alpha by inner predictions on training sections.
    sections_* are lists of per-section arrays (predictions, truth, coords)."""
    best, best_score = 0.0, -np.inf
    for a in grid:
        rs = []
        for P, T, C in zip(sections_pred, sections_true, sections_coords):
            d = np.sqrt(((C[:, None, :] - C[None, :, :]) ** 2).sum(-1))
            W = metrics.gaussian_kernel_weights(d, length_scale)
            S = metrics.spatial_smooth(P, W, a)
            rs.append(np.nanmedian(per_gene_pearson(S, T)))
        score = np.mean(rs) if rs else -np.inf
        if score > best_score:
            best, best_score = a, score
    return best
