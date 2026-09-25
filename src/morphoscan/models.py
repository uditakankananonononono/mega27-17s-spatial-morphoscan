"""Model pipeline for MorphoScan Track A (expression) and Track B (tumor scan).
All selection (genes, scalers, lambda, alpha) happens inside training folds only."""
from __future__ import annotations
import numpy as np
from . import metrics


def library_normalize_log1p(counts):
    """Library-size normalize each spot to 1e4 then log1p."""
    counts = np.asarray(counts, float)
    lib = counts.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    return np.log1p(counts / lib * 1e4)


def select_top_genes(normed_train, n_genes):
    """Top-n genes by mean normalized expression on TRAIN data only."""
    order = np.argsort(-normed_train.mean(axis=0))
    return order[:n_genes]


class Standardizer:
    def fit(self, X):
        self.mu = np.asarray(X, float).mean(axis=0)
        self.sd = np.asarray(X, float).std(axis=0)
        self.sd[self.sd == 0] = 1.0
        return self

    def transform(self, X):
        return (np.asarray(X, float) - self.mu) / self.sd


def fit_ridge_bank(X, Y, lam):
    """Fit one ridge per output column. Returns (n_feat+1, n_out) beta bank."""
    return np.stack([metrics.ridge_fit(X, Y[:, j], lam) for j in range(Y.shape[1])], axis=1)


def predict_bank(X, bank):
    return np.stack([metrics.ridge_predict(X, bank[:, j]) for j in range(bank.shape[1])], axis=1)


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
            r = [metrics.pearson_r(P[:, j], Y[te, j]) for j in range(Y.shape[1])]
            rs.append(np.median(r))
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
            r = [metrics.pearson_r(S[:, j], T[:, j]) for j in range(T.shape[1])]
            rs.append(np.median(r))
        score = np.mean(rs) if rs else -np.inf
        if score > best_score:
            best, best_score = a, score
    return best
