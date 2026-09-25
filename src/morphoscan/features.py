"""Handcrafted morphology features from H&E patches.
Implements Beer-Lambert color deconvolution (Ruifrok-Johnston H&E/DAB matrix),
GLCM texture (contrast, homogeneity, energy, correlation), Shannon entropy,
and color moments. No deep network: interpretable by design."""
from __future__ import annotations
import numpy as np

# Ruifrok & Johnston H&E-DAB stain matrix (rows: H, E, DAB), columns: R,G,B
STAIN_MATRIX = np.array([
    [0.650, 0.704, 0.286],
    [0.072, 0.990, 0.105],
    [0.268, 0.570, 0.776],
])


def optical_density(img, i0=255.0, eps=1.0):
    """Beer-Lambert: OD = -ln((I+eps)/I0) per channel."""
    img = np.asarray(img, float)
    return -np.log((img + eps) / i0)


def color_deconvolution(img, matrix=STAIN_MATRIX):
    """Concentration map C = OD @ inv(M). Returns HxWx3 concentrations."""
    od = optical_density(img)
    inv = np.linalg.inv(matrix)
    h, w, _ = od.shape
    return od.reshape(-1, 3) @ inv.T


def glcm(gray, distances=(1,), angles=(0, np.pi / 2), levels=16):
    """Gray-level co-occurrence matrices (symmetric, normalized)."""
    gray = np.asarray(gray)
    q = np.clip((gray.astype(float) / (gray.max() + 1e-9) * (levels - 1)).astype(int), 0, levels - 1)
    h, w = q.shape
    mats = []
    for d in distances:
        for a in angles:
            dx, dy = int(round(d * np.cos(a))), int(round(d * np.sin(a)))
            m = np.zeros((levels, levels))
            xs = slice(max(0, -dx), min(w, w - dx))
            ys = slice(max(0, -dy), min(h, h - dy))
            xs2 = slice(max(0, dx), min(w, w + dx))
            ys2 = slice(max(0, dy), min(h, h + dy))
            a1 = q[ys, xs].ravel(); a2 = q[ys2, xs2].ravel()
            m = np.bincount(a1 * levels + a2, minlength=levels * levels).reshape(levels, levels).astype(float)
            m = m + m.T
            s = m.sum()
            mats.append(m / s if s > 0 else m)
    return mats


def glcm_props(m):
    """Contrast=sum((i-j)^2 p), homogeneity=sum(p/(1+|i-j|)), energy=sum(p^2),
    correlation with means/stds of marginals."""
    n = m.shape[0]
    i, j = np.indices(m.shape)
    contrast = float(((i - j) ** 2 * m).sum())
    homogeneity = float((m / (1.0 + np.abs(i - j))).sum())
    energy = float((m ** 2).sum())
    pi = m.sum(axis=1); pj = m.sum(axis=0)
    mi = (i[:, 0] * pi).sum(); mj = (j[0, :] * pj).sum()
    si = np.sqrt((((i[:, 0] - mi) ** 2) * pi).sum())
    sj = np.sqrt((((j[0, :] - mj) ** 2) * pj).sum())
    if si * sj == 0:
        corr = 0.0
    else:
        corr = float((((i - mi) * (j - mj) * m).sum()) / (si * sj))
    return contrast, homogeneity, energy, corr


def shannon_entropy(probs):
    """H = -sum p log2 p."""
    p = np.asarray(probs, float)
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def patch_features(img, glcm_levels=16):
    """Full feature vector for one RGB patch (HxWx3, uint8)."""
    img = np.asarray(img)
    feats = {}
    # 1) color moments per RGB channel (mean, std, skew, kurtosis) = 12
    for c, name in enumerate("rgb"):
        v = img[..., c].astype(float).ravel()
        mu, sd = v.mean(), v.std() + 1e-9
        feats[f"{name}_mean"] = mu; feats[f"{name}_std"] = sd
        feats[f"{name}_skew"] = float((((v - mu) / sd) ** 3).mean())
        feats[f"{name}_kurt"] = float((((v - mu) / sd) ** 4).mean())
    # 2) stain concentrations: mean/std/frac-positive per stain = 9
    C = color_deconvolution(img)
    h, w = img.shape[:2]
    C = C.reshape(h, w, 3)
    for c, name in enumerate(["hema", "eosin", "dab"]):
        v = C[..., c].ravel()
        feats[f"{name}_mean"] = float(v.mean())
        feats[f"{name}_std"] = float(v.std())
        feats[f"{name}_fracpos"] = float((v > v.mean()).mean())
    # 3) grayscale GLCM texture: 4 props x 2 distances x 2 angles = 16
    gray = (0.299 * img[..., 0] + 0.587 * img[..., 1] + 0.114 * img[..., 2])
    mats = glcm(gray, distances=(1, 2), angles=(0, np.pi / 2), levels=glcm_levels)
    for k, m in enumerate(mats):
        con, hom, ene, cor = glcm_props(m)
        feats[f"glcm{k}_contrast"] = con; feats[f"glcm{k}_homog"] = hom
        feats[f"glcm{k}_energy"] = ene; feats[f"glcm{k}_corr"] = cor
    # 4) entropies: gray + hematoxylin channel = 2
    gh = np.histogram(gray, bins=32, density=True)[0]; gh = gh / (gh.sum() + 1e-12)
    feats["gray_entropy"] = shannon_entropy(gh)
    hh = np.histogram(C[..., 0], bins=32, density=True)[0]; hh = hh / (hh.sum() + 1e-12)
    feats["hema_entropy"] = shannon_entropy(hh)
    # 5) hematoxylin "nuclear density" proxy: fraction of pixels with high H
    feats["hema_dense_frac"] = float((C[..., 0] > np.percentile(C[..., 0], 90)).mean())
    # 6) patch sharpness proxy: gradient magnitude mean
    gy, gx = np.gradient(gray)
    feats["grad_mag_mean"] = float(np.sqrt(gx ** 2 + gy ** 2).mean())
    return feats


FEATURE_NAMES = None


def feature_names():
    """Deterministic feature ordering."""
    names = []
    for n in "rgb":
        names += [f"{n}_mean", f"{n}_std", f"{n}_skew", f"{n}_kurt"]
    for n in ["hema", "eosin", "dab"]:
        names += [f"{n}_mean", f"{n}_std", f"{n}_fracpos"]
    for k in range(4):
        names += [f"glcm{k}_contrast", f"glcm{k}_homog", f"glcm{k}_energy", f"glcm{k}_corr"]
    names += ["gray_entropy", "hema_entropy", "hema_dense_frac", "grad_mag_mean"]
    return names


def patch_vector(img):
    f = patch_features(img)
    return np.array([f[n] for n in feature_names()])
