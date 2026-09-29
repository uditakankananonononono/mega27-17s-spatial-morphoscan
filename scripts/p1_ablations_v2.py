"""Tier-2 #6 morphology ablations on corrected-geometry features (features_cache2).
PRE-DECLARED in this header before results: diagnostic, no gate, no threshold
decision. Groups over the 41 handcrafted features (deterministic ordering in
src/morphoscan/features.py): color = rgb_* (idx 0-11); stain = hema/eosin/dab_*
(12-20); texture = glcm_* (21-36); entropy = gray_entropy + hema_entropy (37-38);
density = hema_dense_frac + grad_mag_mean (39-40). Protocol identical to
Track A M1 (scripts/run_trackA.py): HER2+ patients, LOPO, top-250 genes
train-selected, library-size 1e4 + log1p, train Standardizer, ridge lambda by
inner patient-held-out CV. Report per-group median per-gene Pearson (pooled
fold medians, same aggregation as run_trackA) vs the full-41 baseline, plus
scanner (LogReg C=1.0) patient-held-out AUC per group.
Run: python3 scripts/p1_ablations_v2.py -> results/p1_ablations_v2.json"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from src.morphoscan import models

CACHE = "features_cache2"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
N_GENES = 250
GROUPS = {"full": list(range(41)), "color": list(range(0, 12)), "stain": list(range(12, 21)),
          "texture": list(range(21, 37)), "entropy": [37, 38], "density": [39, 40]}

meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(CACHE, f"{s}.npz"))]
gene_sets = []
for sec in all_secs:
    z = np.load(os.path.join(CACHE, f"{sec}.npz"), allow_pickle=True)
    gene_sets.append(set(str(g) for g in z["genes"])); z.close()
universe = sorted(set.intersection(*gene_sets))
upos = {g: i for i, g in enumerate(universe)}
her2_secs = [s for s in meta[meta["is_her2"]]["section"] if s in all_secs]

def load(secs):
    out = {}
    for sec in secs:
        z = np.load(os.path.join(CACHE, f"{sec}.npz"), allow_pickle=True)
        X = z["X"]; ok = np.isfinite(X).all(axis=1)
        csr = z["counts_sparse"].item()[ok]
        genes = [str(g) for g in z["genes"]]
        labels = z["labels"][ok]
        z.close()
        lib = np.asarray(csr.sum(axis=1)).ravel(); lib[lib == 0] = 1.0
        lut = np.array([upos.get(g, -1) for g in genes])
        m = csr.tocoo(); nl = lut[m.col]; keep = nl >= 0
        N = sparse.csr_matrix((m.data[keep].astype(np.float32), (m.row[keep], nl[keep])),
                              shape=(m.shape[0], len(upos)), dtype=np.float32)
        N.data *= np.repeat((1e4 / lib).astype(np.float32), np.diff(N.indptr))
        np.log1p(N.data, out=N.data)
        out[sec] = dict(X=X[ok], N=N, labels=labels)
        del csr, m
    import gc; gc.collect()
    return out

sections = load(her2_secs)
her2 = meta[meta["is_her2"]]
patients = sorted(her2["patient"].unique())
print(f"{len(sections)} sections, {len(patients)} patients, universe {len(universe)}", flush=True)

def assemble(pats):
    Xs, Ns, Ls, Gs = [], [], [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Xs.append(s["X"]); Ns.append(s["N"]); Ls.append(s["labels"])
            Gs += [pat] * len(s["X"])
    return np.vstack(Xs), sparse.vstack(Ns), np.concatenate(Ls), np.array(Gs)

res = {"item": "Tier-2 #6 morphology ablations (features_cache2, corrected geometry)",
       "groups": {k: len(v) for k, v in GROUPS.items()}, "per_group": {}, "per_fold": {}}
for gname, cols in GROUPS.items():
    fold_med, fold_auc = {}, {}
    t0 = time.time()
    for test_pat in patients:
        train_pats = [p for p in patients if p != test_pat]
        Xtr, Ntr, Ltr, Gtr = assemble(train_pats)
        top = models.select_top_genes(Ntr, N_GENES)
        Ytr = Ntr[:, top].toarray()
        sc = models.Standardizer().fit(Xtr[:, cols])
        Ztr = sc.transform(Xtr[:, cols])
        lam = models.choose_lambda_inner(Ztr, Ytr, Gtr)
        bank = models.fit_ridge_bank(Ztr, Ytr, lam)
        scan = LogisticRegression(max_iter=2000, C=1.0).fit(Ztr, Ltr)
        rs, as_ = [], []
        for sec in her2[her2["patient"] == test_pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Z = sc.transform(s["X"][:, cols])
            rs.append(models.per_gene_pearson(models.predict_bank(Z, bank), s["N"][:, top].toarray()))
            p = scan.predict_proba(Z)[:, 1]
            if s["labels"].any() and (~s["labels"].astype(bool)).any():
                from src.morphoscan import metrics as _m
                as_.append(_m.auc_mann_whitney(p, s["labels"]))
        fold_med[test_pat] = float(np.nanmedian(np.concatenate(rs)))
        fold_auc[test_pat] = float(np.median(as_)) if as_ else None
        del Xtr, Ntr, Ytr, Ztr
        import gc; gc.collect()
    meds = sorted(fold_med.values())
    aucs = sorted(v for v in fold_auc.values() if v is not None)
    res["per_group"][gname] = {"median_r": float(np.median(meds)),
                               "median_scanner_auc": float(np.median(aucs)) if aucs else None}
    res["per_fold"][gname] = {"r": fold_med, "auc": fold_auc}
    print(f"{gname}: median_r {res['per_group'][gname]['median_r']:.4f} "
          f"auc {res['per_group'][gname]['median_scanner_auc']} ({time.time()-t0:.0f}s)", flush=True)
json.dump(res, open("results/p1_ablations_v2.json", "w"), indent=1)
print("saved results/p1_ablations_v2.json", flush=True)
