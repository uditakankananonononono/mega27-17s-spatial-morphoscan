"""Tier-2 #7 technical-artifact control: raw vs Reinhard-normalized handcrafted
features. PRE-DECLARED here before results: diagnostic, no gate. Same Track A
M1 recipe (LOPO, top-250 train-selected, inner-CV lambda) and same scanner
(LogReg C=1.0) run on features_cache2 (raw) vs features_cache3 (normalized);
report both medians and the per-fold paired delta. Single changed variable is
stain normalization, so any delta is attributable to it.
Run: python3 scripts/p1_norm_control.py -> results/p1_norm_control.json"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from src.morphoscan import models, metrics

DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
N_GENES = 250
meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")

def load_cache(cache):
    all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(cache, f"{s}.npz"))]
    gene_sets = []
    for sec in all_secs:
        z = np.load(os.path.join(cache, f"{sec}.npz"), allow_pickle=True)
        gene_sets.append(set(str(g) for g in z["genes"])); z.close()
    universe = sorted(set.intersection(*gene_sets))
    upos = {g: i for i, g in enumerate(universe)}
    out = {}
    for sec in [s for s in meta[meta["is_her2"]]["section"] if s in all_secs]:
        z = np.load(os.path.join(cache, f"{sec}.npz"), allow_pickle=True)
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

her2 = meta[meta["is_her2"]]
patients = sorted(her2["patient"].unique())

def run(sections):
    def assemble(pats):
        Xs, Ns, Ls, Gs = [], [], [], []
        for pat in pats:
            for sec in her2[her2["patient"] == pat]["section"]:
                if sec not in sections: continue
                s = sections[sec]
                Xs.append(s["X"]); Ns.append(s["N"]); Ls.append(s["labels"])
                Gs += [pat] * len(s["X"])
        return np.vstack(Xs), sparse.vstack(Ns), np.concatenate(Ls), np.array(Gs)
    fold_med, fold_auc = {}, {}
    for test_pat in patients:
        train_pats = [p for p in patients if p != test_pat]
        Xtr, Ntr, Ltr, Gtr = assemble(train_pats)
        top = models.select_top_genes(Ntr, N_GENES)
        Ytr = Ntr[:, top].toarray()
        sc = models.Standardizer().fit(Xtr)
        Ztr = sc.transform(Xtr)
        lam = models.choose_lambda_inner(Ztr, Ytr, Gtr)
        bank = models.fit_ridge_bank(Ztr, Ytr, lam)
        scan = LogisticRegression(max_iter=2000, C=1.0).fit(Ztr, Ltr)
        rs, as_ = [], []
        for sec in her2[her2["patient"] == test_pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Z = sc.transform(s["X"])
            rs.append(models.per_gene_pearson(models.predict_bank(Z, bank), s["N"][:, top].toarray()))
            p = scan.predict_proba(Z)[:, 1]
            if s["labels"].any() and (~s["labels"].astype(bool)).any():
                as_.append(metrics.auc_mann_whitney(p, s["labels"]))
        fold_med[test_pat] = float(np.nanmedian(np.concatenate(rs)))
        fold_auc[test_pat] = float(np.median(as_)) if as_ else None
        del Xtr, Ntr, Ytr, Ztr
        import gc; gc.collect()
    return fold_med, fold_auc

res = {"item": "Tier-2 #7 raw vs Reinhard-normalized handcrafted features (single changed variable)",
       "reference": "results/g2x_v2_refstain.json (internal BC23901_C2)"}
for name, cache in [("raw_cache2", "features_cache2"), ("reinhard_cache3", "features_cache3")]:
    t0 = time.time()
    secs = load_cache(cache)
    fm, fa = run(secs)
    res[name] = {"median_r": float(np.median(sorted(fm.values()))),
                 "median_scanner_auc": float(np.median(sorted(v for v in fa.values() if v is not None))),
                 "per_fold_r": fm, "per_fold_auc": fa, "seconds": round(time.time() - t0, 1)}
    print(name, res[name]["median_r"], res[name]["median_scanner_auc"], flush=True)
    del secs
    import gc; gc.collect()
d = {p: res["reinhard_cache3"]["per_fold_r"][p] - res["raw_cache2"]["per_fold_r"][p] for p in patients}
res["paired_delta_r_per_fold"] = d
res["median_paired_delta_r"] = float(np.median(sorted(d.values())))
json.dump(res, open("results/p1_norm_control.json", "w"), indent=1)
print("saved; median paired delta r:", res["median_paired_delta_r"], flush=True)
