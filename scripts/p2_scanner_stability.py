"""Amendment-queue Tier-2 #17: scanner feature-importance stability. Refits the
10 LOPO scanner LogReg classifiers (identical protocol to p2_model.py /
p2_scanner_labelperm.py: Standardizer + LogisticRegression(max_iter=2000, C=1.0)
on train-patient embeddings -> morphology labels) and measures coefficient
stability ACROSS patient-held-out folds: pairwise Spearman of the 768-d
coefficient vectors, sign consistency of the top-|coef| features, plus a
fold-resampling bootstrap (2,000 draws) for the median pairwise correlation.
Writes results/p2_scanner_stability.json. Checkpointed coefficients per fold."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from sklearn.linear_model import LogisticRegression
from scipy.stats import spearmanr
from src.morphoscan import models
import p2_common as C

COEF_CACHE = "results/p2_scanner_coefs.json"
state = json.load(open(COEF_CACHE)) if os.path.exists(COEF_CACHE) else {"coefs": {}}
meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())

for test_pat in patients:
    if test_pat in state["coefs"]: continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xs, Ls = [], []
    for pat in train_pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]; Xs.append(s["X"]); Ls.append(s["labels"])
    Xtr, Ltr = np.vstack(Xs), np.concatenate(Ls)
    import gc; gc.collect()
    sc = models.Standardizer().fit(Xtr)
    scan = LogisticRegression(max_iter=2000, C=1.0)
    scan.fit(sc.transform(Xtr), Ltr.astype(int))
    state["coefs"][test_pat] = scan.coef_[0].tolist()
    json.dump(state, open(COEF_CACHE, "w"))
    print(f"{test_pat}: fit ({time.time()-t0:.0f}s)", flush=True)
    del Xtr, Ltr

pats = sorted(state["coefs"])
B = np.array([state["coefs"][p] for p in pats])
pairs = [(i, j) for i in range(len(pats)) for j in range(i + 1, len(pats))]
rhos = np.array([spearmanr(B[i], B[j]).statistic for i, j in pairs])
K = 50
top_idx = [set(np.argsort(-np.abs(B[i]))[:K]) for i in range(len(pats))]
jacc = np.array([len(top_idx[i] & top_idx[j]) / K for i, j in pairs])
sign = np.array([(np.sign(B[i][list(top_idx[i] & top_idx[j])]) == np.sign(B[j][list(top_idx[i] & top_idx[j])])).mean()
                 if top_idx[i] & top_idx[j] else np.nan for i, j in pairs])
rng = np.random.default_rng(20260928)
boot = []
for _ in range(2000):
    sel = rng.integers(0, len(rhos), len(rhos))
    boot.append(np.median(rhos[sel]))
boot = np.array(boot)
out = {"n_folds": len(pats), "n_features": int(B.shape[1]), "top_k": K,
       "pairwise_spearman_median": float(np.median(rhos)),
       "pairwise_spearman_min": float(rhos.min()), "pairwise_spearman_max": float(rhos.max()),
       "spearman_median_boot95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
       "top50_jaccard_median": float(np.median(jacc)),
       "top50_overlap_sign_consistency_median": float(np.nanmedian(sign))}
print(json.dumps(out, indent=1))
json.dump(out, open("results/p2_scanner_stability.json", "w"), indent=1)
