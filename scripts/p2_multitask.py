"""Amendment-queue Tier-2 #9: multi-task gene prediction vs independent ridge.
Same LOPO folds, same train-selected top-250 genes, same Standardizer as
p2_model.py. Alternatives: (a) MultiTaskElasticNet (joint L1/L2, alpha=1e-3,
l1_ratio=0.5, max_iter=500, tol=2e-3, selection=random - box-feasible settings (2026-09-28:
5000-iter default >35 CPU-min/fold, 1000-iter cyclic >15 CPU-min/fold unfinished
under throttling); n_iter_ recorded per fold for convergence honesty), (b) PLSRegression (25 components), (c) reduced-rank ridge
(rank 25, from the ridge bank SVD). Metric: per-gene median Pearson on held-out
patient, compared to committed ridge median_r_m1 per fold. Checkpointed per
fold. Writes results/p2_multitask.json.
2026-09-29 leg reorder: enet never completed fold 1 across three configs
(5000-iter, 1000-iter cyclic, 500-iter random; 26+ CPU-min on the last, VM
frozen ~95% of wall time overnight). PLS/RRR/ridge-refold now run first with
per-fold checkpoints (pass 1); enet is a separate per-fold-checkpointed pass 2
so the multi-task verdict cannot be held hostage by the one infeasible leg."""
import os, sys, json, time, csv as _csv, gc
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from scipy import sparse
from sklearn.linear_model import MultiTaskElasticNet
from sklearn.cross_decomposition import PLSRegression
from src.morphoscan import models
import p2_common as C

N_GENES = 250
RANK = 25
OUT = "results/p2_multitask.json"
state = json.load(open(OUT)) if os.path.exists(OUT) else {"folds": []}
meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())
folds = {r["patient"]: float(r["lam"]) for r in _csv.DictReader(open("results/p2_folds_partial.csv"))}
ridge_m1 = {r["patient"]: float(r["median_r_m1"]) for r in _csv.DictReader(open("results/p2_folds_partial.csv"))}

def per_gene_median(P, T):
    rs = []
    for j in range(T.shape[1]):
        if T[:, j].std() > 0 and P[:, j].std() > 0:
            rs.append(np.corrcoef(P[:, j], T[:, j])[0, 1])
    return float(np.nanmedian(rs))

done = {f["patient"] for f in state["folds"]}

def build_fold(test_pat):
    train_pats = [p for p in patients if p != test_pat]
    Xs, Ns = [], []
    for pat in train_pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]; Xs.append(s["X"]); Ns.append(s["N"])
    Xtr, Ntr = np.vstack(Xs), sparse.vstack(Ns)
    gc.collect()
    top = models.select_top_genes(Ntr, N_GENES)
    Ytr = Ntr[:, top].toarray()
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    Zte = np.vstack([sc.transform(sections[s]["X"]) for s in test_secs])
    Tte = np.vstack([sections[s]["N"][:, top].toarray() for s in test_secs])
    return Ztr, Ytr, Zte, Tte

# PASS 1: cheap multi-task legs (PLS, reduced-rank ridge, ridge refold sanity)
for test_pat in patients:
    if test_pat in done: continue
    t0 = time.time()
    lam = folds[test_pat]
    Ztr, Ytr, Zte, Tte = build_fold(test_pat)
    rec = {"patient": test_pat, "ridge_m1_committed": ridge_m1[test_pat]}
    pls = PLSRegression(n_components=RANK)
    pls.fit(Ztr, Ytr)
    rec["pls_median_r"] = per_gene_median(pls.predict(Zte), Tte)
    print(f"{test_pat} pls {rec['pls_median_r']:.4f} ({time.time()-t0:.0f}s)", flush=True)
    bank = models.fit_ridge_bank(Ztr, Ytr, lam)
    U, S, Vt = np.linalg.svd(bank[:-1], full_matrices=False)
    k = min(RANK, len(S))
    Brr = U[:, :k] @ np.diag(S[:k]) @ Vt[:k]
    Zte_a = np.hstack([Zte, np.ones((Zte.shape[0], 1))])
    rec["rrr_median_r"] = per_gene_median(Zte_a @ np.vstack([Brr, bank[-1:]]), Tte)
    rec["ridge_m1_refold"] = per_gene_median(models.predict_bank(Zte, bank), Tte)
    print(f"{test_pat} rrr {rec['rrr_median_r']:.4f} ridge-refold {rec['ridge_m1_refold']:.4f} ({time.time()-t0:.0f}s)", flush=True)
    state["folds"].append(rec)
    json.dump(state, open(OUT, "w"))
    del Ztr, Ytr, Zte, Tte
    gc.collect()

# PASS 2: elastic-net leg (expensive coordinate descent; per-fold checkpoint so
# the cheap verdict stands even under box throttling; n_iter_ kept for honesty)
for rec in state["folds"]:
    if "elasticnet_median_r" in rec: continue
    test_pat = rec["patient"]
    t0 = time.time()
    Ztr, Ytr, Zte, Tte = build_fold(test_pat)
    en = MultiTaskElasticNet(alpha=1e-3, l1_ratio=0.5, max_iter=500, tol=2e-3, selection="random", random_state=20260929)
    en.fit(Ztr, Ytr)
    rec["elasticnet_median_r"] = per_gene_median(en.predict(Zte), Tte)
    rec["elasticnet_n_iter"] = int(en.n_iter_)  # convergence honesty
    print(f"{test_pat} enet {rec['elasticnet_median_r']:.4f} n_iter {rec['elasticnet_n_iter']} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(state, open(OUT, "w"))
    del Ztr, Ytr, Zte, Tte
    gc.collect()

fs = state["folds"]
summ = {}
enet_done = [f for f in fs if "elasticnet_median_r" in f]
for m in ["ridge_m1_committed", "ridge_m1_refold", "pls_median_r", "rrr_median_r"]:
    v = [f[m] for f in fs]
    summ[m] = {"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v))}
if enet_done:
    v = [f["elasticnet_median_r"] for f in enet_done]
    summ["elasticnet_median_r"] = {"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v)),
                                   "n_folds": len(enet_done)}
wins = {m: sum(1 for f in fs if f[m] > f["ridge_m1_committed"]) for m in ["pls_median_r", "rrr_median_r"]}
if enet_done:
    wins["elasticnet_median_r"] = sum(1 for f in enet_done if f["elasticnet_median_r"] > f["ridge_m1_committed"])
state["summary"] = summ
state["folds_better_than_ridge"] = wins
state["enet_folds_done"] = len(enet_done)
print(json.dumps({"summary": summ, "folds_better_than_ridge": wins, "enet_folds_done": len(enet_done)}, indent=1))
json.dump(state, open(OUT, "w"))
