"""Amendment-queue Tier-2 #12 (gene-level leg): negative-control gene sets.
Pre-declared design (2026-09-29, before any run): per LOPO fold, rank genes by
train mean normalized expression (same rule as select_top_genes); the selected
set is ranks 1-250. Controls are 5 independent draws of 250 genes sampled
uniformly from ranks 251-2250 (expressed but not selected), seed 20260929.
Same Standardizer, same per-fold lambda (results/p2_folds_partial.csv), same
ridge bank. Expectation under the null: control per-gene held-out correlations
center at ~0; the selected-250 distribution is shifted positive. Writes
results/p2_gene_negcontrol.json. Checkpointed per fold."""
import os, sys, json, time, csv as _csv, gc
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from scipy import sparse
from src.morphoscan import models
import p2_common as C

N_GENES = 250
N_DRAWS = 5
SEED = 20260929
POOL_HI = 2250  # sample control genes from train-expression ranks 251..POOL_HI
OUT = "results/p2_gene_negcontrol.json"
state = json.load(open(OUT)) if os.path.exists(OUT) else {"folds": [], "design": {
    "n_genes": N_GENES, "n_draws": N_DRAWS, "seed": SEED,
    "control_pool": f"train-expression ranks {N_GENES+1}..{POOL_HI}"}}
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

def per_gene_rs(P, T):
    rs = []
    for j in range(T.shape[1]):
        if T[:, j].std() > 0 and P[:, j].std() > 0:
            rs.append(np.corrcoef(P[:, j], T[:, j])[0, 1])
    return np.array(rs)

done = {f["patient"] for f in state["folds"]}
for test_pat in patients:
    if test_pat in done: continue
    t0 = time.time()
    lam = folds[test_pat]
    train_pats = [p for p in patients if p != test_pat]
    Xs, Ns = [], []
    for pat in train_pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]; Xs.append(s["X"]); Ns.append(s["N"])
    Xtr, Ntr = np.vstack(Xs), sparse.vstack(Ns)
    gc.collect()
    means = np.asarray(Ntr.mean(axis=0)).ravel()
    ranked = np.argsort(-means)
    pool = ranked[N_GENES:POOL_HI]
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    Zte = np.vstack([sc.transform(sections[s]["X"]) for s in test_secs])
    draw_medians, draw_q = [], []
    for d in range(N_DRAWS):
        rng = np.random.default_rng(SEED + d)
        ctrl = rng.choice(pool, size=N_GENES, replace=False)
        Ytr = Ntr[:, ctrl].toarray()
        Tte = np.vstack([sections[s]["N"][:, ctrl].toarray() for s in test_secs])
        bank = models.fit_ridge_bank(Ztr, Ytr, lam)
        rs = per_gene_rs(models.predict_bank(Zte, bank), Tte)
        draw_medians.append(float(np.median(rs)))
        draw_q.append([float(np.percentile(rs, 5)), float(np.percentile(rs, 95))])
    state["folds"].append({"patient": test_pat, "draw_medians": draw_medians,
                           "draw_p5_p95": draw_q})
    json.dump(state, open(OUT, "w"))
    print(f"{test_pat} control medians {[round(m,4) for m in draw_medians]} ({time.time()-t0:.0f}s)", flush=True)
    del Xtr, Ntr, Ztr, Zte
    gc.collect()

fs = state["folds"]
all_med = [m for f in fs for m in f["draw_medians"]]
state["summary"] = {"n_folds": len(fs), "n_draws_per_fold": N_DRAWS,
                    "control_median_of_medians": float(np.median(all_med)),
                    "control_min": float(np.min(all_med)), "control_max": float(np.max(all_med))}
print(json.dumps(state["summary"], indent=1))
json.dump(state, open(OUT, "w"))
