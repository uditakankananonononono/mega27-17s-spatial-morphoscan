"""Amendment-queue Tier-2 #3: complete spatial null model for M2 (smoothing).
Pre-declared design (2026-09-29, before any run). Real M2 per fold: ridge P1 on
held-out patient sections smoothed as P2 = (1-a) P1 + a W P1, with committed
per-fold alpha (p2_folds_partial.csv) and W = gaussian kernel over pixel
distances, length scale 1.5*pitch (exact p2_model recipe). Three nulls, B=20
draws each, seed 20260929:
(a) value-shuffle: smooth P1[perm] with the real W - destroys spot-value
    alignment, keeps graph;
(b) graph-shuffle: smooth real P1 with W[:, perm] - keeps each row's weight
    marginals, destroys neighborhoods;
(c) distance-shell permutation: W columns permuted within 5 quantile distance
    shells per row - preserves the distance profile, destroys exact geometry.
Expectation if the smoothing gain is spatially specific: null gains
(r_null - r_m1) center at/below 0 while the real gain is positive. Metric:
per-gene median Pearson, pooled across test sections per fold (same as M2).
Writes results/p2_spatial_null.json. Memory-lean; checkpointed per fold."""
import os, sys, json, time, csv as _csv, gc
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from src.morphoscan import models, metrics
import p2_common as C

N_GENES, B, SEED = 250, 20, 20260929
OUT = "results/p2_spatial_null.json"
state = json.load(open(OUT)) if os.path.exists(OUT) else {"folds": [], "design": {
    "nulls": ["value_shuffle", "graph_shuffle", "distance_shell"], "B": B, "seed": SEED}}
meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())
foldrows = {r["patient"]: r for r in _csv.DictReader(open("results/p2_folds_partial.csv"))}

def med_r(P, T):
    rs = []
    for j in range(T.shape[1]):
        if T[:, j].std() > 0 and P[:, j].std() > 0:
            rs.append(np.corrcoef(P[:, j], T[:, j])[0, 1])
    return float(np.nanmedian(rs))

def shell_perm_W(d, rng):
    n = d.shape[0]
    W = np.zeros_like(d)
    for i in range(n):
        di = d[i]
        shells = np.quantile(di[di > 0], [0.2, 0.4, 0.6, 0.8]) if (di > 0).any() else [1, 2, 3, 4]
        shell_id = np.digitize(di, shells)
        w = np.exp(-(di ** 2) / (2.0 * (1.5 * PITCH) ** 2)); w[i] = 0.0
        w_new = np.zeros(n)
        for sh in np.unique(shell_id):
            idx = np.where(shell_id == sh)[0]
            idx = idx[idx != i]
            if len(idx) == 0: continue
            w_new[idx] = rng.permutation(w[idx])
        W[i] = w_new
    s = W.sum(axis=1, keepdims=True); s[s == 0] = 1.0
    return W / s

done = {f["patient"] for f in state["folds"]}
for test_pat in patients:
    if test_pat in done: continue
    t0 = time.time()
    lam = float(foldrows[test_pat]["lam"]); alpha = float(foldrows[test_pat]["alpha"])
    train_pats = [p for p in patients if p != test_pat]
    train_secs = [sec for pat in train_pats for sec in her2[her2["patient"] == pat]["section"]
                  if sec in sections]
    sum_x, n_tot = None, 0
    for sec in train_secs:
        N = sections[sec]["N"]
        sx = np.asarray(N.sum(axis=0)).ravel()
        sum_x = sx if sum_x is None else sum_x + sx
        n_tot += N.shape[0]
    top = np.argsort(-(sum_x / n_tot))[:N_GENES]
    Xtr = np.vstack([sections[s]["X"] for s in train_secs])
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    Ytr = np.vstack([sections[s]["N"][:, top].toarray() for s in train_secs])
    bank = models.fit_ridge_bank(Ztr, Ytr, lam)
    del Xtr, Ztr, Ytr; gc.collect()
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    r1s, r2s = [], []
    null_gains = {k: [] for k in ["value_shuffle", "graph_shuffle", "distance_shell"]}
    for sec in test_secs:
        s = sections[sec]
        P1 = models.predict_bank(sc.transform(s["X"]), bank)
        T = s["N"][:, top].toarray()
        d = np.sqrt(((s["px"][:, None, :] - s["px"][None, :, :]) ** 2).sum(-1))
        PITCH = float(s["pitch"])
        W = metrics.gaussian_kernel_weights(d, 1.5 * PITCH)
        P2 = metrics.spatial_smooth(P1, W, alpha)
        r1s.append(med_r(P1, T)); r2s.append(med_r(P2, T))
        rng = np.random.default_rng(SEED)
        n = P1.shape[0]
        for b in range(B):
            Pv = metrics.spatial_smooth(P1[rng.permutation(n)], W, alpha)
            null_gains["value_shuffle"].append(med_r(Pv, T) - r1s[-1])
            Wg = W[:, rng.permutation(n)]
            null_gains["graph_shuffle"].append(med_r(metrics.spatial_smooth(P1, Wg, alpha), T) - r1s[-1])
            Ws = shell_perm_W(d, rng)
            null_gains["distance_shell"].append(med_r(metrics.spatial_smooth(P1, Ws, alpha), T) - r1s[-1])
    rec = {"patient": test_pat, "alpha": alpha,
           "r_m1": float(np.mean(r1s)), "r_m2": float(np.mean(r2s)),
           "real_gain": float(np.mean(r2s) - np.mean(r1s)),
           "null_gain_median": {k: float(np.median(v)) for k, v in null_gains.items()},
           "null_gain_p95": {k: float(np.percentile(v, 95)) for k, v in null_gains.items()}}
    state["folds"].append(rec)
    json.dump(state, open(OUT, "w"))
    print(f"{test_pat} real_gain {rec['real_gain']:.4f} nulls "
          + " ".join(f"{k} {v:.4f}" for k, v in rec["null_gain_median"].items())
          + f" ({time.time()-t0:.0f}s)", flush=True)
    del bank; gc.collect()

fs = state["folds"]
summ = {"n_folds": len(fs),
        "real_gain_median": float(np.median([f["real_gain"] for f in fs])),
        "real_gain_folds_positive": sum(1 for f in fs if f["real_gain"] > 0)}
for k in ["value_shuffle", "graph_shuffle", "distance_shell"]:
    med = [f["null_gain_median"][k] for f in fs]
    summ[k] = {"median_of_fold_medians": float(np.median(med)),
               "folds_null_below_real": sum(1 for f in fs if f["null_gain_median"][k] < f["real_gain"])}
state["summary"] = summ
print(json.dumps(summ, indent=1))
json.dump(state, open(OUT, "w"))
