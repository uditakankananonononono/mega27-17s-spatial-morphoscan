"""P3-SN1 (prespec results/p3_spatialnull_prespec.json, commit 051b175).
Spatial-structure-preserving null for the frozen patient-held-out scanner AUC.
Classifier replicates scripts/p2_scanner_labelperm.py exactly; only the null differs."""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from src.morphoscan import models
import p2_common as C

B = 500; SEED = 20261007; TOL = 0.05; MAXTRY = 5000; NRES = 10000; NBOOT = 10000
OUT = "results/p3_spatialnull.json"
rng = np.random.default_rng(SEED)

def auc(y, p):
    r = rankdata(p); n1 = y.sum(); n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

def null_labels(px, y, rng):
    c = px.mean(axis=0); tf = y.mean(); out = []; tries = 0
    tree_cache = None
    while len(out) < B and tries < MAXTRY:
        tries += 1
        th = rng.uniform(0, 2 * np.pi); refl = rng.random() < 0.5
        q = px - c
        if refl: q = q * np.array([-1.0, 1.0])
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        moved = q @ R.T + c
        _, idx = cKDTree(moved).query(px)   # each real spot takes label of nearest moved labeled spot
        yn = y[idx]
        if abs(yn.mean() - tf) <= TOL and 0 < yn.sum() < len(yn):
            out.append(yn)
    return out, tries

meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz")) and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs); upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())
committed = {(r["patient"], r["section"]): float(r["auc"]) for r in _csv.DictReader(open("results/p2_scanner_partial.csv")) if r["auc"] != ""}

def assemble(pats):
    Xs, Ls = [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec in sections:
                Xs.append(sections[sec]["X"]); Ls.append(sections[sec]["labels"])
    return np.vstack(Xs), np.concatenate(Ls)

rows = []; devs = []
for tp in patients:
    t0 = time.time()
    Xtr, Ltr = assemble([p for p in patients if p != tp])
    sc = models.Standardizer().fit(Xtr)
    clf = LogisticRegression(max_iter=2000, C=1.0).fit(sc.transform(Xtr), Ltr.astype(int))
    for sec in her2[her2["patient"] == tp]["section"]:
        if sec not in sections: continue
        s = sections[sec]; y = s["labels"].astype(int)
        if not (y.any() and (~y.astype(bool)).any()): continue
        p = clf.predict_proba(sc.transform(s["X"]))[:, 1]
        a = auc(y, p)
        if (tp, sec) in committed: devs.append(abs(a - committed[(tp, sec)]))
        nl, tries = null_labels(s["px"].astype(float), y, rng)
        na = [auc(l, p) for l in nl]
        rows.append({"patient": tp, "section": sec, "n_spots": int(len(y)), "tumor_frac": float(y.mean()),
                     "obs_auc": a, "n_null": len(na), "tries": tries, "null_aucs": na,
                     "p": float((1 + sum(x >= a for x in na)) / (1 + len(na))) if na else None})
    print(tp, f"{time.time()-t0:.0f}s", flush=True)

res = {"prespec": "results/p3_spatialnull_prespec.json@051b175", "seed": SEED,
       "max_abs_dev_vs_committed": float(max(devs)) if devs else None}
if not devs or max(devs) > 1e-6:
    # LogReg is sensitive to ~1e-8 embedding deltas (see labelperm note); report, decide per prespec
    res["verdict"] = "REPLAY_HALT"
else:
    ok = [r for r in rows if r["n_null"] >= B]
    res["n_sections"] = len(rows); res["n_null_feasible"] = len(ok)
    res["infeasible_sections"] = [r["section"] for r in rows if r["n_null"] < B]
    obs = np.array([r["obs_auc"] for r in ok]); N = np.array([r["null_aucs"] for r in ok])
    obs_med = float(np.median(obs))
    med_null = np.array([np.median(N[np.arange(len(ok)), rng.integers(0, B, len(ok))]) for _ in range(NRES)])
    p_prim = float((1 + (med_null >= obs_med).sum()) / (1 + NRES))
    excess = obs - np.median(N, axis=1)
    pats = np.array([r["patient"] for r in ok]); up = np.unique(pats)
    bs = []
    for _ in range(NBOOT):
        pick = rng.choice(up, len(up)); idx = np.concatenate([np.where(pats == q)[0] for q in pick])
        bs.append(np.median(excess[idx]))
    ps = np.array([r["p"] for r in ok])
    o = np.argsort(ps); adj = np.minimum.accumulate((ps[o] * len(ps) / (np.arange(len(ps)) + 1))[::-1])[::-1]
    q = np.empty(len(ps)); q[o] = np.clip(adj, 0, 1)
    me = float(np.median(excess))
    res.update({"obs_median_auc": obs_med, "null_median_of_medians": float(np.median(med_null)),
                "null_median_p95": float(np.percentile(med_null, 95)), "primary_p": p_prim,
                "pooled_null_auc": {"mean": float(N.mean()), "median": float(np.median(N)), "p95": float(np.percentile(N, 95))},
                "median_excess": me, "excess_ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                "sections_p_lt_0.05": int((ps < 0.05).sum()), "sections_bh_q_lt_0.05": int((q < 0.05).sum())})
    res["verdict"] = "PASS" if (p_prim < 0.01 and me >= 0.15) else ("PARTIAL" if p_prim < 0.01 else "FAIL")
res["sections"] = [{k: v for k, v in r.items() if k != "null_aucs"} | {"null_median": float(np.median(r["null_aucs"])) if r["null_aucs"] else None} for r in rows]
json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "sections"}, indent=1))
