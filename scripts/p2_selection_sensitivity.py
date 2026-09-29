"""Amendment-queue Tier-2 #11: sensitivity of per-gene results to the train-only
gene-selection rule. Pre-declared design (2026-09-29, before any run): four rules -
(a) top-100 by train mean normalized expression, (b) top-250 mean (the committed
reference), (c) top-500 mean, (d) top-250 by train variance. Same LOPO folds, same
Standardizer, same per-fold committed lambda, same ridge bank. Metric: per-gene
median Pearson on the held-out patient. Question: does the selection rule move the
headline median? Reported as measured either way. Writes
results/p2_selection_sensitivity.json. Checkpointed per fold."""
import os, sys, json, time, csv as _csv, gc
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from scipy import sparse
from src.morphoscan import models
import p2_common as C

RULES = {"top100_mean": 100, "top250_mean": 250, "top500_mean": 500, "top250_var": 250}
OUT = "results/p2_selection_sensitivity.json"
state = json.load(open(OUT)) if os.path.exists(OUT) else {"folds": [], "design": {
    "rules": RULES, "note": "train-only selection; committed reference = top250_mean"}}
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

def per_gene_median(P, T):
    rs = []
    for j in range(T.shape[1]):
        if T[:, j].std() > 0 and P[:, j].std() > 0:
            rs.append(np.corrcoef(P[:, j], T[:, j])[0, 1])
    return float(np.nanmedian(rs))

done = {f["patient"] for f in state["folds"]}
for test_pat in patients:
    if test_pat in done: continue
    t0 = time.time()
    lam = folds[test_pat]
    train_pats = [p for p in patients if p != test_pat]
    train_secs = [sec for pat in train_pats for sec in her2[her2["patient"] == pat]["section"]
                  if sec in sections]
    # memory-lean: accumulate gene stats per section, never a full sparse vstack
    sum_x, sum_x2, n_tot = None, None, 0
    for sec in train_secs:
        N = sections[sec]["N"]
        sx = np.asarray(N.sum(axis=0)).ravel()
        sx2 = np.asarray(N.power(2).sum(axis=0)).ravel()
        sum_x = sx if sum_x is None else sum_x + sx
        sum_x2 = sx2 if sum_x2 is None else sum_x2 + sx2
        n_tot += N.shape[0]
    means = sum_x / n_tot
    variances = sum_x2 / n_tot - means ** 2
    sel = {"top100_mean": np.argsort(-means)[:100],
           "top250_mean": np.argsort(-means)[:250],
           "top500_mean": np.argsort(-means)[:500],
           "top250_var": np.argsort(-variances)[:250]}
    Xtr = np.vstack([sections[s]["X"] for s in train_secs])
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    Zte = np.vstack([sc.transform(sections[s]["X"]) for s in test_secs])
    rec = {"patient": test_pat}
    for rule, idx in sel.items():
        Ytr = np.vstack([sections[s]["N"][:, idx].toarray() for s in train_secs])
        Tte = np.vstack([sections[s]["N"][:, idx].toarray() for s in test_secs])
        bank = models.fit_ridge_bank(Ztr, Ytr, lam)
        rec[rule] = per_gene_median(models.predict_bank(Zte, bank), Tte)
        del Ytr, Tte, bank
    ov = len(set(sel["top250_mean"].tolist()) & set(sel["top250_var"].tolist()))
    rec["overlap_mean_vs_var_250"] = ov
    state["folds"].append(rec)
    json.dump(state, open(OUT, "w"))
    print(f"{test_pat} " + " ".join(f"{k} {v:.4f}" for k, v in rec.items() if k.startswith("top"))
          + f" overlap {ov} ({time.time()-t0:.0f}s)", flush=True)
    del Xtr, Ztr, Zte
    gc.collect()

fs = state["folds"]
summ = {}
for rule in RULES:
    v = [f[rule] for f in fs]
    summ[rule] = {"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v))}
summ["overlap_mean_vs_var_250_median"] = float(np.median([f["overlap_mean_vs_var_250"] for f in fs]))
state["summary"] = summ
print(json.dumps(summ, indent=1))
json.dump(state, open(OUT, "w"))
