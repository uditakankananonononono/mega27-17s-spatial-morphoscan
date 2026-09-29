"""Tier-2 #15 scanner calibration diagnostics (pre-declared in this header before
running; descriptive, no gate, no threshold decision depends on it).
Replicates the committed scanner EXACTLY (scripts/p2_model.py: features_cache +
embeddings_cache join, train Standardizer, LogisticRegression C=1.0
max_iter=2000, patient-held-out folds) and collects out-of-fold per-spot
probabilities on held-out sections. Diagnostics: Brier vs prevalence baseline,
10-bin ECE, reliability table. Replay check: per-section AUCs recomputed from
the collected probabilities must match results/p2_scanner_partial.csv within
1e-6 before any calibration number is trusted.
Run: python3 scripts/p2_scanner_calibration.py -> results/p2_scanner_calibration.json"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from src.morphoscan import models

FC = "features_cache"; EC = "embeddings_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")

all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(FC, f"{s}.npz"))
            and os.path.exists(os.path.join(EC, f"{s}.npz"))]
her2_secs = [s for s in meta[meta["is_her2"]]["section"] if s in all_secs]

# universe not needed for scanner (embeddings only) but labels come from features_cache
def load_labels_X(secs):
    out = {}
    for sec in secs:
        fp = os.path.join(FC, f"{sec}.npz"); ep = os.path.join(EC, f"{sec}.npz")
        if not (os.path.exists(fp) and os.path.exists(ep)): continue
        zf = np.load(fp, allow_pickle=True); ze = np.load(ep, allow_pickle=True)
        f_ids = [str(s) for s in zf["spot_ids"]]
        e_ids = [str(s) for s in ze["spot_ids"]]
        e_pos = {s: i for i, s in enumerate(e_ids)}
        order = [e_pos[s] for s in f_ids if s in e_pos]
        keep_f = [i for i, s in enumerate(f_ids) if s in e_pos]
        if len(order) < 20:
            zf.close(); ze.close(); continue
        X = ze["X"][order].astype(np.float32)
        labels = zf["labels"][keep_f]
        zf.close(); ze.close()
        ok = np.isfinite(X).all(axis=1)
        out[sec] = dict(X=X[ok], labels=labels[ok])
    return out

sections = load_labels_X(her2_secs)
her2 = meta[meta["is_her2"]]
patients = sorted(her2["patient"].unique())
print(f"sections {len(sections)} patients {len(patients)}", flush=True)

probs_all, labels_all, per_sec = [], [], {}
for test_pat in patients:
    train_secs = [s for p in patients if p != test_pat
                  for s in her2[her2["patient"] == p]["section"] if s in sections]
    Xtr = np.vstack([sections[s]["X"] for s in train_secs])
    Ltr = np.concatenate([sections[s]["labels"] for s in train_secs])
    sc = models.Standardizer().fit(Xtr)
    scan = LogisticRegression(max_iter=2000, C=1.0)
    scan.fit(sc.transform(Xtr), Ltr)
    for sec in her2[her2["patient"] == test_pat]["section"]:
        if sec not in sections: continue
        s = sections[sec]
        p = scan.predict_proba(sc.transform(s["X"]))[:, 1]
        probs_all.append(p); labels_all.append(s["labels"])
        per_sec[sec] = float(roc_auc_score(s["labels"], p)) if s["labels"].any() and (~s["labels"].astype(bool)).any() else None
    print(f"fold {test_pat} done", flush=True)

P = np.concatenate(probs_all); Y = np.concatenate(labels_all).astype(float)

# replay check vs committed scanner AUCs
committed = {r["section"]: (float(r["auc"]) if r["auc"] else None)
             for r in __import__("csv").DictReader(open("results/p2_scanner_partial.csv"))}
diffs = {s: abs(per_sec.get(s, -1) - a) for s, a in committed.items() if a is not None and per_sec.get(s) is not None}
max_diff = max(diffs.values()) if diffs else None
assert max_diff is not None and max_diff < 1e-6, f"replay mismatch {max_diff}"
print(f"replay OK: {len(diffs)} sections match, max |dAUC| {max_diff:.2e}", flush=True)

prev = Y.mean()
brier = float(np.mean((P - Y) ** 2))
brier_base = float(np.mean((prev - Y) ** 2))
bins = np.linspace(0, 1, 11)
rel = []
ece = 0.0
for lo, hi in zip(bins[:-1], bins[1:]):
    m = (P >= lo) & (P < hi if hi < 1 else P <= hi)
    n = int(m.sum())
    if n == 0: continue
    mp = float(P[m].mean()); emp = float(Y[m].mean())
    rel.append({"bin": [float(lo), float(hi)], "n": n, "mean_pred": mp, "empirical": emp, "gap": mp - emp})
    ece += n / len(P) * abs(mp - emp)
res = {"item": "Tier-2 #15 scanner calibration (descriptive diagnostics)",
       "n_spots": int(len(P)), "n_patients": len(patients), "n_sections": len(per_sec),
       "prevalence": float(prev), "brier": brier, "brier_prevalence_baseline": brier_base,
       "ece_10bin": float(ece), "reliability": rel,
       "replay_check": {"sections_matched": len(diffs), "max_abs_dauc": max_diff},
       "protocol": "out-of-fold per-spot probabilities, exact committed scanner recipe (LogReg C=1.0, train Standardizer, patient-held-out)"}
json.dump(res, open("results/p2_scanner_calibration.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "reliability"}, indent=1))
