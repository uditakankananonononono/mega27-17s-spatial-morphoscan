"""Track B: tumor-region scanning from morphology, leave-one-patient-out.
Logistic regression + HistGradientBoosting comparison. Writes results/trackB_results.json,
saves final scanner model for the CLI (results/scanner_model.npz)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from src.morphoscan import metrics

CACHE = "features_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)

rows = []
for sec in meta["section"]:
    p = os.path.join(CACHE, f"{sec}.npz")
    if not os.path.exists(p): continue
    z = np.load(p, allow_pickle=True)
    X = z["X"]; ok = np.isfinite(X).all(axis=1)
    r = meta[meta["section"] == sec].iloc[0]
    rows.append(pd.DataFrame(dict(section=sec, patient=r["patient"], subtype=r["type"],
                                  y=z["labels"][ok], row=np.arange(ok.sum()))))
    rows[-1]["X"] = list(X[ok])
D = pd.concat(rows, ignore_index=True)
X = np.stack(D["X"].values); y = D["y"].values; groups = D["patient"].values
print(f"spots={len(y)} tumor_frac={y.mean():.2f} patients={len(set(groups))}")

out = {"folds": []}
for pat in sorted(set(groups)):
    tr = groups != pat; te = ~tr
    lr = LogisticRegression(max_iter=2000, class_weight="balanced")
    lr.fit(X[tr], y[tr])
    p = lr.predict_proba(X[te])[:, 1]
    gb = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08)
    gb.fit(X[tr], y[tr])
    pg = gb.predict_proba(X[te])[:, 1]
    out["folds"].append(dict(patient=pat, n=int(te.sum()), tumor_frac=float(y[te].mean()),
                             auc_lr=metrics.auc_mann_whitney(p, y[te]),
                             auc_gb=metrics.auc_mann_whitney(pg, y[te]),
                             brier_lr=metrics.brier_score(p, y[te])))
    print(pat, out["folds"][-1])

aucs = [f["auc_lr"] for f in out["folds"] if not np.isnan(f["auc_lr"])]
aucsg = [f["auc_gb"] for f in out["folds"] if not np.isnan(f["auc_gb"])]
out["median_auc_lr"] = float(np.median(aucs)); out["median_auc_gb"] = float(np.median(aucsg))
out["mean_auc_lr"] = float(np.mean(aucs)); out["mean_auc_gb"] = float(np.mean(aucsg))

# pooled held-out predictions for operating point
probs = np.full(len(y), np.nan)
for pat in sorted(set(groups)):
    tr = groups != pat; te = ~tr
    lr = LogisticRegression(max_iter=2000, class_weight="balanced").fit(X[tr], y[tr])
    probs[te] = lr.predict_proba(X[te])[:, 1]
pred = (probs >= 0.5).astype(int)
tp = int(((pred == 1) & (y == 1)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
fn = int(((pred == 0) & (y == 1)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
out["operating_point"] = dict(sensitivity=tp / max(1, tp + fn), specificity=tn / max(1, tn + fp),
                              brier=metrics.brier_score(probs, y))

# final scanner model on all data for the CLI
mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1
lr = LogisticRegression(max_iter=3000, class_weight="balanced").fit((X - mu) / sd, y)
np.savez("results/scanner_model.npz", coef=lr.coef_.ravel(), intercept=lr.intercept_, mu=mu, sd=sd)
json.dump(out, open("results/trackB_results.json", "w"), indent=1)
print("TRACKB DONE", out["median_auc_lr"], out["median_auc_gb"])
