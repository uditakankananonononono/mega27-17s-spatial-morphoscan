"""Finish Track B from checkpointed folds: pooled held-out operating point (scaled LR),
final scanner model, statsmodels coefficient inference. Folds CSV is authoritative."""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from src.morphoscan import metrics, features as _f

CACHE = "features_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
folds = pd.read_csv("results/trackB_folds_partial.csv")
out = {"folds": folds.to_dict("records")}
aucs = folds["auc_lr"].astype(float).values; aucsg = folds["auc_gb"].astype(float).values
out["median_auc_lr"] = float(np.nanmedian(aucs)); out["median_auc_gb"] = float(np.nanmedian(aucsg))
out["mean_auc_lr"] = float(np.nanmean(aucs)); out["mean_auc_gb"] = float(np.nanmean(aucsg))

meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
Xs, ys, gs = [], [], []
for _, r in meta.iterrows():
    p = os.path.join(CACHE, f"{r['section']}.npz")
    if not os.path.exists(p): continue
    z = np.load(p, allow_pickle=True)
    X = z["X"]; ok = np.isfinite(X).all(axis=1)
    Xs.append(X[ok]); ys.append(z["labels"][ok]); gs += [r["patient"]] * int(ok.sum())
X = np.vstack(Xs); y = np.concatenate(ys); groups = np.array(gs)
print("assembled", X.shape)

probs = np.full(len(y), np.nan)
for pat in sorted(set(groups)):
    tr = groups != pat; te = ~tr
    sc = StandardScaler().fit(X[tr])
    lr = LogisticRegression(max_iter=1000, class_weight="balanced").fit(sc.transform(X[tr]), y[tr])
    probs[te] = lr.predict_proba(sc.transform(X[te]))[:, 1]
    print("pooled", pat)
pred = (probs >= 0.5).astype(int)
tp = int(((pred == 1) & (y == 1)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
fn = int(((pred == 0) & (y == 1)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
out["operating_point"] = dict(sensitivity=tp / max(1, tp + fn), specificity=tn / max(1, tn + fp),
                              brier=metrics.brier_score(probs, y), note="scaled LR pooled held-out")

mu, sd = X.mean(0), X.std(0); sd[sd == 0] = 1
lr = LogisticRegression(max_iter=3000, class_weight="balanced").fit((X - mu) / sd, y)
np.savez("results/scanner_model.npz", coef=lr.coef_.ravel(), intercept=lr.intercept_, mu=mu, sd=sd)

import statsmodels.api as sm
smres = sm.Logit(y, sm.add_constant((X - mu) / sd)).fit(disp=0, maxiter=200)
names = _f.feature_names()
import csv as _c
with open("results/scanner_coefficients.csv", "w", newline="") as fh:
    w = _c.writer(fh); w.writerow(["feature", "coef", "z", "p"])
    for j in np.argsort(-np.abs(smres.params[1:])):
        w.writerow([names[j], float(smres.params[1:][j]), float(smres.tvalues[1:][j]), float(smres.pvalues[1:][j])])

json.dump(out, open("results/trackB_results.json", "w"), indent=1)
print("TRACKB DONE", out["median_auc_lr"], out["median_auc_gb"], out["operating_point"])
