"""Merge enet shard fragments into results/p2_multitask.json and recompute the
summary blocks exactly as scripts/p2_multitask.py lines 108-122. Idempotent:
folds already carrying elasticnet_median_r are left untouched.
Run: python scripts/p2_enet_merge.py"""
import json, os, glob
import numpy as np

OUT = "results/p2_multitask.json"
state = json.load(open(OUT))
inserted, skipped = [], []
for frag in sorted(glob.glob("results/enet_shards/*.json")):
    rec = json.load(open(frag))
    if rec.get("skipped"):
        skipped.append(rec["patient"] + " (fragment says already-in-main)")
        continue
    for f in state["folds"]:
        if f["patient"] == rec["patient"]:
            if "elasticnet_median_r" in f:
                skipped.append(rec["patient"])
            else:
                f["elasticnet_median_r"] = rec["elasticnet_median_r"]
                f["elasticnet_n_iter"] = rec["elasticnet_n_iter"]
                inserted.append(rec["patient"])
            break
fs = state["folds"]
summ = state.get("summary", {})
enet_done = [f for f in fs if "elasticnet_median_r" in f]
for m in ["ridge_m1_committed", "ridge_m1_refold", "pls_median_r", "rrr_median_r"]:
    v = [f[m] for f in fs if m in f]
    if v:
        summ[m] = {"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v))}
if enet_done:
    v = [f["elasticnet_median_r"] for f in enet_done]
    summ["elasticnet_median_r"] = {"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v)),
                                   "n_folds": len(enet_done)}
wins = state.get("folds_better_than_ridge", {})
for m in ["pls_median_r", "rrr_median_r", "ridge_m1_refold"]:
    wins[m] = sum(1 for f in fs if m in f and f[m] > f["ridge_m1_committed"])
if enet_done:
    wins["elasticnet_median_r"] = sum(1 for f in enet_done if f["elasticnet_median_r"] > f["ridge_m1_committed"])
state["summary"] = summ
state["folds_better_than_ridge"] = wins
state["enet_folds_done"] = len(enet_done)
json.dump(state, open(OUT, "w"))
print(f"inserted {inserted} | skipped {skipped} | enet_folds_done {len(enet_done)}/10")
