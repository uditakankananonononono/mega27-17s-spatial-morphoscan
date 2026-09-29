#!/usr/bin/env python3
"""Queue #10 closure: learned image baseline vs handcrafted, same patient-held-out split.

Reads two committed result artifacts (no new model fitting):
  - results/trackA_folds_partial.csv   : Track A M1 ridge on 41 handcrafted morphology
                                         features, corrected-geometry rerun (123b892),
                                         10-patient leave-one-patient-out, 250-gene panel.
  - results/p2_multitask.json          : ridge_m1_committed per fold = Phase-2 ridge on
                                         768-d CTransPath learned embeddings, same 10
                                         patients held out, fold-varying top-250 panel.
Both are median per-gene Pearson r on the held-out patient's spots.
Outputs results/learned_vs_handcrafted.json with the paired comparison.
"""
import json, csv
import numpy as np
from scipy import stats

ta = {r["patient"]: float(r["median_r_m1"])
      for r in csv.DictReader(open("results/trackA_folds_partial.csv"))}
p2 = {f["patient"]: f["ridge_m1_committed"]
      for f in json.load(open("results/p2_multitask.json"))["folds"]}

common = sorted(set(ta) & set(p2))
assert len(common) == 10, common
hc = np.array([ta[p] for p in common])   # handcrafted
le = np.array([p2[p] for p in common])   # learned embeddings
delta = le - hc
W, pval = stats.wilcoxon(le, hc)

out = {
    "description": "Queue #10: learned image baseline (CTransPath 768-d embeddings, Phase-2 ridge) "
                   "vs handcrafted morphology baseline (41 features, Track A M1 ridge), identical "
                   "10-patient leave-one-patient-out splits, median per-gene Pearson r on held-out spots.",
    "sources": {"handcrafted": "results/trackA_folds_partial.csv (corrected-geometry rerun 123b892)",
                "learned": "results/p2_multitask.json ridge_m1_committed (rebuilt canonical embeddings)"},
    "patients": common,
    "per_fold": [{"patient": p, "handcrafted_r": float(h), "learned_r": float(l), "delta": float(d)}
                 for p, h, l, d in zip(common, hc, le, delta)],
    "median_handcrafted": float(np.median(hc)),
    "median_learned": float(np.median(le)),
    "median_paired_delta": float(np.median(delta)),
    "folds_learned_better": int((delta > 0).sum()),
    "wilcoxon_W": float(W),
    "wilcoxon_p": float(pval),
    "gate_context": {"G2_R_gate": 0.08,
                     "both_below_gate": bool(np.median(hc) < 0.08 and np.median(le) < 0.08)},
    "caveats": ["Phase-2 gene panel is fold-varying top-250 by train-fold mean expression; "
                "Track A uses its committed 250-gene panel - panels overlap but are not identical.",
                "Elastic-net l1_ratio differs slightly (Track A alpha=0.4, Phase-2 alpha=0.5); "
                "both are ridge-dominant linear probes at lambda=100.",
                "No new model fitting in this script; it compares two committed, replay-verified artifacts."],
}
json.dump(out, open("results/learned_vs_handcrafted.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ["median_handcrafted","median_learned","median_paired_delta",
                                       "folds_learned_better","wilcoxon_p"]}, indent=1))
