"""Amendment-queue Tier-1 #1/#5: patient-level bootstrap 95% CIs for the Phase-2
headline metrics (G2-P median pathway r, per-gene median r m1/m2, immune AUC)
straight from the committed 10-patient fold table results/p2_folds_partial.csv.
10,000 nonparametric bootstrap resamples over patients; percentile CIs.
Light compute; deterministic seed. Writes results/p2_bootstrap_ci.json."""
import json, csv
import numpy as np

rows = list(csv.DictReader(open("results/p2_folds_partial.csv")))
assert len(rows) == 10, f"expected 10 committed folds, found {len(rows)}"
rng = np.random.default_rng(20260927)
metrics = ["median_r_m1", "median_r_m2", "median_r_pathway", "immune_auc_median"]
vals = {m: np.array([float(r[m]) for r in rows]) for m in metrics}
out = {"n_patients": len(rows), "n_boot": 10000, "seed": 20260927,
       "method": "patient-level nonparametric bootstrap, percentile 95% CI, statistic = median over patients",
       "metrics": {}}
for m in metrics:
    v = vals[m]
    idx = rng.integers(0, len(v), size=(10000, len(v)))
    boot = np.median(v[idx], axis=1)
    out["metrics"][m] = {
        "point_median": float(np.median(v)),
        "ci95_lo": float(np.percentile(boot, 2.5)),
        "ci95_hi": float(np.percentile(boot, 97.5)),
        "min": float(v.min()), "max": float(v.max())}
print(json.dumps(out, indent=2))
json.dump(out, open("results/p2_bootstrap_ci.json", "w"), indent=2)
