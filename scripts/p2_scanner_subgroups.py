"""Amendment-queue Tier-2 #14: scanner subgroup analyses (judge weakness 11).
Stratifies the committed scanner AUCs (results/p2_scanner_partial.csv, refreshed
to rebuilt-embedding canonical output) by: molecular subtype (metadata type),
tumor-fraction tertile, section-size tertile, and scanner-difficulty tertile
(|tumor_frac - 0.5|). Median AUC + IQR per stratum. Writes
results/p2_scanner_subgroups.json."""
import json
import numpy as np
import pandas as pd

scan = pd.read_csv("results/p2_scanner_partial.csv")
scan = scan[scan["auc"] != ""].copy()
scan["auc"] = scan["auc"].astype(float)
meta = pd.read_csv("data/her2st/metadata.csv")
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
df = scan.merge(meta[["section", "type"]], on="section", how="left")

def tertile(s):
    q = s.quantile([1/3, 2/3]).values
    return pd.cut(s, [-np.inf, q[0], q[1], np.inf], labels=["low", "mid", "high"])

df["frac_tertile"] = tertile(df["tumor_frac"])
df["size_tertile"] = tertile(df["n_spots"])
df["difficulty"] = (df["tumor_frac"] - 0.5).abs()
df["difficulty_tertile"] = tertile(df["difficulty"])

def summarize(col):
    out = {}
    for k, g in df.groupby(col, observed=True):
        out[str(k)] = {"n_sections": int(len(g)), "auc_median": float(g["auc"].median()),
                       "auc_iqr": [float(g["auc"].quantile(0.25)), float(g["auc"].quantile(0.75))],
                       "auc_min": float(g["auc"].min()), "auc_max": float(g["auc"].max())}
    return out

out = {"n_sections": int(len(df)),
       "by_subtype": summarize("type"),
       "by_tumor_fraction_tertile": summarize("frac_tertile"),
       "by_section_size_tertile": summarize("size_tertile"),
       "by_scanner_difficulty_tertile": summarize("difficulty_tertile"),
       "note": "difficulty = |tumor_frac - 0.5| (low = hard, near-balanced labels; high = easy, extreme fractions)"}
print(json.dumps(out, indent=1))
json.dump(out, open("results/p2_scanner_subgroups.json", "w"), indent=1)
