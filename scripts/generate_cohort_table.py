"""Appendix table: one row per analyzed section (cohort manifest)."""
import os, json
import numpy as np, pandas as pd
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")
enr = json.load(open("results/enrichment_top25.json"))

with open("paper/cohort_table.tex", "w") as f:
    f.write("\\begin{longtable}{lllccc}\\toprule\nPatient & Section & Subtype & Spots & Tumor frac. & HER2+ \\\\\n\\midrule\n\\endfirsthead\n\\toprule Patient & Section & Subtype & Spots & Tumor & HER2+ \\\\ \\midrule \\endhead\n")
    for _, r in meta.iterrows():
        sec = r["section"]; p = f"features_cache/{sec}.npz"
        if not os.path.exists(p): continue
        z = np.load(p, allow_pickle=True)
        n = int(z["X"].shape[0]); lab = z["labels"]
        tf = float(np.mean([str(x).lower().startswith("t") for x in lab])) if len(lab) else 0.0
        her2 = "yes" if bool(r["is_her2"]) else "no"
        st = str(r["type"]).replace("_", "\\_"); sc = sec.replace("_", "\\_")
        f.write(f"{r['patient']} & {sc} & {st} & {n} & {tf:.2f} & {her2} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

with open("paper/enrichment_table.tex", "w") as f:
    f.write("\\begin{longtable}{llcc}\\toprule\nSource & Term & $p$ & Genes \\\\\n\\midrule\n\\endfirsthead\n\\toprule Source & Term & $p$ & Genes \\\\ \\midrule \\endhead\n")
    for t in enr["terms"][:15]:
        term = t["term"].replace("_", "\\_")
        f.write(f"{t['source']} & {term} & {t['p']:.2e} & {t['intersection']} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")
print("cohort + enrichment tables written")
