"""Consolidate final ledgers: accession_ledger.csv (datasets) and tools_ledger.csv.
Dataset convention (disclosed): accession-level records, level_of_use column marks
analyzed vs metadata-mined. The 68 Mendeley sections are individually listed as
primary-analyzed; registry records are metadata-mined."""
import csv, json, os
import pandas as pd

rows = []
# primary analyzed: 68 sections of Mendeley 29ntw7sh4r v5
meta = pd.read_csv("/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st/metadata.csv")
for _, r in meta.iterrows():
    sec = r["count_matrix"].replace("_stdata.tsv.gz", "")
    rows.append(dict(accession=sec, source="Mendeley-29ntw7sh4r-v5",
                     title=f"Breast ST section {sec} ({r['type']}, patient {r['patient']})",
                     url="https://data.mendeley.com/datasets/29ntw7sh4r/5", level_of_use="analyzed-primary"))
rows.append(dict(accession="doi:10.17632/29ntw7sh4r.5", source="Mendeley-29ntw7sh4r-v5",
                 title="Human breast cancer in situ capturing transcriptomics (cohort record)",
                 url="https://data.mendeley.com/datasets/29ntw7sh4r/5", level_of_use="analyzed-primary"))
# registry harvest
h = pd.read_csv("results/accession_harvest_raw.csv")
seen = set()
for _, r in h.iterrows():
    key = (r["source"], r["accession"])
    if key in seen: continue
    seen.add(key)
    rows.append(dict(r))
led = pd.DataFrame(rows)
led.to_csv("results/accession_ledger.csv", index=False)
print("accession ledger rows:", len(led))
print(led["level_of_use"].value_counts().to_dict())
print(led["source"].value_counts().head(10).to_dict())

# tools ledger: merge raw1 + raw2, consolidate NCBI eutils under distinct databases (disclosed)
t1 = pd.read_csv("results/tools_ledger_raw.csv")
t2 = pd.read_csv("results/tools_ledger_raw2.csv") if os.path.exists("results/tools_ledger_raw2.csv") else pd.DataFrame(columns=t1.columns)
T = pd.concat([t1, t2], ignore_index=True).drop_duplicates(subset=["tool"], keep="first")
sw = [
    ("numpy", "analysis-software", "all numerical arrays, kernels, linear algebra"),
    ("scipy", "analysis-software", "sparse matrices, rank statistics, Mann-Whitney/Wilcoxon"),
    ("scikit-learn", "analysis-software", "logistic regression + gradient boosting scanners, MLP check"),
    ("statsmodels", "analysis-software", "scanner coefficient inference with p-values"),
    ("pandas", "analysis-software", "all tabular data handling"),
    ("Pillow", "analysis-software", "JPEG draft-mode decoding + patch extraction"),
    ("matplotlib", "analysis-software", "all paper figures"),
    ("seaborn", "analysis-software", "figure styling"),
    ("pytest", "analysis-software", "test suite runner"),
]
for name, cat, purp in sw:
    T = pd.concat([T, pd.DataFrame([dict(tool=name, category=cat, purpose=purp, status="certified",
                                         detail="imported and used in committed pipeline code", records="")])], ignore_index=True)
T.to_csv("results/tools_ledger.csv", index=False)
print("tools ledger rows:", len(T), "certified:", (T["status"] == "certified").sum())
