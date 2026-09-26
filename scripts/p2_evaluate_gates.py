"""Phase-2 locked-gate evaluation from checkpoint partials (PREREG-2).
Reads results/p2_*_partial.* and Phase-1 results/trackA_per_gene.csv (M1).
Prints a verdict table; writes results/p2_gate_verdicts.json.
Honesty rules: numbers come only from these files, caveats attached.
"""
import json, csv, os
import numpy as np
from scipy import stats

R = "results"
def rows(p):
    with open(os.path.join(R, p)) as f:
        return list(csv.DictReader(f))

folds = rows("p2_folds_partial.csv")
pg = rows("p2_per_gene_partial.csv")
pp = rows("p2_per_pathway_partial.csv")
perm = json.load(open(os.path.join(R, "p2_perm_partial.json")))
scan = rows("p2_scanner_partial.csv")
p1 = rows("trackA_per_gene.csv")

def f(x):
    try: return float(x)
    except: return np.nan

# ---- G2-R: median held-out per-gene Pearson >= 0.08; beats M1 in >=80% of genes
gene_r = {}
for r in pg:
    gene_r.setdefault(r["gene"], []).append(f(r["r_m1"]))
p2_med = {g: float(np.nanmedian(v)) for g, v in gene_r.items()}
med_gene = float(np.nanmedian(list(p2_med.values())))
# Phase-1 M1 per-gene (trackA_per_gene.csv columns? inspect header at runtime)
p1cols = list(p1[0].keys())
gcol = next(c for c in p1cols if "gene" in c.lower())
m1col = next((c for c in p1cols if c.lower() in ("r_m1","m1","pearson_m1","m1_r")), None)
p1_med = {}
for r in p1:
    p1_med.setdefault(r[gcol], []).append(f(r[m1col]))
p1_med = {g: float(np.nanmedian(v)) for g, v in p1_med.items()}
common = sorted(set(p2_med) & set(p1_med))
beat = [g for g in common if p2_med[g] > p1_med[g]]
frac_beat = len(beat) / max(len(common), 1)
d = np.array([p2_med[g] - p1_med[g] for g in common])
wil_p = float(stats.wilcoxon(d).pvalue) if len(d) >= 10 else np.nan
g2r_pass = med_gene >= 0.08 and frac_beat >= 0.80 and wil_p < 0.05

# ---- G2-P: median pathway Pearson >= 0.20; perm p<0.01; immune AUROC >= 0.75
pw_r = {}
for r in pp:
    if r["r"] != "": pw_r.setdefault(r["pathway"], []).append(f(r["r"]))
pw_med = {k: float(np.nanmedian(v)) for k, v in pw_r.items()}
med_pw = float(np.nanmedian(list(pw_med.values())))
perm_pw = [f(r["pathway_median"]) for r in perm]
real_fold_pw = [f(r["median_r_pathway"]) for r in folds]
real_pw = float(np.nanmedian(real_fold_pw))
perm_p = (1 + sum(1 for x in perm_pw if x >= real_pw)) / (1 + len(perm_pw))
imm = [f(r["immune_auc_median"]) for r in folds if r["immune_auc_median"]]
imm_med = float(np.nanmedian(imm)) if imm else np.nan
g2p_pass = med_pw >= 0.20 and perm_p < 0.01 and imm_med >= 0.75

# ---- G2-B: vs ST-Net 0.19 (characterization per round-2 pre-committed caveat)
g2b_gap = med_gene - 0.19

# ---- G2-A: numeric null |median r| < 0.02
perm_gene = [f(r["gene_median"]) for r in perm]
g2a_gene = float(np.nanmedian(perm_gene))
g2a_pw = float(np.nanmedian(perm_pw))
g2a_pass = abs(g2a_gene) < 0.02 and abs(g2a_pw) < 0.02

# ---- G2-S: scanner patient-held-out AUC >= 0.85 (median per-patient AUC)
scan_auc = [f(r["auc"]) for r in scan if r["auc"]]
scan_med = float(np.nanmedian(scan_auc))
g2s_pass = scan_med >= 0.85

out = {
 "n_folds": len(folds),
 "G2-R": {"median_per_gene_r": med_gene, "frac_genes_beat_M1": frac_beat,
          "n_genes_common": len(common), "wilcoxon_p": wil_p,
          "pass": bool(g2r_pass),
          "caveat": "Phase-2 top-250 are fold-varying; comparison on genes present in both runs"},
 "G2-P": {"median_pathway_r": med_pw, "perm_p": perm_p, "n_perm": len(perm_pw),
          "immune_auc_median": imm_med, "pass": bool(g2p_pass)},
 "G2-B": {"median_per_gene_r": med_gene, "stnet": 0.19, "gap": g2b_gap,
          "note": "round-2 pre-committed: benchmark characterization, not pass/fail"},
 "G2-A": {"perm_gene_median": g2a_gene, "perm_pathway_median": g2a_pw,
          "pass": bool(g2a_pass), "threshold": "|r| < 0.02"},
 "G2-S": {"scanner_auc_median": scan_med, "pass": bool(g2s_pass)},
 "D1": {"genes_ge_015": int(sum(1 for v in p2_med.values() if v >= 0.15)),
        "note": "enrichment (confirmatory antigen-processing + exploratory) evaluated separately"},
}
json.dump(out, open(os.path.join(R, "p2_gate_verdicts.json"), "w"), indent=1)
for k, v in out.items():
    print(k, json.dumps(v) if isinstance(v, dict) else v)
