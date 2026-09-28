"""Phase-2 figure regeneration (post-rebuild values). Overwrites fig1/2/3/5/6;
fig4_transport stays Phase-1 (labeled in text). Sources: trackA_results.json,
p2_per_gene_partial.csv, p2_perm10k.json + p2_perm_partial.json,
p2_scanner_partial.csv, p2_scanner_labelperm.json, p2_perm10k_verdict.json."""
import json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "results/"; F = "figures/"
BLUE, GRAY, RED = "#1f4e91", "#888888", "#b03030"

# fig1: per-fold medians
d = json.load(open(R+"trackA_results.json"))["folds"]
pats = [f["patient"].replace("BC","") for f in d]
m1 = [f["median_r_m1"] for f in d]; m2 = [f["median_r_m2"] for f in d]
x = np.arange(len(d))
fig, ax = plt.subplots(figsize=(8,3.2))
ax.plot(x, m1, "o-", color=GRAY, label="M1 ridge (embeddings)")
ax.plot(x, m2, "s-", color=BLUE, label="M2 + spatial smoothing")
ax.axhline(0.19, ls="--", color=RED, lw=1, label="ST-Net published ~0.19")
ax.set_xticks(x); ax.set_xticklabels(pats, rotation=45, fontsize=8)
ax.set_ylabel("median per-gene r"); ax.set_xlabel("held-out patient")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(F+"fig1_folds.pdf"); plt.close(fig)

# fig2: per-gene distribution
rows = list(csv.DictReader(open(R+"p2_per_gene_partial.csv")))
g1, g2 = {}, {}
for r in rows:
    if r["r_m1"] and r["r_m2"]:
        g1.setdefault(r["gene"], []).append(float(r["r_m1"]))
        g2.setdefault(r["gene"], []).append(float(r["r_m2"]))
v1 = np.array([np.median(v) for v in g1.values()]); v2 = np.array([np.median(v) for v in g2.values()])
fig, ax = plt.subplots(figsize=(5.2,3.4))
ax.hist(v1, bins=40, alpha=0.6, color=GRAY, label="M1")
ax.hist(v2, bins=40, alpha=0.6, color=BLUE, label="M2")
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("per-gene held-out r (median over sections)"); ax.set_ylabel("genes")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(F+"fig2_genes.pdf"); plt.close(fig)

# fig3: perm10k pathway null
app = json.load(open(R+"p2_perm10k.json"))["rows"]
null = [r["pathway_median"] for r in app]
try:
    com = json.load(open(R+"p2_perm_partial.json"))
    rows200 = com["rows"] if isinstance(com, dict) and "rows" in com else com
    null = [r["pathway_median"] if isinstance(r, dict) else r for r in rows200] + null
except Exception as e:
    print("committed-200 skipped:", e)
verd = json.load(open(R+"p2_perm10k_verdict.json"))
obs = verd["real_median_pathway_r"]
fig, ax = plt.subplots(figsize=(5.2,3.4))
ax.hist(null, bins=80, color=GRAY, alpha=0.8, label=f"null (n={len(null):,})")
ax.axvline(obs, color=RED, lw=1.6, label=f"observed {obs:.3f} (p={verd['perm_p_10k']:.4f})")
ax.set_xlabel("median pathway r under gene-axis permutation"); ax.set_ylabel("permutations")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(F+"fig3_permutation.pdf"); plt.close(fig)

# fig5: scanner AUCs + label-perm null
sc = [r for r in csv.DictReader(open(R+"p2_scanner_partial.csv")) if r.get("auc")]
aucs = sorted(float(r["auc"]) for r in sc)
lp = json.load(open(R+"p2_scanner_labelperm.json"))
nulls = [a for f in lp["folds"] for s in f["sections"] for a in s["perm_aucs"]]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9,3.4))
a1.bar(range(len(aucs)), aucs, color=BLUE)
a1.axhline(0.85, ls="--", color=RED, lw=1, label="G2-S gate 0.85")
a1.axhline(0.5, color="k", lw=0.7)
a1.set_xlabel(f"section (n={len(aucs)}, sorted)"); a1.set_ylabel("held-out AUC"); a1.legend(fontsize=8)
a2.hist(nulls, bins=40, color=GRAY, alpha=0.8, label=f"label-perm null (n={len(nulls)})")
a2.axvline(0.926, color=RED, lw=1.6, label="real median 0.926")
a2.set_xlabel("AUC"); a2.legend(fontsize=8)
fig.tight_layout(); fig.savefig(F+"fig5_scanner_auc.pdf"); plt.close(fig)

# fig6: M2 vs M1 scatter
fig, ax = plt.subplots(figsize=(4.6,4.2))
ax.scatter(v1, v2, s=6, alpha=0.4, color=BLUE)
lo, hi = min(v1.min(), v2.min()), max(v1.max(), v2.max())
ax.plot([lo,hi],[lo,hi], color="k", lw=0.8)
ax.set_xlabel("M1 per-gene r"); ax.set_ylabel("M2 per-gene r")
fig.tight_layout(); fig.savefig(F+"fig6_scatter.pdf"); plt.close(fig)
print("figures regenerated:", len(null), "null perms,", len(aucs), "scanner sections")
