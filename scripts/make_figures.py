"""Generate all paper figures from committed results."""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
plt.rcParams.update({"font.size": 9, "figure.dpi": 150})

FIG = "figures"; os.makedirs(FIG, exist_ok=True)
A = json.load(open("results/trackA_results.json"))
B = json.load(open("results/trackB_results.json"))
G = pd.read_csv("results/trackA_per_gene.csv")

# F1: per-fold median r M0/M1/M2
F = pd.DataFrame(A["folds"])
fig, ax = plt.subplots(figsize=(6.2, 3.2))
x = np.arange(len(F)); w = 0.27
ax.bar(x - w, F["median_r_m0"], w, label="M0 train-mean", color="#bbbbbb")
ax.bar(x, F["median_r_m1"], w, label="M1 ridge probe", color="#2f6fb2")
ax.bar(x + w, F["median_r_m2"], w, label="M2 + spatial smoothing", color="#e0812f")
ax.axhline(0.19, ls="--", lw=1, color="k", label="ST-Net published ~0.19")
ax.set_xticks(x); ax.set_xticklabels(F["patient"], rotation=45, ha="right", fontsize=7)
ax.set_ylabel("median per-gene Pearson r"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{FIG}/fig1_folds.pdf"); plt.close(fig)

# F2: per-gene r distribution M1 vs M2
gg = G.groupby("gene")[["r_m1", "r_m2"]].median()
fig, ax = plt.subplots(figsize=(4.2, 3.2))
ax.hist(gg["r_m1"].dropna(), bins=40, alpha=0.6, label="M1", color="#2f6fb2")
ax.hist(gg["r_m2"].dropna(), bins=40, alpha=0.6, label="M2", color="#e0812f")
ax.axvline(gg["r_m1"].median(), color="#2f6fb2", ls="--", lw=1)
ax.axvline(gg["r_m2"].median(), color="#e0812f", ls="--", lw=1)
ax.set_xlabel("per-gene Pearson r (held-out patients)"); ax.set_ylabel("genes"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/fig2_genes.pdf"); plt.close(fig)

# F3: permutation null vs observed
fig, ax = plt.subplots(figsize=(4.2, 3.0))
ax.hist(A["perm_null_median_r"], bins=20, color="#999999")
ax.axvline(A["median_r_m1"], color="r", lw=1.5, label=f"observed M1 ({A['median_r_m1']:.3f})")
ax.set_xlabel("median per-gene r"); ax.set_ylabel("permutations"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/fig3_permutation.pdf"); plt.close(fig)

# F4: transport by subtype
T = pd.DataFrame(A["transport"])
fig, ax = plt.subplots(figsize=(5.0, 3.0))
sns.boxplot(data=T, x="subtype", y="median_r", ax=ax, color="#88b4de")
ax.axhline(A["median_r_m1"], ls="--", lw=1, color="k", label="HER2+ within-cohort median")
ax.set_ylabel("median per-gene r"); ax.set_xlabel("held-out subtype cohort")
ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(f"{FIG}/fig4_transport.pdf"); plt.close(fig)

# F5: Track B per-patient AUC
Fb = pd.DataFrame(B["folds"])
fig, ax = plt.subplots(figsize=(6.2, 3.0))
x = np.arange(len(Fb))
ax.bar(x - 0.2, Fb["auc_lr"], 0.4, label="logistic", color="#2f6fb2")
ax.bar(x + 0.2, Fb["auc_gb"], 0.4, label="gradient boosting", color="#57a35a")
ax.axhline(0.8, ls="--", lw=1, color="k", label="G-B1 gate 0.80")
ax.set_xticks(x); ax.set_xticklabels(Fb["patient"], rotation=45, ha="right", fontsize=6)
ax.set_ylabel("patient-held-out AUC"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{FIG}/fig5_scanner_auc.pdf"); plt.close(fig)

# F6: M1 vs M2 scatter per gene
fig, ax = plt.subplots(figsize=(3.6, 3.4))
ax.scatter(gg["r_m1"], gg["r_m2"], s=4, alpha=0.4, color="#444444")
lims = [min(gg.min()), max(gg.max())]
ax.plot(lims, lims, "r--", lw=1)
ax.set_xlabel("M1 per-gene r"); ax.set_ylabel("M2 per-gene r")
fig.tight_layout(); fig.savefig(f"{FIG}/fig6_scatter.pdf"); plt.close(fig)
print("figures done")
