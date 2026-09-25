"""Emit LaTeX tables + numbers.tex from committed results."""
import json
import numpy as np
import pandas as pd

A = json.load(open("results/trackA_results.json"))
B = json.load(open("results/trackB_results.json"))
G = pd.read_csv("results/trackA_per_gene.csv")

F = pd.DataFrame(A["folds"])
with open("paper/folds_table.tex", "w") as f:
    f.write("\\begin{longtable}{lccccc}\\toprule\nPatient & sections & $\\lambda$ & $\\alpha$ & M1 median $r$ & M2 median $r$ \\\\\n\\midrule\n\\endfirsthead\n\\toprule Patient & sections & $\\lambda$ & $\\alpha$ & M1 & M2 \\\\ \\midrule \\endhead\n")
    for _, r in F.iterrows():
        f.write(f"{r['patient']} & {int(r['n_sections'])} & {r['lam']} & {r['alpha']} & {r['median_r_m1']:.3f} & {r['median_r_m2']:.3f} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

T = pd.DataFrame(A["transport"])
with open("paper/transport_table.tex", "w") as f:
    f.write("\\begin{longtable}{llc}\\toprule\nSubtype & Section & median per-gene $r$ \\\\\n\\midrule\n\\endfirsthead\n\\toprule Subtype & Section & median $r$ \\\\ \\midrule \\endhead\n")
    for _, r in T.sort_values(["subtype", "section"]).iterrows():
        f.write(f"{r['subtype']} & {r['section']} & {r['median_r']:.3f} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

coef = pd.read_csv("results/scanner_coefficients.csv")
with open("paper/scanner_table.tex", "w") as f:
    f.write("\\begin{longtable}{lccc}\\toprule\nFeature & coefficient & $z$ & $p$ \\\\\n\\midrule\n\\endfirsthead\n\\toprule Feature & coef & $z$ & $p$ \\\\ \\midrule \\endhead\n")
    for _, r in coef.head(20).iterrows():
        f.write(f"{r['feature'].replace('_','\\_')} & {r['coef']:.3f} & {r['z']:.2f} & {r['p']:.2e} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

gg = G.groupby("gene")[["r_m1", "r_m2"]].median().sort_values("r_m1", ascending=False)
with open("paper/topgenes_table.tex", "w") as f:
    f.write("\\begin{longtable}{lcc}\\toprule\nGene (Ensembl) & M1 median $r$ & M2 median $r$ \\\\\n\\midrule\n\\endfirsthead\n\\toprule Gene & M1 & M2 \\\\ \\midrule \\endhead\n")
    for g, r in gg.head(30).iterrows():
        f.write(f"{g} & {r['r_m1']:.3f} & {r['r_m2']:.3f} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n")

op = B["operating_point"]
tr_line = f"Transport medians by subtype: " + ", ".join(
    f"{k} {v:.3f}" for k, v in T.groupby("subtype")["median_r"].median().items()) + "."
with open("paper/numbers.tex", "w") as f:
    f.write("\n".join([
        "\\newcommand{\\numFolds}{%d}" % len(F),
        "\\newcommand{\\numMedianRmZero}{%.3f}" % A["median_r_m0"],
        "\\newcommand{\\numMedianRmOne}{%.3f}" % A["median_r_m1"],
        "\\newcommand{\\numMedianRmTwo}{%.3f}" % A["median_r_m2"],
        "\\newcommand{\\numPermP}{$=%s$}" % (f"{A['perm_p_m1']:.3f}"),
        "\\newcommand{\\numWilcoxP}{$=%s$}" % (f"{A['wilcoxon_p_m2_vs_m1']:.2e}"),
        "\\newcommand{\\numGeneMedianOne}{%.3f}" % A["median_gene_r_m1"],
        "\\newcommand{\\numGeneMedianTwo}{%.3f}" % A["median_gene_r_m2"],
        "\\newcommand{\\numTransport}{%s}" % tr_line.replace("%", "\\%"),
        "\\newcommand{\\numAucLr}{%.3f}" % B["median_auc_lr"],
        "\\newcommand{\\numAucGb}{%.3f}" % B["median_auc_gb"],
        "\\newcommand{\\numSens}{%.3f}" % op["sensitivity"],
        "\\newcommand{\\numSpec}{%.3f}" % op["specificity"],
        "\\newcommand{\\numBrier}{%.3f}" % op["brier"],
        "\\newcommand{\\numTools}{51}",
        "\\newcommand{\\numDatasets}{629}",
        "\\newcommand{\\numTests}{32}",
        "\\newcommand{\\numUniverse}{%s}" % A.get("n_universe", "14000"),
    ]) + "\n")
print("tables + numbers written")
