"""Generate scanner appendix LaTeX tables from committed results (paper growth).
Reads results/p2_scanner_partial.csv, p2_scanner_subgroups.json, p2_scanner_stability.json,
p2_scanner_labelperm.json, scanner_coefficients.csv. Writes paper/tab_scanner_*.tex."""
import json, csv

B = chr(92)           # backslash
RE = " " + B + B + chr(10)   # LaTeX row end + newline
def esc(s): return s.replace("_", B + "_")
def f3(x): return f"{x:.3f}"

# 1. per-section AUC table
rows = list(csv.DictReader(open("results/p2_scanner_partial.csv")))
hdr = B + "toprule Patient & Section & $n$ spots & Tumor frac. & AUC" + RE + B + "midrule"
with open("paper/tab_scanner_sections.tex", "w") as out:
    out.write(B + "begin{longtable}{llrrr}" + chr(10))
    out.write(B + "caption{Phase-2 scanner: per-section patient-held-out AUC. Sections with undefined AUC "
              "(tumor fraction 1.0: BC24105" + B + "_C1, BC24105" + B + "_C2) are excluded; BC24105" + B +
              "_D1 is retained and reported as a named outlier.}" + B + "label{tab:scannersec}" + RE)
    out.write(hdr + B + "endfirsthead" + chr(10) + hdr + B + "endhead" + chr(10))
    for r in sorted(rows, key=lambda r: (r["patient"], r["section"])):
        auc = r["auc"]
        aucs = f3(float(auc)) if auc not in ("", "nan", None) else "---"
        out.write(f"{esc(r['patient'])} & {esc(r['section'])} & {r['n_spots']} & {f3(float(r['tumor_frac']))} & {aucs}" + RE)
    out.write(B + "bottomrule" + chr(10) + B + "end{longtable}" + chr(10))

# 2. subgroup table
sg = json.load(open("results/p2_scanner_subgroups.json"))
with open("paper/tab_scanner_subgroups.tex", "w") as out:
    out.write(B + "begin{table}[h]" + B + "centering" + B + "footnotesize" + chr(10))
    out.write(B + "caption{Phase-2 scanner subgroup performance (27 evaluable sections).}" + B + "label{tab:scannersub}" + chr(10))
    out.write(B + "begin{tabular}{lrrrr}" + chr(10) + B + "toprule Subgroup & $n$ & Median AUC & IQR & Range" + RE + B + "midrule" + chr(10))
    blocks = [("Subtype", sg["by_subtype"]), ("Tumor-fraction tertile", sg["by_tumor_fraction_tertile"])]
    for gi, (grp, blk) in enumerate(blocks):
        for k, v in blk.items():
            label = k.replace("_", " ")
            out.write(f"{label} & {v['n_sections']} & {f3(v['auc_median'])} & "
                      f"[{f3(v['auc_iqr'][0])}, {f3(v['auc_iqr'][1])}] & [{f3(v['auc_min'])}, {f3(v['auc_max'])}]" + RE)
        if gi == 0: out.write(B + "midrule" + chr(10))
    out.write(B + "bottomrule" + chr(10) + B + "end{tabular}" + chr(10) + B + "end{table}" + chr(10))

# 3. stability table
st = json.load(open("results/p2_scanner_stability.json"))
boot = st.get("spearman_median_boot95")
boots = f"[{f3(boot[0])}, {f3(boot[1])}]" if boot else "---"
with open("paper/tab_scanner_stability.tex", "w") as out:
    out.write(B + "begin{table}[h]" + B + "centering" + B + "footnotesize" + chr(10))
    out.write(B + "caption{Phase-2 scanner coefficient stability across the 10 LOPO folds "
              "(768-dimensional embedding-logistic coefficient vectors).}" + B + "label{tab:scannerstab}" + chr(10))
    out.write(B + "begin{tabular}{lr}" + chr(10) + B + "toprule Metric & Value" + RE + B + "midrule" + chr(10))
    out.write(f"Pairwise Spearman, median & {f3(st['pairwise_spearman_median'])}" + RE)
    out.write(f"Pairwise Spearman, bootstrap 95\\% CI & {boots}" + RE)
    out.write(f"Pairwise Spearman, range & [{f3(st['pairwise_spearman_min'])}, {f3(st['pairwise_spearman_max'])}]" + RE)
    out.write(f"Top-50 feature Jaccard, median & {f3(st['top50_jaccard_median'])}" + RE)
    out.write(f"Top-50 sign consistency, median & {f3(st['top50_overlap_sign_consistency_median'])}" + RE)
    out.write(B + "bottomrule" + chr(10) + B + "end{tabular}" + chr(10) + B + "end{table}" + chr(10))

# 4. label-permutation per-section table
lp = json.load(open("results/p2_scanner_labelperm.json"))
sec_recs = []
for fold in lp["folds"]:
    for s in fold["sections"]:
        pa = sorted(s["perm_aucs"]); n = len(pa)
        med = pa[n // 2] if n % 2 else 0.5 * (pa[n // 2 - 1] + pa[n // 2])
        p95 = pa[int(0.95 * (n - 1))]
        exceed = sum(1 for a in pa if a >= s["true_auc"])
        sec_recs.append((s["section"], s["true_auc"], med, p95, exceed, n))
hdr4 = B + "toprule Section & True AUC & Null median & Null p95 & Exceedances" + RE + B + "midrule"
with open("paper/tab_scanner_labelperm.tex", "w") as out:
    out.write(B + "begin{longtable}{lrrrr}" + chr(10))
    out.write(B + "caption{Phase-2 scanner within-section label-permutation control (20 permutations per section): "
              "true AUC vs the permutation null. Exceedances counts permuted fits reaching the true AUC.}" +
              B + "label{tab:scannerperm}" + RE)
    out.write(hdr4 + B + "endfirsthead" + chr(10) + hdr4 + B + "endhead" + chr(10))
    for sec, t, med, p95, ex, n in sorted(sec_recs):
        out.write(f"{esc(sec)} & {f3(t)} & {f3(med)} & {f3(p95)} & {ex}/{n}" + RE)
    out.write(B + "bottomrule" + chr(10) + B + "end{longtable}" + chr(10))

print("wrote 4 tables:", len(rows), "sections,", len(sec_recs), "perm rows")
