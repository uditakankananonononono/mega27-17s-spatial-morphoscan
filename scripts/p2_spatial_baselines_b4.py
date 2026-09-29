"""Addendum to #4 B4: corr(Moran's I, committed-model r_m2) per fold.
Joins per-(patient, section, gene) r_m2 from results/p2_per_gene_partial.csv
with freshly recomputed per-gene Moran's I (same recipe as
p2_spatial_baselines.py). Question: does the morphology model's per-gene
predictability track spatial autocorrelation the way the oracle's does (0.957)?
Writes the merged leg into results/p2_spatial_baselines.json as b4_model_leg."""
import os, sys, json, csv, collections
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from src.morphoscan import metrics
import p2_common as C

meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())

# model r_m2 per (patient, gene), mean over that patient's sections
model_r = collections.defaultdict(lambda: collections.defaultdict(list))
for r in csv.DictReader(open("results/p2_per_gene_partial.csv")):
    if r["r_m2"]:
        model_r[r["patient"]][r["gene"]].append(float(r["r_m2"]))

state = json.load(open("results/p2_spatial_baselines.json"))
leg = {}
for test_pat in patients:
    train_pats = [p for p in patients if p != test_pat]
    sum_x, n_tot = None, 0
    for pat in train_pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            N = sections[sec]["N"]
            sx = np.asarray(N.sum(axis=0)).ravel()
            sum_x = sx if sum_x is None else sum_x + sx
            n_tot += N.shape[0]
    top = np.argsort(-(sum_x / n_tot))[:250]
    genes = [universe[i] for i in top]
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    mi_per_gene = collections.defaultdict(list)
    for sec in test_secs:
        s = sections[sec]
        T = s["N"][:, top].toarray()
        px = s["px"]
        d = np.sqrt(((px[:, None, :] - px[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(d, np.inf)
        W = metrics.gaussian_kernel_weights(np.where(np.isinf(d), 0.0, d), 1.5 * float(s["pitch"]))
        for j, g in enumerate(genes):
            mi_per_gene[g].append(metrics.morans_i(T[:, j], W))
    mis, mrs = [], []
    for g in genes:
        if g in model_r[test_pat] and mi_per_gene[g]:
            mis.append(float(np.mean(mi_per_gene[g])))
            mrs.append(float(np.mean(model_r[test_pat][g])))
    mis, mrs = np.array(mis), np.array(mrs)
    c = float(np.corrcoef(mis, mrs)[0, 1]) if len(mis) > 10 else None
    leg[test_pat] = {"n_genes": len(mis), "corr_morans_vs_model_r2": c}
    print(f"{test_pat} n {len(mis)} corr(I, model_r2) {c:.3f}", flush=True)

state["b4_model_leg"] = leg
vals = [v["corr_morans_vs_model_r2"] for v in leg.values() if v["corr_morans_vs_model_r2"] is not None]
state["summary"]["corr_morans_vs_model_r2_median"] = float(np.median(vals))
json.dump(state, open("results/p2_spatial_baselines.json", "w"))
print(json.dumps({"corr_morans_vs_model_r2_median": float(np.median(vals))}))
