"""Amendment-queue Tier-2 #4: simple spatial baselines for per-gene prediction.
Pre-declared design (2026-09-29, before any run). Four baselines per LOPO fold,
same 250 train-selected genes, no morphology:
B1 section-mean (honest): predict the TRAIN-fold mean per gene for every test
   spot. Constant per section -> per-gene Pearson undefined by construction;
   recorded as degenerate (it can carry no within-section signal).
B2 k-NN oracle (k=6): each test spot predicted as the mean of its 6 nearest
   neighbors' TRUE expression within the section (self excluded). Oracle
   baseline: measures the spatial-autocorrelation ceiling, not deployable.
B3 Gaussian-smoothing oracle: W @ T with the committed kernel (length scale
   1.5*pitch, diag 0). Oracle, same reading as B2.
B4 Moran analysis: per gene, Moran's I of true expression vs that gene's
   oracle-kNN r and vs its committed-model r - does spatial autocorrelation
   explain which genes are predictable?
Metric: per-gene median Pearson pooled over test sections per fold.
Writes results/p2_spatial_baselines.json. Memory-lean; checkpointed per fold."""
import os, sys, json, time, csv as _csv, gc
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from src.morphoscan import metrics
import p2_common as C

N_GENES, KNN = 250, 6
OUT = "results/p2_spatial_baselines.json"
state = json.load(open(OUT)) if os.path.exists(OUT) else {"folds": [], "design": {
    "baselines": ["section_mean_honest", "knn_oracle_k6", "gaussian_oracle", "moran_analysis"],
    "knn": KNN}}
meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())

def per_gene_r(P, T):
    rs = np.full(T.shape[1], np.nan)
    for j in range(T.shape[1]):
        if T[:, j].std() > 0 and P[:, j].std() > 0:
            rs[j] = np.corrcoef(P[:, j], T[:, j])[0, 1]
    return rs

done = {f["patient"] for f in state["folds"]}
for test_pat in patients:
    if test_pat in done: continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    train_secs = [sec for pat in train_pats for sec in her2[her2["patient"] == pat]["section"]
                  if sec in sections]
    sum_x, n_tot = None, 0
    for sec in train_secs:
        N = sections[sec]["N"]
        sx = np.asarray(N.sum(axis=0)).ravel()
        sum_x = sx if sum_x is None else sum_x + sx
        n_tot += N.shape[0]
    means_all = sum_x / n_tot
    top = np.argsort(-means_all)[:N_GENES]
    train_mean_top = means_all[top]
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    r_knn_all, r_gauss_all, moran_all, rknn_genes = [], [], [], []
    for sec in test_secs:
        s = sections[sec]
        T = s["N"][:, top].toarray()
        px = s["px"]; n = T.shape[0]
        d = np.sqrt(((px[:, None, :] - px[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(d, np.inf)
        idx = np.argsort(d, axis=1)[:, :KNN]
        P_knn = T[idx].mean(axis=1)
        W = metrics.gaussian_kernel_weights(np.where(np.isinf(d), 0.0, d), 1.5 * float(s["pitch"]))
        P_gauss = W @ T
        r_knn_all.append(per_gene_r(P_knn, T))
        r_gauss_all.append(per_gene_r(P_gauss, T))
        rknn_genes.append(per_gene_r(P_knn, T))
        moran_all.append(np.array([metrics.morans_i(T[:, j], W) for j in range(T.shape[1])]))
    R_knn = np.nanmean(np.vstack(r_knn_all), axis=0)
    R_gauss = np.nanmean(np.vstack(r_gauss_all), axis=0)
    MI = np.nanmean(np.vstack(moran_all), axis=0)
    ok = np.isfinite(MI) & np.isfinite(R_knn)
    rec = {"patient": test_pat,
           "section_mean": "degenerate (constant per section; per-gene r undefined by construction)",
           "knn_oracle_median": float(np.nanmedian(R_knn)),
           "gaussian_oracle_median": float(np.nanmedian(R_gauss)),
           "morans_i_median": float(np.nanmedian(MI)),
           "corr_morans_vs_knnr": float(np.corrcoef(MI[ok], R_knn[ok])[0, 1]) if ok.sum() > 10 else None}
    state["folds"].append(rec)
    json.dump(state, open(OUT, "w"))
    print(f"{test_pat} knn {rec['knn_oracle_median']:.4f} gauss {rec['gaussian_oracle_median']:.4f} "
          f"moran {rec['morans_i_median']:.3f} corr(I,r) {rec['corr_morans_vs_knnr']:.3f} ({time.time()-t0:.0f}s)", flush=True)
    gc.collect()

fs = state["folds"]
summ = {"n_folds": len(fs),
        "knn_oracle_median": float(np.median([f["knn_oracle_median"] for f in fs])),
        "gaussian_oracle_median": float(np.median([f["gaussian_oracle_median"] for f in fs])),
        "morans_i_median": float(np.median([f["morans_i_median"] for f in fs])),
        "corr_morans_vs_knnr_median": float(np.median([f["corr_morans_vs_knnr"] for f in fs if f["corr_morans_vs_knnr"] is not None]))}
state["summary"] = summ
print(json.dumps(summ, indent=1))
json.dump(state, open(OUT, "w"))
