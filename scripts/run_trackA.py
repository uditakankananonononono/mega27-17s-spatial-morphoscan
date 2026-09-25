"""Track A: morphology -> spatial expression, leave-one-patient-out on HER2+ patients.
Locked protocol per PREREG.md. Writes results/trackA_results.json + per-gene table."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from src.morphoscan import models, metrics

CACHE = "features_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
N_GENES = 250
rng = np.random.default_rng(20260925)
log = open("/tmp/trackA.log", "w")

def say(m): print(m); log.write(m + "\n"); log.flush()

meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")

sections = {}
for sec in meta["section"]:
    p = os.path.join(CACHE, f"{sec}.npz")
    if not os.path.exists(p):
        continue
    z = np.load(p, allow_pickle=True)
    X = z["X"]; bad = ~np.isfinite(X).all(axis=1)
    counts = z["counts_sparse"].item().toarray() if hasattr(z["counts_sparse"], "item") else z["counts_sparse"]
    sections[sec] = dict(X=X[~bad], counts=counts[~bad], px=z["px"][~bad],
                         labels=z["labels"][~bad], genes=z["genes"], pitch=float(z["pitch"][0]))
say(f"loaded {len(sections)} sections")

her2 = meta[meta["is_her2"]]
patients = sorted(her2["patient"].unique())
say(f"HER2+ patients: {len(patients)}")

def assemble(pats):
    Xs, Cs, Pxs, Gs, Secs = [], [], [], [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Xs.append(s["X"]); Cs.append(s["counts"]); Pxs.append(s["px"]); Secs += [sec] * len(s["X"])
    X = np.vstack(Xs); C = np.vstack(Cs)
    return X, C, Pxs, Secs

fold_rows = []
per_gene_records = []
mlp_note = "not run"
for test_pat in patients:
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xtr, Ctr, _, SecTr = assemble(train_pats)
    test_secs = her2[her2["patient"] == test_pat]["section"]
    test_secs = [s for s in test_secs if s in sections]
    if not test_secs: continue
    Ntr = models.library_normalize_log1p(Ctr)
    genes = sections[test_secs[0]]["genes"]
    top = models.select_top_genes(Ntr, N_GENES)
    Ytr = Ntr[:, top]
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    lam = models.choose_lambda_inner(Ztr, Ytr, np.array([s.split("_")[0] for s in SecTr]))
    bank = models.fit_ridge_bank(Ztr, Ytr, lam)
    # alpha: inner LOPO predictions on train patients
    inner_pred, inner_true, inner_px = [], [], []
    for ip in train_pats:
        ip_train = [p for p in train_pats if p != ip]
        Xi, Ci, _, _ = assemble(ip_train)
        Ni = models.library_normalize_log1p(Ci)
        topi = models.select_top_genes(Ni, N_GENES)
        sci = models.Standardizer().fit(Xi)
        banki = models.fit_ridge_bank(sci.transform(Xi), Ni[:, topi], lam)
        for sec in her2[her2["patient"] == ip]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Pi = models.predict_bank(sci.transform(s["X"]), banki)
            Ni_sec = models.library_normalize_log1p(s["counts"])[:, topi]
            inner_pred.append(Pi); inner_true.append(Ni_sec); inner_px.append(s["px"])
    alpha = models.choose_alpha_inner(inner_pred, inner_true, inner_px, length_scale=1.5 * sections[test_secs[0]]["pitch"])
    # test predictions per section
    r1_all, r2_all, r0_all = [], [], []
    for sec in test_secs:
        s = sections[sec]
        Z = sc.transform(s["X"])
        P1 = models.predict_bank(Z, bank)
        T = models.library_normalize_log1p(s["counts"])[:, top]
        d = np.sqrt(((s["px"][:, None, :] - s["px"][None, :, :]) ** 2).sum(-1))
        W = metrics.gaussian_kernel_weights(d, 1.5 * s["pitch"])
        P2 = metrics.spatial_smooth(P1, W, alpha)
        P0 = np.tile(Ytr.mean(axis=0), (T.shape[0], 1))
        r1_all.append(models.per_gene_pearson(P1, T))
        r2_all.append(models.per_gene_pearson(P2, T))
        r0_all.append(models.per_gene_pearson(P0, T))
        for j, g in enumerate(genes[top]):
            per_gene_records.append(dict(patient=test_pat, section=sec, gene=g, r_m1=r1_all[-1][j], r_m2=r2_all[-1][j]))
    m1 = float(np.nanmedian(np.concatenate(r1_all))); m2 = float(np.nanmedian(np.concatenate(r2_all)))
    m0 = float(np.nanmedian(np.concatenate(r0_all)))
    fold_rows.append(dict(patient=test_pat, n_sections=len(test_secs), lam=lam, alpha=alpha,
                          median_r_m0=m0, median_r_m1=m1, median_r_m2=m2))
    say(f"fold {test_pat}: secs={len(test_secs)} lam={lam} alpha={alpha} r0={m0:.3f} r1={m1:.3f} r2={m2:.3f} ({time.time()-t0:.0f}s)")

F = pd.DataFrame(fold_rows)
G = pd.DataFrame(per_gene_records)
res = {"folds": fold_rows,
       "median_r_m0": float(F["median_r_m0"].median()),
       "median_r_m1": float(F["median_r_m1"].median()),
       "median_r_m2": float(F["median_r_m2"].median()),
       "mean_r_m1": float(F["median_r_m1"].mean()), "mean_r_m2": float(F["median_r_m2"].mean())}

# G-A1: permutation test on M1 (permute gene columns jointly, refit with frozen per-fold lam)
say("permutation test (50x) ...")
null = []
for b in range(50):
    rs = []
    for fr in fold_rows:
        test_pat = fr["patient"]
        train_pats = [p for p in patients if p != test_pat]
        Xtr, Ctr, _, _ = assemble(train_pats)
        Ntr = models.library_normalize_log1p(Ctr)
        perm = rng.permutation(Ntr.shape[1])[:N_GENES]
        Ytr = Ntr[:, perm]
        sc = models.Standardizer().fit(Xtr)
        bank = models.fit_ridge_bank(sc.transform(Xtr), Ytr, fr["lam"])
        test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
        rrs = []
        genes = sections[test_secs[0]]["genes"]
        for sec in test_secs:
            s = sections[sec]
            P = models.predict_bank(sc.transform(s["X"]), bank)
            T = models.library_normalize_log1p(s["counts"])[:, perm]
            rrs.append(models.per_gene_pearson(P, T))
        rs.append(np.nanmedian(np.concatenate(rrs)))
    null.append(float(np.median(rs)))
res["perm_null_median_r"] = [float(x) for x in null]
res["perm_p_m1"] = metrics.permutation_pvalue(res["median_r_m1"], np.array(null))

# G-A2: paired Wilcoxon on per-gene r (M2 vs M1), per gene averaged over folds
gg = G.groupby("gene")[["r_m1", "r_m2"]].median()
res["wilcoxon_p_m2_vs_m1"] = metrics.wilcoxon_p(gg["r_m2"].values, gg["r_m1"].values)
res["median_gene_r_m1"] = float(gg["r_m1"].median()); res["median_gene_r_m2"] = float(gg["r_m2"].median())
res["n_genes_eval"] = int(len(gg))

# G-A4 transport: train all HER2+, test each non-HER2 patient
say("transport to non-HER2 ...")
Xh, Ch, _, SecH = assemble(patients)
Nh = models.library_normalize_log1p(Ch)
toph = models.select_top_genes(Nh, N_GENES)
sch = models.Standardizer().fit(Xh)
lamh = models.choose_lambda_inner(sch.transform(Xh), Nh[:, toph], np.array([s.split("_")[0] for s in SecH]))
bankh = models.fit_ridge_bank(sch.transform(Xh), Nh[:, toph], lamh)
transport = []
for _, row in meta[~meta["is_her2"]].iterrows():
    sec = row["section"]
    if sec not in sections: continue
    s = sections[sec]
    P = models.predict_bank(sch.transform(s["X"]), bankh)
    T = models.library_normalize_log1p(s["counts"])[:, toph]
    r = models.per_gene_pearson(P, T)
    transport.append(dict(section=sec, subtype=row["type"], patient=row["patient"], median_r=float(np.nanmedian(r))))
res["transport"] = transport
res["genes_top"] = [str(g) for g in genes[toph]]

G.to_csv("results/trackA_per_gene.csv", index=False)
json.dump(res, open("results/trackA_results.json", "w"), indent=1)
say("TRACKA DONE " + json.dumps({k: v for k, v in res.items() if not isinstance(v, (list, dict))}))
log.close()
