"""Track A: morphology -> spatial expression, leave-one-patient-out on HER2+ patients.
Locked protocol per PREREG.md. Sections carry different Ensembl universes: counts are
aligned by gene-ID intersection and kept sparse (CSR) until a 250-gene slice is taken.
Incremental: resumes from results/trackA_folds_partial.csv / trackA_perm_partial.txt."""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from src.morphoscan import models, metrics

CACHE = "features_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
N_GENES = 250
rng = np.random.default_rng(20260925)
log = open("/tmp/trackA.log", "a")
def say(m): print(m); log.write(m + "\n"); log.flush()

meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")

universe = None
upos = None

def build_universe(secs):
    gene_sets = []
    for sec in secs:
        z = np.load(os.path.join(CACHE, f"{sec}.npz"), allow_pickle=True)
        gene_sets.append(set(str(g) for g in z["genes"]))
        z.close()
    return sorted(set.intersection(*gene_sets))

def load_sections(secs):
    """Load + normalize only the requested sections (memory-lean)."""
    out = {}
    for sec in secs:
        p = os.path.join(CACHE, f"{sec}.npz")
        if not os.path.exists(p): continue
        z = np.load(p, allow_pickle=True)
        X = z["X"]; ok = np.isfinite(X).all(axis=1)
        csr = z["counts_sparse"].item()[ok]
        genes = [str(g) for g in z["genes"]]
        px = z["px"][ok]; labels = z["labels"][ok]; pitch = float(z["pitch"][0])
        z.close()
        lib = np.asarray(csr.sum(axis=1)).ravel(); lib[lib == 0] = 1.0
        lut = np.array([upos.get(g, -1) for g in genes])
        m = csr.tocoo(); nl = lut[m.col]; keep = nl >= 0
        N = sparse.csr_matrix((m.data[keep].astype(np.float32), (m.row[keep], nl[keep])),
                              shape=(m.shape[0], len(upos)), dtype=np.float32)
        N.data *= np.repeat((1e4 / lib).astype(np.float32), np.diff(N.indptr))
        np.log1p(N.data, out=N.data)
        out[sec] = dict(X=X[ok], N=N, px=px, labels=labels, pitch=pitch)
        del csr, m, N
    import gc; gc.collect()
    return out

all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(CACHE, f"{s}.npz"))]
her2_secs = [s for s in meta[meta["is_her2"]]["section"] if s in all_secs]
universe = build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
say(f"gene universe (intersection): {len(universe)}")
sections = load_sections(her2_secs)
say(f"loaded {len(sections)} HER2+ sections")

her2 = meta[meta["is_her2"]]
patients = sorted(her2["patient"].unique())
say(f"HER2+ patients: {len(patients)}")

def assemble(pats):
    Xs, Ns = [], []
    secs_of = []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Xs.append(s["X"]); Ns.append(s["N"]); secs_of += [pat] * len(s["X"])
    return np.vstack(Xs), sparse.vstack(Ns), np.array(secs_of)

PARTIAL = "results/trackA_folds_partial.csv"
fold_rows = list(_csv.DictReader(open(PARTIAL))) if os.path.exists(PARTIAL) else []
done_pats = {r["patient"] for r in fold_rows}
PG = "results/trackA_per_gene_partial.csv"
per_gene_records = list(_csv.DictReader(open(PG))) if os.path.exists(PG) else []

for test_pat in patients:
    if test_pat in done_pats:
        say(f"skip {test_pat} (done)"); continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xtr, Ntr, groups_tr = assemble(train_pats)
    top = models.select_top_genes(Ntr, N_GENES)
    Ytr = Ntr[:, top].toarray()
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    lam = models.choose_lambda_inner(Ztr, Ytr, groups_tr)
    bank = models.fit_ridge_bank(Ztr, Ytr, lam)
    # alpha via inner LOPO on train patients
    inner_pred, inner_true, inner_px = [], [], []
    for ip in train_pats:
        Xi, Ni, _ = assemble([p for p in train_pats if p != ip])
        topi = models.select_top_genes(Ni, N_GENES)
        sci = models.Standardizer().fit(Xi)
        banki = models.fit_ridge_bank(sci.transform(Xi), Ni[:, topi].toarray(), lam)
        for sec in her2[her2["patient"] == ip]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            inner_pred.append(models.predict_bank(sci.transform(s["X"]), banki))
            inner_true.append(s["N"][:, topi].toarray())
            inner_px.append(s["px"])
    alpha = models.choose_alpha_inner(inner_pred, inner_true, inner_px,
                                      length_scale=1.5 * np.median([sections[s]["pitch"] for s in sections]))
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    r1_all, r2_all, r0_all = [], [], []
    for sec in test_secs:
        s = sections[sec]
        P1 = models.predict_bank(sc.transform(s["X"]), bank)
        T = s["N"][:, top].toarray()
        d = np.sqrt(((s["px"][:, None, :] - s["px"][None, :, :]) ** 2).sum(-1))
        W = metrics.gaussian_kernel_weights(d, 1.5 * s["pitch"])
        P2 = metrics.spatial_smooth(P1, W, alpha)
        P0 = np.tile(Ytr.mean(axis=0), (T.shape[0], 1))
        r1 = models.per_gene_pearson(P1, T); r2 = models.per_gene_pearson(P2, T)
        r1_all.append(r1); r2_all.append(r2); r0_all.append(np.nan_to_num(models.per_gene_pearson(P0, T), nan=0.0))
        for j in range(len(top)):
            per_gene_records.append(dict(patient=test_pat, section=sec, gene=universe[top[j]],
                                         r_m1=float(r1[j]), r_m2=float(r2[j])))
    m1 = float(np.nanmedian(np.concatenate(r1_all))); m2 = float(np.nanmedian(np.concatenate(r2_all)))
    m0 = float(np.nanmedian(np.concatenate(r0_all)))
    fold_rows.append(dict(patient=test_pat, n_sections=len(test_secs), lam=lam, alpha=alpha,
                          median_r_m0=m0, median_r_m1=m1, median_r_m2=m2))
    say(f"fold {test_pat}: secs={len(test_secs)} lam={lam} alpha={alpha} r0={m0:.3f} r1={m1:.3f} r2={m2:.3f} ({time.time()-t0:.0f}s)")
    with open(PARTIAL, "w", newline="") as _f:
        _w = _csv.DictWriter(_f, fieldnames=list(fold_rows[0].keys())); _w.writeheader(); _w.writerows(fold_rows)
    with open(PG, "w", newline="") as _f:
        _w = _csv.DictWriter(_f, fieldnames=list(per_gene_records[0].keys())); _w.writeheader(); _w.writerows(per_gene_records)

F = pd.DataFrame(fold_rows)
for _c in ("n_sections", "lam", "alpha", "median_r_m0", "median_r_m1", "median_r_m2"):
    F[_c] = pd.to_numeric(F[_c], errors="coerce")
G = pd.DataFrame(per_gene_records)
G["r_m1"] = pd.to_numeric(G["r_m1"], errors="coerce")
G["r_m2"] = pd.to_numeric(G["r_m2"], errors="coerce")
res = {"folds": fold_rows, "n_universe": len(universe),
       "median_r_m0": float(F["median_r_m0"].fillna(0.0).median()),
       "median_r_m1": float(F["median_r_m1"].median()),
       "median_r_m2": float(F["median_r_m2"].median()),
       "mean_r_m1": float(F["median_r_m1"].mean()), "mean_r_m2": float(F["median_r_m2"].mean())}

say("permutation test (50x) ...")
PERM = "results/trackA_perm_partial.txt"
null = [float(x) for x in open(PERM).read().split()] if os.path.exists(PERM) else []
for b in range(len(null), 50):
    rs = []
    for fr in fold_rows:
        test_pat = fr["patient"]
        train_pats = [p for p in patients if p != test_pat]
        Xtr, Ntr, _ = assemble(train_pats)
        perm = rng.permutation(Ntr.shape[1])[:N_GENES]
        Ytr = Ntr[:, perm].toarray()
        sc = models.Standardizer().fit(Xtr)
        bank = models.fit_ridge_bank(sc.transform(Xtr), Ytr, float(fr["lam"]))
        rrs = []
        for sec in [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]:
            s = sections[sec]
            P = models.predict_bank(sc.transform(s["X"]), bank)
            T = s["N"][:, perm].toarray()
            rrs.append(models.per_gene_pearson(P, T))
        rs.append(np.nanmedian(np.concatenate(rrs)))
    null.append(float(np.median(rs)))
    open(PERM, "w").write(" ".join(f"{x:.6f}" for x in null))
    if (b + 1) % 10 == 0: say(f"perm {b+1}/50")
res["perm_null_median_r"] = [float(x) for x in null]
res["perm_p_m1"] = metrics.permutation_pvalue(res["median_r_m1"], np.array(null))

gg = G.groupby("gene")[["r_m1", "r_m2"]].median()
res["wilcoxon_p_m2_vs_m1"] = metrics.wilcoxon_p(gg["r_m2"].values, gg["r_m1"].values)
res["median_gene_r_m1"] = float(gg["r_m1"].median()); res["median_gene_r_m2"] = float(gg["r_m2"].median())
res["n_genes_eval"] = int(len(gg))

say("transport to non-HER2 ...")
Xh, Nh, gh = assemble(patients)
toph = models.select_top_genes(Nh, N_GENES)
sch = models.Standardizer().fit(Xh)
lamh = models.choose_lambda_inner(sch.transform(Xh), Nh[:, toph].toarray(), gh)
bankh = models.fit_ridge_bank(sch.transform(Xh), Nh[:, toph].toarray(), lamh)
transport = []
nonher2_secs = [s for s in meta[~meta["is_her2"]]["section"] if s in all_secs]
sections = load_sections(nonher2_secs)
say(f"loaded {len(sections)} non-HER2 sections for transport")
for _, row in meta[~meta["is_her2"]].iterrows():
    sec = row["section"]
    if sec not in sections: continue
    s = sections[sec]
    P = models.predict_bank(sch.transform(s["X"]), bankh)
    T = s["N"][:, toph].toarray()
    r = models.per_gene_pearson(P, T)
    transport.append(dict(section=sec, subtype=row["type"], patient=row["patient"], median_r=float(np.nanmedian(r))))
res["transport"] = transport
res["genes_top"] = [universe[i] for i in toph]

G.to_csv("results/trackA_per_gene.csv", index=False)
json.dump(res, open("results/trackA_results.json", "w"), indent=1)
say("TRACKA DONE " + json.dumps({k: v for k, v in res.items() if not isinstance(v, (list, dict))}))
log.close()
