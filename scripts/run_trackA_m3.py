"""M3 MLP nonlinear check - PREREG-declared model (M0/M1/M2/M3 ladder, no gate).
M3 was marked "not run" in the original pipeline; this script fulfills it.
PRE-DECLARED SPEC (committed before any M3 result is seen; descriptive, no gate):
  data/split/genes identical to committed Track A (run_trackA.py): HER2+
  patients, LOPO, top-250 genes selected on train folds only, library-size 1e4
  + log1p, same features_cache X as committed M1/M2 (same representation, so
  the comparison is apples-to-apples; Phase-1 geometry caveat applies equally).
  model: sklearn.neural_network.MLPRegressor, hidden (256, 128), relu, adam,
  lr_init 1e-3, batch_size 512, max_iter 150, early_stopping=True with
  validation_fraction 0.1 drawn from TRAIN spots only (n_iter_no_change 10),
  random_state 20260929; inputs standardized by train Standardizer.
  output: median + mean per-gene Pearson on held-out patients (same protocol as
  M1/M2), per-fold rows checkpointed to results/trackA_m3_partial.csv,
  final results/trackA_m3.json with spec + params count + per-fold timings.
Run: setsid nohup python3 scripts/run_trackA_m3.py >/tmp/m3.log 2>&1 &
"""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neural_network import MLPRegressor
from src.morphoscan import models

CACHE = "features_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
N_GENES = 250
log = open("/tmp/m3.log", "a")
def say(m): print(m); log.write(m + "\n"); log.flush()

meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(CACHE, f"{s}.npz"))]
her2_secs = [s for s in meta[meta["is_her2"]]["section"] if s in all_secs]

gene_sets = []
for sec in all_secs:
    z = np.load(os.path.join(CACHE, f"{sec}.npz"), allow_pickle=True)
    gene_sets.append(set(str(g) for g in z["genes"])); z.close()
universe = sorted(set.intersection(*gene_sets))
upos = {g: i for i, g in enumerate(universe)}

def load_sections(secs):
    out = {}
    for sec in secs:
        p = os.path.join(CACHE, f"{sec}.npz")
        if not os.path.exists(p): continue
        z = np.load(p, allow_pickle=True)
        X = z["X"]; ok = np.isfinite(X).all(axis=1)
        csr = z["counts_sparse"].item()[ok]
        genes = [str(g) for g in z["genes"]]
        z.close()
        lib = np.asarray(csr.sum(axis=1)).ravel(); lib[lib == 0] = 1.0
        lut = np.array([upos.get(g, -1) for g in genes])
        m = csr.tocoo(); nl = lut[m.col]; keep = nl >= 0
        N = sparse.csr_matrix((m.data[keep].astype(np.float32), (m.row[keep], nl[keep])),
                              shape=(m.shape[0], len(upos)), dtype=np.float32)
        N.data *= np.repeat((1e4 / lib).astype(np.float32), np.diff(N.indptr))
        np.log1p(N.data, out=N.data)
        out[sec] = dict(X=X[ok], N=N)
        del csr, m
    import gc; gc.collect()
    return out

sections = load_sections(her2_secs)
her2 = meta[meta["is_her2"]]
patients = sorted(her2["patient"].unique())
say(f"loaded {len(sections)} sections, {len(patients)} patients, universe {len(universe)}")

def assemble(pats):
    Xs, Ns = [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Xs.append(s["X"]); Ns.append(s["N"])
    return np.vstack(Xs), sparse.vstack(Ns)

PARTIAL = "results/trackA_m3_partial.csv"
done = {r["patient"] for r in _csv.DictReader(open(PARTIAL))} if os.path.exists(PARTIAL) else set()
fh = open(PARTIAL, "a", newline=""); wr = None
fold_rows, all_rs = [], []
for test_pat in patients:
    if test_pat in done: continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xtr, Ntr = assemble(train_pats)
    top = models.select_top_genes(Ntr, N_GENES)
    Ytr = Ntr[:, top].toarray()
    sc = models.Standardizer().fit(Xtr)
    mlp = MLPRegressor(hidden_layer_sizes=(256, 128), activation="relu", solver="adam",
                       learning_rate_init=1e-3, batch_size=512, max_iter=150,
                       early_stopping=True, validation_fraction=0.1, n_iter_no_change=10,
                       random_state=20260929)
    mlp.fit(sc.transform(Xtr), Ytr)
    rs = []
    for sec in her2[her2["patient"] == test_pat]["section"]:
        if sec not in sections: continue
        s = sections[sec]
        pred = mlp.predict(sc.transform(s["X"]))
        true = s["N"][:, top].toarray()
        rs.append(models.per_gene_pearson(pred, true))
    rs = np.concatenate(rs).astype(float)
    med = float(np.nanmedian(rs))
    row = {"patient": test_pat, "n_genes": len(rs), "median_r_m3": med,
           "mean_r_m3": float(np.nanmean(rs)), "n_iter": int(mlp.n_iter_),
           "seconds": round(time.time() - t0, 1)}
    if wr is None:
        wr = _csv.DictWriter(fh, fieldnames=list(row)); 
        if os.path.getsize(PARTIAL) == 0: wr.writeheader()
    wr.writerow(row); fh.flush()
    fold_rows.append(row); all_rs.append(rs)
    say(f"{test_pat}: median {med:.4f} iters {mlp.n_iter_} {row['seconds']}s")
    del Xtr, Ntr, Ytr, mlp
    import gc; gc.collect()

prev = {r["patient"]: float(r["median_r_m3"]) for r in _csv.DictReader(open(PARTIAL))}
meds = sorted(prev.values())
R = np.concatenate(all_rs) if all_rs else np.array([np.nan])
n_params = sum(w.size for w in mlp.coefs_) + sum(b.size for b in mlp.intercepts_)
res = {"item": "M3 MLP nonlinear check (PREREG model ladder, descriptive, no gate)",
       "spec": {"hidden": [256, 128], "activation": "relu", "solver": "adam",
                "lr_init": 1e-3, "batch_size": 512, "max_iter": 150,
                "early_stopping": True, "val_fraction_train_only": 0.1,
                "random_state": 20260929, "n_params": int(n_params)},
       "median_r_m3_all_folds": float(np.median(meds)), "per_fold": prev,
       "protocol": "identical folds/features/genes protocol as committed M1/M2 (run_trackA.py)"}
json.dump(res, open("results/trackA_m3.json", "w"), indent=1)
say(json.dumps(res))
