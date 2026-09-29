"""enet SHARD leg: computes ONE fold's elastic-net leg, byte-identical in
statistics to scripts/p2_multitask.py PASS 2 (same build_fold, same
MultiTaskElasticNet config incl. random_state=20260929). Parallelization of
execution only - user directive 2026-09-29 1:54 PM via parent. Writes a
per-fold fragment results/enet_shards/<patient>.json; merge step
(p2_enet_merge.py) folds fragments into results/p2_multitask.json.
Run: python scripts/p2_enet_shard.py <patient>"""
import os, sys, json, time, csv as _csv, gc
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from scipy import sparse
from sklearn.linear_model import MultiTaskElasticNet
from src.morphoscan import models
import p2_common as C

N_GENES = 250
OUT = "results/p2_multitask.json"
test_pat = sys.argv[1]

state = json.load(open(OUT)) if os.path.exists(OUT) else {"folds": []}
for f in state.get("folds", []):
    if f["patient"] == test_pat and "elasticnet_median_r" in f:
        print(f"{test_pat} already has elasticnet_median_r={f['elasticnet_median_r']:.6f} in main JSON - shard no-op")
        os.makedirs("results/enet_shards", exist_ok=True)
        json.dump({"patient": test_pat, "skipped": "already in main JSON"},
                  open(f"results/enet_shards/{test_pat}.json", "w"))
        sys.exit(0)

meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())
assert test_pat in patients, f"unknown patient {test_pat}"

def per_gene_median(P, T):
    rs = []
    for j in range(T.shape[1]):
        if T[:, j].std() > 0 and P[:, j].std() > 0:
            rs.append(np.corrcoef(P[:, j], T[:, j])[0, 1])
    return float(np.nanmedian(rs))

# build_fold copied VERBATIM from scripts/p2_multitask.py
train_pats = [p for p in patients if p != test_pat]
Xs, Ns = [], []
for pat in train_pats:
    for sec in her2[her2["patient"] == pat]["section"]:
        if sec not in sections: continue
        s = sections[sec]; Xs.append(s["X"]); Ns.append(s["N"])
Xtr, Ntr = np.vstack(Xs), sparse.vstack(Ns)
gc.collect()
top = models.select_top_genes(Ntr, N_GENES)
Ytr = Ntr[:, top].toarray()
sc = models.Standardizer().fit(Xtr)
Ztr = sc.transform(Xtr)
test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
Zte = np.vstack([sc.transform(sections[s]["X"]) for s in test_secs])
Tte = np.vstack([sections[s]["N"][:, top].toarray() for s in test_secs])

t0 = time.time()
en = MultiTaskElasticNet(alpha=1e-3, l1_ratio=0.5, max_iter=500, tol=2e-3, selection="random", random_state=20260929)
en.fit(Ztr, Ytr)
rec = {"patient": test_pat,
       "elasticnet_median_r": per_gene_median(en.predict(Zte), Tte),
       "elasticnet_n_iter": int(en.n_iter_),
       "wall_seconds": round(time.time() - t0, 1),
       "enet_config": "alpha=1e-3 l1_ratio=0.5 max_iter=500 tol=2e-3 selection=random random_state=20260929 (identical to p2_multitask.py PASS 2)"}
os.makedirs("results/enet_shards", exist_ok=True)
json.dump(rec, open(f"results/enet_shards/{test_pat}.json", "w"), indent=1)
print(f"SHARD DONE {test_pat} enet {rec['elasticnet_median_r']:.4f} n_iter {rec['elasticnet_n_iter']} ({rec['wall_seconds']}s)", flush=True)
