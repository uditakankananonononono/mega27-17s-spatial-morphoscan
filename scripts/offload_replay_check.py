"""Offload replay check (run on the remote box BEFORE any heavy leg):
re-runs pass-1 fold BC23287 (PLS/RRR/ridge-refold) and asserts equality with
the committed results/p2_multitask.json values to 1e-6. Guards environment,
data-integrity, and harness identity on the offloaded machine. Exit 0 = match."""
import json, sys
sys.path.insert(0, "scripts")
sys.path.insert(0, ".")
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from src.morphoscan import models
import p2_common as C
import os, csv as _csv

TEST_PAT = "BC23287"
N_GENES, RANK = 250, 25
meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections([s for s in her2["section"] if s in all_secs], upos)
patients = sorted(her2["patient"].unique())
fr = {r["patient"]: r for r in _csv.DictReader(open("results/p2_folds_partial.csv"))}
lam = float(fr[TEST_PAT]["lam"])
train_secs = [sec for pat in patients if pat != TEST_PAT
              for sec in her2[her2["patient"] == pat]["section"] if sec in sections]
Xtr = np.vstack([sections[s]["X"] for s in train_secs])
Ntr_list = [sections[s]["N"] for s in train_secs]
from scipy import sparse
Ntr = sparse.vstack(Ntr_list)
top = models.select_top_genes(Ntr, N_GENES)
Ytr = Ntr[:, top].toarray()
sc = models.Standardizer().fit(Xtr)
Ztr = sc.transform(Xtr)
test_secs = [s for s in her2[her2["patient"] == TEST_PAT]["section"] if s in sections]
Zte = np.vstack([sc.transform(sections[s]["X"]) for s in test_secs])
Tte = np.vstack([sections[s]["N"][:, top].toarray() for s in test_secs])

def med(P, T):
    rs = [np.corrcoef(P[:, j], T[:, j])[0, 1] for j in range(T.shape[1])
          if T[:, j].std() > 0 and P[:, j].std() > 0]
    return float(np.nanmedian(rs))

pls = PLSRegression(n_components=RANK); pls.fit(Ztr, Ytr)
bank = models.fit_ridge_bank(Ztr, Ytr, lam)
U, S, Vt = np.linalg.svd(bank[:-1], full_matrices=False)
k = min(RANK, len(S))
Brr = U[:, :k] @ np.diag(S[:k]) @ Vt[:k]
Zte_a = np.hstack([Zte, np.ones((Zte.shape[0], 1))])
got = {"pls_median_r": med(pls.predict(Zte), Tte),
       "rrr_median_r": med(Zte_a @ np.vstack([Brr, bank[-1:]]), Tte),
       "ridge_m1_refold": med(models.predict_bank(Zte, bank), Tte)}
committed = json.load(open("results/p2_multitask.json"))
ref = [f for f in committed["folds"] if f["patient"] == TEST_PAT][0]
ok = True
for m, v in got.items():
    diff = abs(v - ref[m])
    print(f"{m}: remote {v:.8f} committed {ref[m]:.8f} diff {diff:.2e}")
    ok = ok and diff < 1e-6
print("REPLAY", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
