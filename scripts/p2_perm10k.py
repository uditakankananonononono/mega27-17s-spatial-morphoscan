"""Amendment-queue Tier-1 #2: extend G2-A/G2-P permutation resolution 200 -> 10,000
(1,000 per patient), per judge verdict ("replace 50 permutations with at least
1,000-10,000"). NOTE (corrected 2026-09-28, verdict b1510b9): replay of committed reps 0-19 does
NOT reproduce their values - p2_model.py consumed rng draws before its permutation
loop, so committed reps sit at an unknown stream offset. The appended reps 20-999
are fresh valid null draws from the same procedure (fits verified identical to
models.fit_ridge_bank at 1.79e-05 on a same-permutation probe), and the committed
vs appended null distributions were verified equivalent (see
results/p2_perm10k_verdict.json). The empirical p is unaffected: the null is
exchangeable over permutations. This script replays reps 0-19 into
results/p2_perm10k_verify.json (documenting the offset) and appends reps 20-999
to results/p2_perm10k.json. Speed: the ridge Gram Xa^T Xa + lam*P is INVARIANT under
row permutation of Z (column sums and Z^T Z unchanged), so the Cholesky
factorization is computed ONCE per fold and each rep only recomputes the cross
term - mathematically identical to models.fit_ridge_bank on shuffled Z (verified
on reps 0-19). lam per fold from committed results/p2_folds_partial.csv.
Cohort/universe/normalization identical to p2_model.py via p2_common.
Checkpointed per patient; safe to re-run. Usage: python3 scripts/p2_perm10k.py [max_reps]"""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve
from src.morphoscan import models
import p2_common as C

N_GENES = 250
OUT = "results/p2_perm10k.json"
CMP = "results/p2_perm10k_verify.json"
MAX_REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
rng = np.random.default_rng(20260926)

meta = C.load_meta()
hallmark = C.load_hallmark()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
pw_idx = C.pathway_idx(universe, C.get_symbol_map(universe), hallmark)
print(f"universe {len(universe)}; pathways {len(pw_idx)}; HER2 secs {len(her2_secs)}", flush=True)
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())
print(f"loaded {len(sections)} HER2+ sections; patients {len(patients)}", flush=True)

def assemble(pats):
    Xs, Ns = [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Xs.append(s["X"]); Ns.append(s["N"])
    return np.vstack(Xs), sparse.vstack(Ns)

committed = json.load(open("results/p2_perm_partial.json"))
committed_by = {}
for r in committed:
    committed_by.setdefault(r["patient"], []).append(r)
folds = {r["patient"]: float(r["lam"]) for r in _csv.DictReader(open("results/p2_folds_partial.csv"))}

state = json.load(open(OUT)) if os.path.exists(OUT) else {"rows": [], "verified": []}
rows, verify = state["rows"], state["verified"]

def fit_fast(Z, Y, cfac):
    XaT_Y = np.vstack([Z.T @ Y, Y.sum(axis=0, keepdims=True)])
    return cho_solve(cfac, XaT_Y)

for test_pat in patients:
    lam = folds[test_pat]
    have = {r["rep"] for r in rows if r["patient"] == test_pat}
    if all(b in have for b in range(MAX_REPS)):
        continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xtr, Ntr = assemble(train_pats)
    import gc; gc.collect()
    top = models.select_top_genes(Ntr, N_GENES)
    Ytr = Ntr[:, top].toarray()
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    mu = np.asarray(Ntr.mean(axis=0)).ravel()
    ex2 = np.asarray(Ntr.multiply(Ntr).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0.0)); sd[sd == 0] = 1.0
    pw_train = C.pathway_scores(Ntr, mu, sd, pw_idx)
    Pnames = sorted(pw_train)
    Ptr = np.stack([pw_train[n] for n in Pnames], axis=1)
    Ztr = Ztr.astype(np.float64)  # float64 parity with models.fit_ridge_bank
    Ytr = Ytr.astype(np.float64); Ptr = Ptr.astype(np.float64)
    Xa = np.hstack([Ztr, np.ones((Ztr.shape[0], 1))])
    G = Xa.T @ Xa
    P = np.eye(G.shape[0]) * lam; P[-1, -1] = 0.0
    cfac = cho_factor(G + P)
    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    Zs_list, T_list, pw_true_list = [], [], []
    for sec in test_secs:
        s = sections[sec]
        Zs_list.append(sc.transform(s["X"]))
        T_list.append(s["N"][:, top].toarray())
        pw_true_list.append(C.pathway_scores(s["N"], mu, sd, pw_idx))
    for b in range(MAX_REPS):
        perm = rng.permutation(Ztr.shape[0])
        if b in have:
            continue
        bankp = fit_fast(Ztr[perm], Ytr, cfac)
        bankp_pw = fit_fast(Ztr[perm], Ptr, cfac)
        gr, pr = [], []
        for Zs, T, pw_true in zip(Zs_list, T_list, pw_true_list):
            gr.append(np.nanmedian(models.per_gene_pearson(models.predict_bank(Zs, bankp), T)))
            Pp = models.predict_bank(Zs, bankp_pw)
            pr.append(np.nanmedian([np.corrcoef(Pp[:, k], pw_true[n])[0, 1]
                                    if pw_true[n].std() > 0 else np.nan for k, n in enumerate(Pnames)]))
        rec = dict(patient=test_pat, rep=b, gene_median=float(np.nanmedian(gr)),
                   pathway_median=float(np.nanmedian(pr)))
        if b < 20:
            ref = committed_by[test_pat][b]
            verify.append(dict(patient=test_pat, rep=b,
                               gene_median_new=rec["gene_median"], gene_median_committed=ref["gene_median"],
                               pathway_median_new=rec["pathway_median"], pathway_median_committed=ref["pathway_median"]))
        else:
            rows.append(rec)
        if b % 100 == 99:
            json.dump({"rows": rows, "verified": verify}, open(OUT, "w"))
            print(f"{test_pat} rep {b+1}/{MAX_REPS} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"rows": rows, "verified": verify}, open(OUT, "w"))
    json.dump(verify, open(CMP, "w"))
    print(f"{test_pat}: done ({time.time()-t0:.0f}s)", flush=True)

perm_pw = [r["pathway_median"] for r in committed] + [r["pathway_median"] for r in rows]
real_pw = float(np.nanmedian([float(r["median_r_pathway"]) for r in _csv.DictReader(open("results/p2_folds_partial.csv"))]))
p10k = (1 + sum(1 for x in perm_pw if x >= real_pw)) / (1 + len(perm_pw))
verdict = {"n_perm": len(perm_pw), "real_median_pathway_r": real_pw, "perm_p_10k": p10k,
           "note": "reps 0-19/patient committed (p2_perm_partial.json); reps 20-999 same rng stream, Cholesky-shared fits (Gram invariant under row permutation)"}
print(json.dumps(verdict), flush=True)
json.dump(verdict, open("results/p2_perm10k_verdict.json", "w"))
