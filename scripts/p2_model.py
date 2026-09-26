"""Phase 2: CTransPath-embedding ridge, LOPO on HER2+ patients, per PREREG-2.

Identical cohort, split, normalization, and no-leak discipline as Track A
(scripts/run_trackA.py): gene universe = Ensembl intersection across ALL
sections; counts from features_cache (geometry-independent, unaffected by the
Phase-1 bug); embeddings + corrected pixel coords from embeddings_cache.
Targets: top-250 genes (train-selected, same selector) + Hallmark pathway
scores (mean of z-scored member-gene log-normalized expression, z-parameters
from TRAIN folds only). Checkpointed per fold. Gates evaluated at the end.
"""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from src.morphoscan import models, metrics

FC = "features_cache"
EC = "embeddings_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
N_GENES = 250
N_PERM = 20
rng = np.random.default_rng(20260926)
log = open("/tmp/p2_model.log", "a")
def say(m): print(m); log.write(m + "\n"); log.flush()

meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
meta["is_her2"] = meta["type"].str.startswith("HER2")

# ---- Hallmark-50 (Enrichr mirror of MSigDB_Hallmark_2020; 49 sets present) ----
hallmark = {}
for line in open("data/hallmark50.gmt"):
    parts = line.rstrip("\n").split("\t")
    if len(parts) > 2:
        hallmark[parts[0]] = [g for g in parts[2:] if g]
say(f"hallmark sets: {len(hallmark)}")

IMMUNE_PROGRAMS = ["Interferon Alpha Response", "Interferon Gamma Response",
                   "Inflammatory Response", "Complement", "Allograft Rejection"]

# ---- Ensembl -> symbol map (cached; MyGene batch POST) ----
MAP_PATH = "data/ensg_symbol_map.json"

def get_symbol_map(universe):
    if os.path.exists(MAP_PATH):
        return json.load(open(MAP_PATH))
    import urllib.request
    out = {}
    B = 1000
    for i in range(0, len(universe), B):
        chunk = universe[i:i+B]
        req = urllib.request.Request(
            "https://mygene.info/v3/query",
            data=("q=" + ",".join(chunk) + "&scopes=ensembl.gene&fields=symbol&species=human").encode(),
            headers={"User-Agent": "mega27-research/1.0", "Content-Type": "application/x-www-form-urlencoded"})
        with urllib.request.urlopen(req, timeout=60) as r:
            for hit in json.loads(r.read()):
                if "symbol" in hit and not hit.get("notfound"):
                    out[hit["query"]] = hit["symbol"]
        time.sleep(0.3)
    json.dump(out, open(MAP_PATH, "w"))
    return out

def build_universe(secs):
    gene_sets = []
    for sec in secs:
        z = np.load(os.path.join(FC, f"{sec}.npz"), allow_pickle=True)
        gene_sets.append(set(str(g) for g in z["genes"]))
        z.close()
    return sorted(set.intersection(*gene_sets))

def load_sections(secs, upos):
    """Join features_cache counts/labels with embeddings_cache X/px on spot_ids."""
    out = {}
    for sec in secs:
        fp = os.path.join(FC, f"{sec}.npz"); ep = os.path.join(EC, f"{sec}.npz")
        if not (os.path.exists(fp) and os.path.exists(ep)): continue
        zf = np.load(fp, allow_pickle=True); ze = np.load(ep, allow_pickle=True)
        f_ids = [str(s) for s in zf["spot_ids"]]
        e_ids = [str(s) for s in ze["spot_ids"]]
        e_pos = {s: i for i, s in enumerate(e_ids)}
        order = [e_pos[s] for s in f_ids if s in e_pos]
        keep_f = [i for i, s in enumerate(f_ids) if s in e_pos]
        if len(order) < 20:
            zf.close(); ze.close(); continue
        X = ze["X"][order].astype(np.float32)
        px = np.stack([ze["sx"][order], ze["sy"][order]], axis=1)
        csr = zf["counts_sparse"].item()[keep_f]
        genes = [str(g) for g in zf["genes"]]
        labels = zf["labels"][keep_f]
        pitch = float(zf["pitch"][0])
        zf.close(); ze.close()
        ok = np.isfinite(X).all(axis=1)
        X, px, csr, labels = X[ok], px[ok], csr[ok], labels[ok]
        lib = np.asarray(csr.sum(axis=1)).ravel(); lib[lib == 0] = 1.0
        lut = np.array([upos.get(g, -1) for g in genes])
        m = csr.tocoo(); nl = lut[m.col]; keep = nl >= 0
        N = sparse.csr_matrix((m.data[keep].astype(np.float32), (m.row[keep], nl[keep])),
                              shape=(m.shape[0], len(upos)), dtype=np.float32)
        N.data *= np.repeat((1e4 / lib).astype(np.float32), np.diff(N.indptr))
        np.log1p(N.data, out=N.data)
        out[sec] = dict(X=X, N=N, px=px, labels=labels, pitch=pitch)
        del csr, m, N
    import gc; gc.collect()
    return out

all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(FC, f"{s}.npz"))
            and os.path.exists(os.path.join(EC, f"{s}.npz"))]
her2_secs = [s for s in meta[meta["is_her2"]]["section"] if s in all_secs]
universe = build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
say(f"gene universe (intersection): {len(universe)}; sections with both caches: {len(all_secs)}, HER2+: {len(her2_secs)}")

sym_map = get_symbol_map(universe)
# pathway -> universe indices (via symbols)
pw_idx = {}
sym_of = {g: sym_map.get(g) for g in universe}
for name, syms in hallmark.items():
    sset = set(syms)
    idx = [upos[g] for g in universe if sym_of[g] in sset]
    if len(idx) >= 5:
        pw_idx[name] = sorted(set(idx))
say(f"pathways with >=5 mapped member genes: {len(pw_idx)}")

sections = load_sections(her2_secs, upos)
say(f"loaded {len(sections)} HER2+ sections")

her2 = meta[meta["is_her2"]]
patients = sorted(her2["patient"].unique())
say(f"HER2+ patients: {len(patients)}")

def assemble(pats):
    Xs, Ns, Ls, Ps = [], [], [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Xs.append(s["X"]); Ns.append(s["N"]); Ls.append(s["labels"]); Ps.append(s["px"])
    return np.vstack(Xs), sparse.vstack(Ns), np.concatenate(Ls), np.vstack(Ps)

def pathway_scores(N, mu, sd):
    """N: spots x universe log-normed CSR. mu/sd: per-gene train params (dense).
    Slices member-gene columns only - never densifies the full matrix."""
    out = {}
    for name, idx in pw_idx.items():
        Z = N[:, idx].toarray()
        Z = (Z - mu[idx]) / sd[idx]
        out[name] = Z.mean(axis=1)
    return out

# checkpointing
PF = "results/p2_folds_partial.csv"
fold_rows = list(_csv.DictReader(open(PF))) if os.path.exists(PF) else []
done_pats = {r["patient"] for r in fold_rows}
PERM_PATH = "results/p2_perm_partial.json"
perm_rows = json.load(open(PERM_PATH)) if os.path.exists(PERM_PATH) else []
PG = "results/p2_per_gene_partial.csv"
per_gene_records = list(_csv.DictReader(open(PG))) if os.path.exists(PG) else []
PP = "results/p2_per_pathway_partial.csv"
per_pw_records = list(_csv.DictReader(open(PP))) if os.path.exists(PP) else []
SCAN = "results/p2_scanner_partial.csv"
scan_records = list(_csv.DictReader(open(SCAN))) if os.path.exists(SCAN) else []

def checkpoint():
    if fold_rows:
        with open(PF, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=list(fold_rows[0].keys())); w.writeheader(); w.writerows(fold_rows)
    json.dump(perm_rows, open(PERM_PATH, "w"))
    if per_gene_records:
        with open(PG, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=list(per_gene_records[0].keys())); w.writeheader(); w.writerows(per_gene_records)
    if per_pw_records:
        with open(PP, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=list(per_pw_records[0].keys())); w.writeheader(); w.writerows(per_pw_records)
    if scan_records:
        with open(SCAN, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=list(scan_records[0].keys())); w.writeheader(); w.writerows(scan_records)

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

for test_pat in patients:
    if test_pat in done_pats:
        say(f"skip {test_pat} (done)"); continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xtr, Ntr, Ltr, _ = assemble(train_pats)
    top = models.select_top_genes(Ntr, N_GENES)
    Ytr = Ntr[:, top].toarray()
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    lam = models.choose_lambda_inner(Ztr, Ytr, np.repeat(range(len(train_pats)),
        [sum(sections[s]["X"].shape[0] for s in her2[her2["patient"] == p]["section"] if s in sections) for p in train_pats]))
    bank = models.fit_ridge_bank(Ztr, Ytr, lam)

    # pathway z-params from TRAIN only
    # sparse-safe mean/std over train spots (log1p values; zeros included)
    n_tr = Ntr.shape[0]
    mu = np.asarray(Ntr.mean(axis=0)).ravel()
    ex2 = np.asarray(Ntr.multiply(Ntr).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0.0)); sd[sd == 0] = 1.0
    pw_train = pathway_scores(Ntr, mu, sd)
    Pnames = sorted(pw_train)
    Ptr = np.stack([pw_train[n] for n in Pnames], axis=1)
    bank_pw = models.fit_ridge_bank(Ztr, Ptr, lam)

    # spatial alpha via inner LOPO on train patients (gene task, same as Track A)
    inner_pred, inner_true, inner_px = [], [], []
    for ip in train_pats:
        others = [p for p in train_pats if p != ip]
        if not others: continue
        Xi, Ni, _, _ = assemble(others)
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

    # scanner: patient-held-out tumor-vs-other logistic on embeddings
    scan = LogisticRegression(max_iter=2000, C=1.0)
    scan.fit(Ztr, Ltr)

    test_secs = [s for s in her2[her2["patient"] == test_pat]["section"] if s in sections]
    r1_all, r2_all = [], []
    pw_r = {n: [] for n in Pnames}
    immune_auc = []
    for sec in test_secs:
        s = sections[sec]
        Zs = sc.transform(s["X"])
        P1 = models.predict_bank(Zs, bank)
        T = s["N"][:, top].toarray()
        d = np.sqrt(((s["px"][:, None, :] - s["px"][None, :, :]) ** 2).sum(-1))
        W = metrics.gaussian_kernel_weights(d, 1.5 * s["pitch"])
        P2 = metrics.spatial_smooth(P1, W, alpha)
        r1 = models.per_gene_pearson(P1, T); r2 = models.per_gene_pearson(P2, T)
        r1_all.append(r1); r2_all.append(r2)
        for j in range(len(top)):
            per_gene_records.append(dict(patient=test_pat, section=sec, gene=universe[top[j]],
                                         r_m1=float(r1[j]), r_m2=float(r2[j])))
        pw_true = pathway_scores(s["N"], mu, sd)
        Pp = models.predict_bank(Zs, bank_pw)
        for k, n in enumerate(Pnames):
            t = pw_true[n]; p = Pp[:, k]
            r = np.corrcoef(p, t)[0, 1] if t.std() > 0 and p.std() > 0 else np.nan
            pw_r[n].append(r)
            per_pw_records.append(dict(patient=test_pat, section=sec, pathway=n, r=float(r) if r == r else ""))
            if n in IMMUNE_PROGRAMS and t.std() > 0:
                hi = t > np.median(t)
                if hi.any() and (~hi).any():
                    immune_auc.append((n, sec, float(roc_auc_score(hi.astype(int), p))))
        p_scan = scan.predict_proba(Zs)[:, 1]
        scan_records.append(dict(patient=test_pat, section=sec, n_spots=len(s["labels"]),
                                 tumor_frac=float(s["labels"].mean()),
                                 auc=float(roc_auc_score(s["labels"], p_scan)) if s["labels"].any() and (~s["labels"].astype(bool)).any() else ""))

    # G2-A / G2-P permutation: shuffled-embedding controls (gene + pathway)
    for b in range(N_PERM):
        perm = rng.permutation(Ztr.shape[0])
        bankp = models.fit_ridge_bank(Ztr[perm], Ytr, lam)
        bankp_pw = models.fit_ridge_bank(Ztr[perm], Ptr, lam)
        gr, pr = [], []
        for sec in test_secs:
            s = sections[sec]
            Zs = sc.transform(s["X"])
            T = s["N"][:, top].toarray()
            gr.append(np.nanmedian(models.per_gene_pearson(models.predict_bank(Zs, bankp), T)))
            Pp = models.predict_bank(Zs, bankp_pw)
            pw_true = pathway_scores(s["N"], mu, sd)
            pr.append(np.nanmedian([np.corrcoef(Pp[:, k], pw_true[n])[0, 1]
                                    if pw_true[n].std() > 0 else np.nan for k, n in enumerate(Pnames)]))
        perm_rows.append(dict(patient=test_pat, rep=b, gene_median=float(np.nanmedian(gr)),
                              pathway_median=float(np.nanmedian(pr))))

    m1 = float(np.nanmedian(np.concatenate(r1_all))); m2 = float(np.nanmedian(np.concatenate(r2_all)))
    med_pw = {n: float(np.nanmedian(pw_r[n])) for n in Pnames}
    fold_rows.append(dict(patient=test_pat, n_sections=len(test_secs), lam=lam, alpha=alpha,
                          median_r_m1=m1, median_r_m2=m2,
                          median_r_pathway=float(np.nanmedian(list(med_pw.values()))),
                          immune_auc_median=float(np.median([a for _, _, a in immune_auc])) if immune_auc else ""))
    say(f"fold {test_pat}: secs={len(test_secs)} lam={lam} alpha={alpha} r1={m1:.3f} r2={m2:.3f} "
        f"pw={fold_rows[-1]['median_r_pathway']:.3f} ({time.time()-t0:.0f}s)")
    checkpoint()

say("ALL FOLDS DONE")
checkpoint()
