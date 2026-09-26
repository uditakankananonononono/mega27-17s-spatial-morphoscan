"""HGB baseline (PREREG-2 judge round-2 adoption A7): HistGradientBoostingRegressor
on the same CTransPath embeddings, same LOPO folds, same train-only pathway
z-params, scoped to the 49 Hallmark pathway targets. Strong tabular baseline vs
the Phase-2 ridge on identical data/splits. Checkpointed per fold."""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import HistGradientBoostingRegressor
sys.path.insert(0, "scripts")
import p2_common as C

OUT = "results/p2_hgb_pathway_partial.csv"
log = open("/tmp/p2_hgb.log", "a")
def say(m): print(m); log.write(m + "\n"); log.flush()

meta = C.load_meta()
hallmark = C.load_hallmark()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
pw_idx = C.pathway_idx(universe, C.get_symbol_map(universe), hallmark)
say(f"universe {len(universe)}; pathways {len(pw_idx)}; HER2 secs {len(her2_secs)}")
sections = C.load_sections(her2_secs, upos)
say(f"loaded {len(sections)} HER2+ sections")
patients = sorted(her2["patient"].unique())

def assemble(pats):
    Xs, Ns = [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]
            Xs.append(s["X"]); Ns.append(s["N"])
    return np.vstack(Xs), sparse.vstack(Ns)

records = list(_csv.DictReader(open(OUT))) if os.path.exists(OUT) else []
done_pats = {r["patient"] for r in records}

def checkpoint():
    if records:
        with open(OUT, "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=list(records[0].keys())); w.writeheader(); w.writerows(records)

from src.morphoscan import models
for test_pat in patients:
    if test_pat in done_pats:
        say(f"skip {test_pat} (done)"); continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xtr, Ntr = assemble(train_pats)
    import gc; gc.collect()
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr).astype(np.float32)
    # train-only pathway z-params (identical to p2_model.py)
    mu = np.asarray(Ntr.mean(axis=0)).ravel()
    ex2 = np.asarray(Ntr.multiply(Ntr).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0.0)); sd[sd == 0] = 1.0
    pw_train = C.pathway_scores(Ntr, mu, sd, pw_idx)
    Pnames = sorted(pw_train)
    Ytr = np.stack([pw_train[n] for n in Pnames], axis=1)
    models_hgb = []
    for k, n in enumerate(Pnames):
        h = HistGradientBoostingRegressor(max_iter=100, learning_rate=0.08,  # halved 2026-09-27: fold time ~80min at 200 iters on this box; HGB config not gate-locked, documented in JUDGE_ROUNDS round-2 adoption A7
                                          early_stopping=False, random_state=20260927)
        h.fit(Ztr, Ytr[:, k])
        models_hgb.append(h)
    del Xtr, Ntr, Ztr, Ytr; gc.collect()
    for sec in her2[her2["patient"] == test_pat]["section"]:
        if sec not in sections: continue
        s = sections[sec]
        Zs = sc.transform(s["X"]).astype(np.float32)
        pw_true = C.pathway_scores(s["N"], mu, sd, pw_idx)
        for k, n in enumerate(Pnames):
            t = pw_true[n]; p = models_hgb[k].predict(Zs)
            r = np.corrcoef(p, t)[0, 1] if t.std() > 0 and p.std() > 0 else float("nan")
            records.append(dict(patient=test_pat, section=sec, pathway=n,
                                r_hgb=f"{r:.6f}" if r == r else ""))
    del models_hgb; gc.collect()
    checkpoint()
    say(f"fold {test_pat}: done ({time.time()-t0:.0f}s)")

say("HGB baseline complete; comparing vs ridge")
df = pd.DataFrame(records)
df["r_hgb"] = pd.to_numeric(df["r_hgb"], errors="coerce")
ridge = pd.read_csv("results/p2_per_pathway_partial.csv")
ridge["r"] = pd.to_numeric(ridge["r"], errors="coerce")
m = df.merge(ridge, on=["patient", "section", "pathway"], how="inner").dropna(subset=["r_hgb", "r"])
med_hgb = m.groupby(["patient", "section"])["r_hgb"].median().median()
med_ridge = m.groupby(["patient", "section"])["r"].median().median()
from scipy import stats
w = stats.wilcoxon(m["r_hgb"], m["r"])
verdict = {"n_pairs": int(len(m)), "median_hgb": float(med_hgb), "median_ridge": float(med_ridge),
           "wilcoxon_p": float(w.pvalue), "hgb_wins_frac": float((m["r_hgb"] > m["r"]).mean())}
json.dump(verdict, open("results/p2_hgb_verdict.json", "w"), indent=1)
say(json.dumps(verdict))
say("DONE")
