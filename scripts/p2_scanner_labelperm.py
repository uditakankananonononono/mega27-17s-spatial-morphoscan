"""Amendment-queue Tier-1 #13: label-permutation negative control for the
patient-held-out tumor-vs-other scanner (p2_scanner_partial.csv). Replicates
p2_model.py's scan classifier exactly (LogisticRegression(max_iter=2000, C=1.0)
on Standardized train-patient embeddings -> morphology labels), then on each
held-out section computes AUC vs (a) true labels and (b) 20 within-section
label permutations. Pre-declared pass: true AUCs reproduce the committed
scanner CSV (within 1e-6) and permuted-label AUCs concentrate at 0.5
(p95 < 0.6). If permuted AUCs stay high, the signal is an artifact. Writes
results/p2_scanner_labelperm.json."""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from src.morphoscan import models
import p2_common as C

N_PERM_LABEL = 20
rng = np.random.default_rng(20260928)
meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())
committed = {(r["patient"], r["section"]): float(r["auc"]) for r in _csv.DictReader(open("results/p2_scanner_partial.csv")) if r["auc"] != ""}
OUT = "results/p2_scanner_labelperm.json"
state = json.load(open(OUT)) if os.path.exists(OUT) else {"folds": []}
done = {f["patient"] for f in state["folds"]}

def assemble(pats):
    Xs, Ls = [], []
    for pat in pats:
        for sec in her2[her2["patient"] == pat]["section"]:
            if sec not in sections: continue
            s = sections[sec]; Xs.append(s["X"]); Ls.append(s["labels"])
    return np.vstack(Xs), np.concatenate(Ls)

for test_pat in patients:
    if test_pat in done: continue
    t0 = time.time()
    train_pats = [p for p in patients if p != test_pat]
    Xtr, Ltr = assemble(train_pats)
    import gc; gc.collect()
    sc = models.Standardizer().fit(Xtr)
    Ztr = sc.transform(Xtr)
    scan = LogisticRegression(max_iter=2000, C=1.0)
    scan.fit(Ztr, Ltr.astype(int))
    fold_rec = {"patient": test_pat, "sections": []}
    for sec in her2[her2["patient"] == test_pat]["section"]:
        if sec not in sections: continue
        s = sections[sec]
        y = s["labels"].astype(int)
        if not (y.any() and (~y.astype(bool)).any()): continue
        p_scan = scan.predict_proba(sc.transform(s["X"]))[:, 1]
        true_auc = float(roc_auc_score(y, p_scan))
        perm_aucs = [float(roc_auc_score(rng.permutation(y), p_scan)) for _ in range(N_PERM_LABEL)]
        fold_rec["sections"].append({
            "section": sec, "n_spots": int(len(y)), "tumor_frac": float(y.mean()),
            "true_auc": true_auc,
            "committed_auc": committed.get((test_pat, sec)),
            "perm_aucs": perm_aucs})
    state["folds"].append(fold_rec)
    json.dump(state, open(OUT, "w"))
    print(f"{test_pat}: {len(fold_rec['sections'])} sections ({time.time()-t0:.0f}s)", flush=True)

all_true, all_perm, devs = [], [], []
for f in state["folds"]:
    for sr in f["sections"]:
        all_true.append(sr["true_auc"]); all_perm.extend(sr["perm_aucs"])
        if sr["committed_auc"] is not None:
            devs.append(abs(sr["true_auc"] - sr["committed_auc"]))
verdict = {"n_sections": len(all_true), "n_label_perms_per_section": N_PERM_LABEL,
           "true_auc_median": float(np.median(all_true)),
           "max_abs_dev_vs_committed_scanner_csv": float(max(devs)) if devs else None,
           "perm_auc_median": float(np.median(all_perm)),
           "perm_auc_p95": float(np.percentile(all_perm, 95)),
           "criterion": "true reproduces committed (<=1e-6) and perm p95 < 0.6",
           "pass": bool((max(devs) if devs else 1) <= 1e-6 and np.percentile(all_perm, 95) < 0.6)}
print(json.dumps(verdict, indent=2))
state["verdict"] = verdict
json.dump(state, open(OUT, "w"))
