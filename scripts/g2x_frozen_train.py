"""G2-X frozen model: all-HER2+ CTransPath-ridge pathway model, per
results/g2x_protocol_prespec.json v1.1 (c1b0f44 + correction commit).
Frozen artifact committed BEFORE any external embedding extraction.
Run: python scripts/g2x_frozen_train.py -> results/g2x_frozen_model.npz"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
from src.morphoscan import models
import p2_common as C

meta = C.load_meta()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
sections = C.load_sections(her2_secs, upos)
patients = sorted(her2["patient"].unique())
print(f"HER2+ sections: {len(sections)} | patients: {len(patients)} | universe: {len(universe)}")

Xs, Ns, grp = [], [], []
for i, p in enumerate(patients):
    for sec in her2[her2["patient"] == p]["section"]:
        if sec not in sections: continue
        s = sections[sec]
        Xs.append(s["X"]); Ns.append(s["N"])
        grp += [i] * s["X"].shape[0]
X = np.vstack(Xs); N = sparse_vstack = __import__("scipy.sparse", fromlist=["vstack"]).vstack(Ns)
grp = np.array(grp)
print(f"spots: {X.shape[0]}")

hallmark = C.load_hallmark()
sym_map = C.get_symbol_map(universe)
pw_idx = C.pathway_idx(universe, sym_map, hallmark)
print(f"pathways with >=5 mapped members: {len(pw_idx)}")

sc = models.Standardizer().fit(X)
Z = sc.transform(X)
mu = np.asarray(N.mean(axis=0)).ravel()
ex2 = np.asarray(N.multiply(N).mean(axis=0)).ravel()
sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0.0)); sd[sd == 0] = 1.0
pw = C.pathway_scores(N, mu, sd, pw_idx)
Pnames = sorted(pw)
P = np.stack([pw[n] for n in Pnames], axis=1)

lam = models.choose_lambda_inner(Z, P, grp)
print(f"chosen lambda (pathway task, full-cohort inner CV): {lam}")
bank = models.fit_ridge_bank(Z, P, lam)

# sanity: in-sample pathway r (record only, NOT the gate - gate is external)
Pred = models.predict_bank(Z, bank)
r_ins = {}
for k, n in enumerate(Pnames):
    t, p = P[:, k], Pred[:, k]
    r_ins[n] = float(np.corrcoef(p, t)[0, 1]) if t.std() > 0 and p.std() > 0 else float("nan")
med = float(np.nanmedian(list(r_ins.values())))
print(f"in-sample median pathway r (sanity, not gate): {med:.4f}")

np.savez("results/g2x_frozen_model.npz",
         bank=bank, lam=np.array([lam]), Pnames=np.array(Pnames),
         pw_idx_flat=np.array([i for n in Pnames for i in pw_idx[n]]),
         pw_idx_counts=np.array([len(pw_idx[n]) for n in Pnames]),
         universe=np.array(universe),
         sc_mean=sc.mu, sc_std=sc.sd, mu=mu, sd=sd)
json.dump({"protocol": "results/g2x_protocol_prespec.json v1.1",
           "n_sections": len(sections), "n_patients": len(patients), "n_spots": int(X.shape[0]),
           "n_pathways": len(Pnames), "lam": float(lam),
           "in_sample_median_pathway_r_sanity": med},
          open("results/g2x_frozen_model_meta.json", "w"), indent=1)
print("froze results/g2x_frozen_model.npz + meta")
