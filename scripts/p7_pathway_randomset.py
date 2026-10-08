"""P7-PW1: random-gene-set null for hallmark pathway spatial autocorrelation (Moran's I).
Prespec: results/p7_pathway_randomset_prespec.json (commit 1877e0a, before this script).
Sections: the 29 HER2+ sections in results/p2_niches.json (metadata.csv absent from rebuilt workspace).
Run: python3 scripts/p7_pathway_randomset.py [--smoke] -> results/p7_pathway_randomset.json"""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "scripts")
import numpy as np
import p2_common as C
SEED = 20261008; NNULL = 200; NB = 10000
rng = np.random.default_rng(SEED)
def morans_i(x, W):
    n = len(x); xc = x - x.mean(); den = (xc ** 2).sum()
    if den == 0: return np.nan
    return (n / W.sum()) * (W * np.outer(xc, xc)).sum() / den
secs = [r["section"] for r in json.load(open("results/p2_niches.json"))]
hallmark = C.load_hallmark()
all_secs = [s for s in secs if os.path.exists(os.path.join(C.FC, f"{s}.npz")) and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
pw_idx = C.pathway_idx(universe, C.get_symbol_map(universe), hallmark)
sections = C.load_sections(all_secs, upos)
print(f"sections {len(sections)}/{len(secs)} universe {len(universe)} pathways {len(pw_idx)}", flush=True)
if "--smoke" in sys.argv: sys.exit(0)
res = {}
for sec, s in sections.items():
    n = s["N"].shape[0]
    if n < 60: continue
    d = np.sqrt(((s["px"][:, None, :] - s["px"][None, :, :]) ** 2).sum(-1))
    W = np.exp(-(d ** 2) / (2 * (1.5 * s["pitch"]) ** 2)); np.fill_diagonal(W, 0.0)
    N = s["N"]; mu = np.asarray(N.mean(axis=0)).ravel()
    ex2 = np.asarray(N.multiply(N).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0.0)); sd[sd == 0] = 1.0
    Nc = N.tocsc()
    q = np.quantile(mu, [.2, .4, .6, .8]); qb = np.searchsorted(q, mu, side="right")
    pools = [np.where(qb == b)[0] for b in range(5)]
    for name, idx in pw_idx.items():
        sc = lambda ii: ((Nc[:, ii].toarray() - mu[ii]) / sd[ii]).mean(axis=1)
        t = sc(idx)
        if t.std() == 0: continue
        obs = morans_i(t, W)
        null = np.empty(NNULL)
        for b in range(NNULL):
            rs = np.array([rng.choice(pools[qb[g]]) for g in idx])
            tt = sc(rs); null[b] = morans_i(tt, W) if tt.std() > 0 else np.nan
        null = null[np.isfinite(null)]
        res.setdefault(name, []).append(dict(section=sec, obs=float(obs), null_median=float(np.median(null)),
                                             p=float((1 + (null >= obs).sum()) / (len(null) + 1))))
    print(sec, "done", flush=True)
out = {"prespec": "results/p7_pathway_randomset_prespec.json", "seed": SEED, "pathways": {}}
brng = np.random.default_rng(SEED + 1)
for name, L in res.items():
    eff = np.array([z["obs"] - z["null_median"] for z in L]); fr = float(np.mean([z["p"] < 0.05 for z in L]))
    bs = [np.median(brng.choice(eff, len(eff))) for _ in range(NB)]
    out["pathways"][name] = {"n_sections": len(L), "frac_p05": fr, "median_excess": float(np.median(eff)),
        "excess_ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
        "survives": bool(fr >= 0.6 and np.median(eff) > 0)}
P = out["pathways"]
def find(k): return next((n for n in P if k.lower() in n.lower()), None)
top = [find("Myc Targets V1"), find("Mitotic Spindle"), find("G2-M")]; ctl = find("Pancreas")
out["top3"] = top; out["negative_control"] = ctl
out["n_survive"] = int(sum(v["survives"] for v in P.values())); out["n_pathways"] = len(P)
out["D2_claim"] = "SURVIVES" if (all(P[t]["survives"] for t in top if t) and ctl and not P[ctl]["survives"]) else "NOT SUPPORTED"
json.dump(out, open("results/p7_pathway_randomset.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("top3", "negative_control", "n_survive", "n_pathways", "D2_claim")}))
for t in top + [ctl]:
    if t: print(t, P[t])
