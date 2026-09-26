"""D2 (PREREG-2 judge round-2 adoption): morphology-defined molecular niches.
Per HER2+ section: freeze k=8, MiniBatchKMeans on CTransPath embeddings ->
niche labels; spatial coherence = Moran's I of niche indicators with gaussian
weights (length 1.5 x pitch), null = 200 coordinate-preserving label shuffles.
Molecular characterization: median true hallmark pathway scores per niche.
Also Moran's I of each true pathway score itself vs the same null (are the
molecular programs spatially organized, not just the morphology?)."""
import os, sys, json, time, csv as _csv
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
sys.path.insert(0, "scripts")
import p2_common as C
from sklearn.cluster import MiniBatchKMeans
from scipy import stats

K = 8            # frozen cluster count (round-2 log)
N_NULL = 200     # frozen null shuffles
rng = np.random.default_rng(20260927)

def morans_i(x, W):
    n = len(x); xc = x - x.mean()
    denom = (xc ** 2).sum()
    if denom == 0: return np.nan
    num = (W * np.outer(xc, xc)).sum()
    S0 = W.sum()
    return (n / S0) * num / denom

meta = C.load_meta()
hallmark = C.load_hallmark()
all_secs = [s for s in meta["section"] if os.path.exists(os.path.join(C.FC, f"{s}.npz"))
            and os.path.exists(os.path.join(C.EC, f"{s}.npz"))]
her2 = meta[meta["is_her2"]]
her2_secs = [s for s in her2["section"] if s in all_secs]
universe = C.build_universe(all_secs)
upos = {g: i for i, g in enumerate(universe)}
pw_idx = C.pathway_idx(universe, C.get_symbol_map(universe), hallmark)
sections = C.load_sections(her2_secs, upos)
print(f"loaded {len(sections)} HER2+ sections; pathways {len(pw_idx)}")

# whole-cohort pathway z-params (descriptive spatial analysis, not prediction;
# no train/test leak concept applies to unsupervised niche discovery)
rows = []
for sec, s in sections.items():
    n = s["X"].shape[0]
    if n < 60: continue
    d = np.sqrt(((s["px"][:, None, :] - s["px"][None, :, :]) ** 2).sum(-1))
    W = np.exp(-(d ** 2) / (2 * (1.5 * s["pitch"]) ** 2))
    np.fill_diagonal(W, 0.0)
    km = MiniBatchKMeans(n_clusters=K, random_state=20260927, batch_size=2048, n_init=3)
    lab = km.fit_predict(s["X"])
    Is = []
    for k in range(K):
        ind = (lab == k).astype(float)
        if ind.sum() < 5 or ind.sum() > n - 5: continue
        i_obs = morans_i(ind, W)
        null = np.empty(N_NULL)
        for b in range(N_NULL):
            null[b] = morans_i(ind[rng.permutation(n)], W)
        p_emp = (1 + (null >= i_obs).sum()) / (N_NULL + 1)
        Is.append((i_obs, p_emp))
    # pathway-score spatial autocorrelation (z-params from this section; descriptive)
    N = s["N"]; mu = np.asarray(N.mean(axis=0)).ravel()
    ex2 = np.asarray(N.multiply(N).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0.0)); sd[sd == 0] = 1.0
    pw = C.pathway_scores(N, mu, sd, pw_idx)
    pw_I = {}
    for name, t in pw.items():
        if t.std() == 0: continue
        i_obs = morans_i(t, W)
        null = np.empty(50)
        for b in range(50):
            null[b] = morans_i(t[rng.permutation(n)], W)
        pw_I[name] = dict(I=float(i_obs), p=float((1 + (null >= i_obs).sum()) / 51))
    # molecular profile per niche: median pathway score
    niche_prof = {}
    for k in range(K):
        m_ = lab == k
        if m_.sum() < 5: continue
        top = sorted(((name, float(np.median(t[m_]) - np.median(t))) for name, t in pw.items()),
                     key=lambda z: -abs(z[1]))[:5]
        niche_prof[int(k)] = dict(n_spots=int(m_.sum()), top_pathways=top)
    rows.append(dict(section=sec, n_spots=int(n),
                     niche_I_median=float(np.median([z[0] for z in Is])) if Is else None,
                     niche_I_max_p=float(max(z[1] for z in Is)) if Is else None,
                     n_niches=int(len(Is)),
                     pw_I= pw_I, niche_profiles=niche_prof))
    print(f"{sec}: niche I_med={rows[-1]['niche_I_median']:.3f} p_max={rows[-1]['niche_I_max_p']:.4f}")

json.dump(rows, open("results/p2_niches.json", "w"), indent=1)
medI = np.median([r["niche_I_median"] for r in rows if r["niche_I_median"] is not None])
frac_sig = np.mean([r["niche_I_max_p"] < 0.05 for r in rows if r["niche_I_max_p"] is not None])
pw_all = {}
for r in rows:
    for name, v in r["pw_I"].items():
        pw_all.setdefault(name, []).append(v)
pw_summary = {name: dict(I_median=float(np.median([v["I"] for v in vs])),
                         frac_p_lt_05=float(np.mean([v["p"] < 0.05 for v in vs])))
              for name, vs in pw_all.items()}
summary = dict(k=K, n_null=N_NULL, n_sections=len(rows),
               niche_I_median_across_sections=float(medI),
               frac_sections_all_niches_sig=float(frac_sig),
               pathway_spatial=pw_summary)
json.dump(summary, open("results/p2_niches_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "pathway_spatial"}, indent=1))
top = sorted(pw_summary.items(), key=lambda z: -z[1]["I_median"])[:10]
for name, v in top:
    print(f"  {name}: I={v['I_median']:.3f} sig_frac={v['frac_p_lt_05']:.2f}")
print("DONE")
