"""Shared loaders for Phase-2 analyses (p2_hgb_baseline, p2_niches).
Extracted verbatim from scripts/p2_model.py so every Phase-2 analysis uses the
identical cohort, universe, symbol map, and no-leak normalization."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from scipy import sparse

FC = "features_cache"; EC = "embeddings_cache"
DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
MAP_PATH = "data/ensg_symbol_map.json"

def load_meta():
    meta = pd.read_csv(os.path.join(DATA, "metadata.csv"))
    meta["section"] = meta["count_matrix"].str.replace("_stdata.tsv.gz", "", regex=False)
    meta["is_her2"] = meta["type"].str.startswith("HER2")
    return meta

def load_hallmark():
    hallmark = {}
    for line in open("data/hallmark50.gmt"):
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 2:
            hallmark[parts[0]] = [g for g in parts[2:] if g]
    return hallmark

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

def pathway_idx(universe, sym_map, hallmark):
    upos = {g: i for i, g in enumerate(universe)}
    pw_idx = {}
    sym_of = {g: sym_map.get(g) for g in universe}
    for name, syms in hallmark.items():
        sset = set(syms)
        idx = [upos[g] for g in universe if sym_of[g] in sset]
        if len(idx) >= 5:
            pw_idx[name] = sorted(set(idx))
    return pw_idx

def pathway_scores(N, mu, sd, pw_idx):
    out = {}
    for name, idx in pw_idx.items():
        Z = N[:, idx].toarray()
        Z = (Z - mu[idx]) / sd[idx]
        out[name] = Z.mean(axis=1)
    return out
