"""DeepSpot-M head-to-head on the G2-X external section.
Prespec: results/deepspotm_h2h_prespec.json (committed 50b8c85 BEFORE this run).
Same spots/tiles/geometry as G2-X (0.75x pitch, 224x224); model:
deepspotm==1.0.0, source='scgpt' (declared primary); query genes = symbols of
universe ENSGs in the G2-X pathway index, intersected with model.gene_names;
pathway scores z-scored with the FROZEN internal mu/sd (declared harmonization);
endpoint = median per-pathway Pearson vs observed on the same 4,325 spots.
Characterization, not a gate. Checkpointed per 250 spots.
Env: G2X_DATA (default /home/sandbox/g2x). Run: python scripts/deepspotm_h2h.py"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
from PIL import Image
import torch
import tifffile, h5py
from scipy import sparse

G2X = os.environ.get("G2X_DATA", "/home/sandbox/g2x")
CKPT = os.environ.get("H2H_CKPT", "results/deepspotm_h2h_ckpt.npz")
OUT = "results/deepspotm_h2h.json"
Image.MAX_IMAGE_PIXELS = None
torch.set_num_threads(int(os.environ.get("P2_THREADS", "4")))

# frozen internal pathway definitions (same artifact as G2-X)
fz = np.load("results/g2x_frozen_model.npz", allow_pickle=True)
Pnames = list(fz["Pnames"])
universe = [str(g) for g in fz["universe"]]
mu, sd = fz["mu"], fz["sd"]
counts = fz["pw_idx_counts"]; flat = fz["pw_idx_flat"]
pw_idx, o = {}, 0
for n, c in zip(Pnames, counts):
    pw_idx[n] = list(flat[o:o + c]); o += c
sym_map = json.load(open("data/ensg_symbol_map.json"))

# spots + geometry (identical rule to scripts/g2x_external_eval.py)
pos = {}
for line in open(os.path.join(G2X, "spatial/tissue_positions_list.csv")):
    p = line.rstrip("\n").split(",")
    if p[0] == "barcode" or len(p) < 6: continue
    if p[1] == "1":
        pos[p[0]] = (float(p[5]), float(p[4]))
bcs = sorted(pos)
xy = np.array([pos[b] for b in bcs])
from scipy.spatial import cKDTree
d, _ = cKDTree(xy).query(xy, k=2)
pitch = float(np.median(d[:, 1])); half = 0.375 * pitch
print(f"spots {len(bcs)} pitch {pitch:.1f} half {half:.1f}", flush=True)

from deepspotm import DeepSpotM
model, image_processor = DeepSpotM.from_pretrained("ratschlab/DeepSpotM", source="scgpt")
model = model.eval()
vocab = set(model.gene_names)
# query genes: symbols of pathway-member universe genes present in the model vocab
pw_syms = {}   # pathway -> list of (symbol, universe_index)
all_syms = set()
for n, idxs in pw_idx.items():
    rows = []
    for i in idxs:
        s = sym_map.get(universe[i])
        if s and s in vocab:
            rows.append((s, i)); all_syms.add(s)
    if len(rows) >= 5:
        pw_syms[n] = rows
genes = sorted(all_syms)
print(f"model vocab {len(vocab)} | query genes {len(genes)} | pathways >=5 queryable members {len(pw_syms)}/{len(Pnames)}", flush=True)

store = tifffile.imread(os.path.join(G2X, "image.tif"), aszarr=True)
import zarr
img = zarr.open(store, mode="r")
H, W = img.shape[0], img.shape[1]

start, PRED = 0, []
if os.path.exists(CKPT):
    z = np.load(CKPT); PRED = [z["P"]]; start = z["P"].shape[0]
t0 = time.time()
batch, batch_i = [], []
def flush():
    if not batch: return
    pv = torch.stack(batch)
    with torch.no_grad():
        out = model.predict_genes(pv, genes)
    PRED.append(out.float().numpy().astype(np.float32))
    batch.clear(); batch_i.clear()
for i in range(start, len(bcs)):
    x, y = pos[bcs[i]]
    x0, x1 = int(round(x - half)), int(round(x + half))
    y0, y1 = int(round(y - half)), int(round(y + half))
    tile = np.zeros((y1 - y0, x1 - x0, 3), np.uint8)
    sx0, sy0 = max(0, x0), max(0, y0); sx1, sy1 = min(W, x1), min(H, y1)
    if sx1 > sx0 and sy1 > sy0:
        patch = img[sy0:sy1, sx0:sx1]
        if patch.ndim == 2: patch = np.repeat(patch[..., None], 3, axis=2)
        tile[sy0 - y0:sy0 - y0 + patch.shape[0], sx0 - x0:sx0 - x0 + patch.shape[1]] = patch[..., :3]
    batch.append(image_processor(Image.fromarray(tile).resize((224, 224), Image.BILINEAR)))
    batch_i.append(i)
    if len(batch) >= 8: flush()
    if (i + 1) % 250 == 0 or i + 1 == len(bcs):
        flush()
        np.savez(CKPT, P=np.concatenate(PRED, axis=0))
        print(f"pred {i+1}/{len(bcs)} ({time.time()-t0:.0f}s)", flush=True)
flush()
P = np.concatenate(PRED, axis=0)  # spots x genes
gpos = {g: j for j, g in enumerate(genes)}

# observed expression (10x h5 is CSC: features x barcodes)
h5 = h5py.File(os.path.join(G2X, "Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5"), "r")
m = h5["matrix"]
feat_ids = [g.decode().split(".")[0] for g in m["features"]["id"][:]]
hbc = [b.decode() for b in m["barcodes"][:]]
N_all = sparse.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=(len(feat_ids), len(hbc))).T.tocsr()
upos = {g: i for i, g in enumerate(universe)}
lut = np.array([upos.get(g, -1) for g in feat_ids])
mm = N_all.tocoo(); nl = lut[mm.col]; keep = nl >= 0
N = sparse.csr_matrix((mm.data[keep].astype(np.float32), (mm.row[keep], nl[keep])),
                      shape=(N_all.shape[0], len(universe)))
lib = np.asarray(N_all.sum(axis=1)).ravel(); lib[lib == 0] = 1.0
N.data *= np.repeat((1e4 / lib).astype(np.float32), np.diff(N.indptr))
np.log1p(N.data, out=N.data)
bpos = {b: i for i, b in enumerate(hbc)}
N = N[[bpos[b] for b in bcs]]

per_pw = {}
for n, rows in pw_syms.items():
    t_idx = [i for _, i in rows]
    T = N[:, t_idx].toarray(); T = (T - mu[t_idx]) / sd[t_idx]
    t = T.mean(axis=1)
    p_idx = [gpos[s] for s, _ in rows]
    p = P[:, p_idx].mean(axis=1)
    per_pw[n] = float(np.corrcoef(p, t)[0, 1]) if t.std() > 0 and p.std() > 0 else float("nan")
vals = np.array(list(per_pw.values()), float)
med = float(np.nanmedian(vals))
res = {"prespec": "results/deepspotm_h2h_prespec.json",
       "model": "deepspotm==1.0.0 ratschlab/DeepSpotM source=scgpt",
       "n_spots": int(len(bcs)), "n_query_genes": len(genes), "n_pathways": len(per_pw),
       "median_pathway_r": med, "per_pathway_r": per_pw,
       "frac_pathways_r_pos": float(np.mean(vals > 0)),
       "comparison": {"our_frozen_model_g2x": -0.1547,
                      "deepspotm_published_zero_shot_band": 0.19},
       "note": "characterization benchmark, not a gate (prespec); PolyForm NC - research use, attribution in paper"}
json.dump(res, open(OUT, "w"), indent=1)
print(f"H2H DeepSpot-M median pathway r = {med:.4f} (ours -0.1547, their published band 0.19) (wall {time.time()-t0:.0f}s)", flush=True)
