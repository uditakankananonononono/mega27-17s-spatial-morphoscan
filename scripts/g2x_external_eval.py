"""G2-X external eval: frozen model (results/g2x_frozen_model.npz, committed
e9f6ee4 BEFORE this ran) on 10x Visium Parent_Visium_Human_BreastCancer.
Per results/g2x_protocol_prespec.json v1.1: spot-centered tiles, side =
0.75 x median NN pitch in full-res pixels (IDENTICAL geometry rule to
scripts/p2_extract_ctranspath.py: half = 0.375 x median NN pitch), resized
224x224, CTransPath embeddings (same weights/recipe), frozen Standardizer +
ridge pathway bank, UNSMOOTHED predictions, external pathway scores from the
h5 with frozen universe/mu/sd. Gate (PREREG-2): median pathway Pearson >= 0.10.
Checkpointed per 500 spots. Env: G2X_DATA (default /home/sandbox/g2x).
Run: python scripts/g2x_external_eval.py -> results/g2x_external_eval.json"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, os.environ.get("CTRAN_DIR", "/home/sandbox/ctranspath_src"))
import numpy as np
from PIL import Image
import torch
from ctran import ConvStem
import timm
import tifffile
import h5py
from scipy import sparse
from src.morphoscan import models

G2X = os.environ.get("G2X_DATA", "/home/sandbox/g2x")
CKPT = os.environ.get("G2X_CKPT", "results/g2x_embed_ckpt.npz")
OUT = "results/g2x_external_eval.json"
Image.MAX_IMAGE_PIXELS = None
MEAN = np.array([0.485, 0.456, 0.406], np.float32)
STD = np.array([0.229, 0.224, 0.225], np.float32)
torch.set_num_threads(int(os.environ.get("P2_THREADS", "2")))

# --- frozen model (loaded first; fail fast if artifact missing)
fz = np.load("results/g2x_frozen_model.npz", allow_pickle=True)
bank, Pnames = fz["bank"], list(fz["Pnames"])
universe = [str(g) for g in fz["universe"]]
sc_mu, sc_sd, mu, sd = fz["sc_mean"], fz["sc_std"], fz["mu"], fz["sd"]
counts = fz["pw_idx_counts"]; flat = fz["pw_idx_flat"]
pw_idx, o = {}, 0
for n, c in zip(Pnames, counts):
    pw_idx[n] = list(flat[o:o + c]); o += c

# --- spot positions + pitch (geometry rule identical to internal extraction)
pos = {}
for line in open(os.path.join(G2X, "spatial/tissue_positions_list.csv")):
    p = line.rstrip("\n").split(",")
    if p[0] == "barcode" or len(p) < 6: continue
    if p[1] == "1":
        pos[p[0]] = (float(p[5]), float(p[4]))  # (x=col, y=row) full-res
bcs = sorted(pos)
xy = np.array([pos[b] for b in bcs])
from scipy.spatial import cKDTree
d, _ = cKDTree(xy).query(xy, k=2)
pitch = float(np.median(d[:, 1]))
half = 0.375 * pitch
print(f"in-tissue spots: {len(bcs)} | median NN pitch: {pitch:.1f}px | tile half: {half:.1f}px", flush=True)

# --- CTransPath (recipe identical to p2_extract_ctranspath.py)
model = timm.create_model("swin_tiny_patch4_window7_224", pretrained=False)
model.patch_embed = ConvStem(embed_dim=96, norm_layer=torch.nn.LayerNorm)
sdict = torch.load(os.environ.get("CTRAN_PTH", "/home/sandbox/models/ctranspath.pth"), map_location="cpu", weights_only=False, mmap=False)
sdict = sdict["model"] if "model" in sdict else sdict
sdict = {k: v.clone() for k, v in sdict.items()}
missing, unexpected = model.load_state_dict(sdict, strict=False)
del sdict
assert all(m.startswith("head.") for m in missing) and not unexpected, (missing, unexpected)
model.eval()

store = tifffile.imread(os.path.join(G2X, "image.tif"), aszarr=True)
import zarr
img = zarr.open(store, mode="r")
H, W = img.shape[0], img.shape[1]
print(f"image {W}x{H}", flush=True)

# --- embeddings, checkpointed per 500 spots
start = 0
if os.path.exists(CKPT):
    z = np.load(CKPT)
    Xs, done_bcs = list(z["Xs"]), list(z["bcs"])
    start = len(done_bcs)
else:
    Xs, done_bcs = [], []
t0 = time.time()
MAX_SPOTS = int(os.environ.get("MAX_SPOTS", "0")) or len(bcs)
for i in range(start, min(MAX_SPOTS, len(bcs))):
    x, y = pos[bcs[i]]
    x0, x1 = int(round(x - half)), int(round(x + half))
    y0, y1 = int(round(y - half)), int(round(y + half))
    tile = np.zeros((y1 - y0, x1 - x0, 3), np.uint8)
    sx0, sy0 = max(0, x0), max(0, y0)
    sx1, sy1 = min(W, x1), min(H, y1)
    if sx1 > sx0 and sy1 > sy0:
        patch = img[sy0:sy1, sx0:sx1]
        if patch.ndim == 2: patch = np.repeat(patch[..., None], 3, axis=2)
        tile[sy0 - y0:sy0 - y0 + patch.shape[0], sx0 - x0:sx0 - x0 + patch.shape[1]] = patch[..., :3]
    im = np.asarray(Image.fromarray(tile).resize((224, 224), Image.BILINEAR), np.float32) / 255.0
    t = torch.from_numpy(((im - MEAN) / STD).transpose(2, 0, 1)[None]).float()
    with torch.no_grad():
        ff = model.forward_features(t)  # identical to p2_extract_ctranspath.py
        f = ff.mean(dim=(1, 2)) if ff.dim() == 4 else (ff.mean(dim=1) if ff.dim() == 3 else ff)
        Xs.append(f[0].numpy().astype(np.float32))
    done_bcs.append(bcs[i])
    if (i + 1) % 500 == 0 or i + 1 == len(bcs):
        np.savez(CKPT, Xs=np.stack(Xs), bcs=np.array(done_bcs))
        print(f"embed {i+1}/{len(bcs)} ({time.time()-t0:.0f}s)", flush=True)
X = np.stack(Xs)
if MAX_SPOTS < len(bcs):
    print(f"SMOKE MODE: {MAX_SPOTS} spots embedded OK - skipping eval/final JSON", flush=True)
    sys.exit(0)

# --- external expression -> pathway scores with frozen universe/mu/sd
h5 = h5py.File(os.path.join(G2X, "Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5"), "r")
m = h5["matrix"]
feat_ids = [g.decode().split(".")[0] for g in m["features"]["id"][:]]
hbc = [b.decode() for b in m["barcodes"][:]]
N_all = sparse.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=(len(feat_ids), len(hbc))).T.tocsr()  # 10x stores CSC (features x barcodes); -> spots x genes
upos = {g: i for i, g in enumerate(universe)}
lut = np.array([upos.get(g, -1) for g in feat_ids])
mm = N_all.tocoo(); nl = lut[mm.col]; keep = nl >= 0
N = sparse.csr_matrix((mm.data[keep].astype(np.float32), (mm.row[keep], nl[keep])),
                      shape=(N_all.shape[0], len(universe)))
lib = np.asarray(N_all.sum(axis=1)).ravel(); lib[lib == 0] = 1.0
N.data *= np.repeat((1e4 / lib).astype(np.float32), np.diff(N.indptr))
np.log1p(N.data, out=N.data)
bpos = {b: i for i, b in enumerate(hbc)}
rows = [bpos[b] for b in done_bcs]
N = N[rows]

def pathway_scores(Nm, pw_idx):
    out = {}
    for n, idx in pw_idx.items():
        Z = Nm[:, idx].toarray()
        Z = (Z - mu[idx]) / sd[idx]
        out[n] = Z.mean(axis=1)
    return out

pw_true = pathway_scores(N, pw_idx)
Z = (X - sc_mu) / sc_sd
Pp = models.predict_bank(Z, bank)
per_pw = {}
for k, n in enumerate(Pnames):
    t, p = pw_true[n], Pp[:, k]
    per_pw[n] = float(np.corrcoef(p, t)[0, 1]) if t.std() > 0 and p.std() > 0 else float("nan")
vals = np.array(list(per_pw.values()), float)
med = float(np.nanmedian(vals))
verdict = "PASS" if med >= 0.10 else "FAIL"
res = {"gate": "G2-X (PREREG-2): frozen model, external Visium breast, median pathway Pearson >= 0.10",
       "frozen_model": "results/g2x_frozen_model.npz (committed before external use)",
       "external_cohort": "10x Parent_Visium_Human_BreastCancer v1.2.0 (recipe dfb54b8)",
       "n_spots": int(len(done_bcs)), "n_pathways": len(Pnames),
       "pitch_px_fullres": pitch, "tile_half_px": half,
       "median_pathway_r": med, "per_pathway_r": per_pw,
       "frac_pathways_r_pos": float(np.mean(vals > 0)),
       "verdict": verdict,
       "declared_gaps": "single external section; different tumor, not HER2+-selected; spot-level regional co-variation, not patient-level generalization"}
json.dump(res, open(OUT, "w"), indent=1)
print(f"G2-X median pathway r = {med:.4f} -> {verdict} (wall {time.time()-t0:.0f}s)", flush=True)
