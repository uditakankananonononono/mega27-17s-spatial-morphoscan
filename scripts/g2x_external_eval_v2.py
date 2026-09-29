"""G2-X v2: stain-normalized external re-eval (prespec results/g2x_v2_prespec.json,
def3b27). IDENTICAL to scripts/g2x_external_eval.py except ONE input-processing
step: each external tile is Reinhard-normalized in CIELAB space to reference
statistics computed from internal BC23901_C2 tiles (p2 geometry, resize 224,
stats over all tile pixels). Reference stats committed to
results/g2x_v2_refstain.json on first run. Same frozen model, same external
section/spots/tiles, same UNSMOOTHED pathway eval, same gate (median r >= 0.10).
v1 FAIL (34afa2f) stays on record; v2 carries full lineage disclosure.
Env: G2X_DATA, HER2ST_RAW, G2X_V2_CKPT, MAX_SPOTS smoke mode.
Run: python scripts/g2x_external_eval_v2.py -> results/g2x_external_eval_v2.json"""
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
from skimage.color import rgb2lab, lab2rgb
from src.morphoscan import models

G2X = os.environ.get("G2X_DATA", "/home/sandbox/g2x")
HER2ST = os.environ.get("HER2ST_RAW", "/home/sandbox/her2st_raw")
CKPT = os.environ.get("G2X_V2_CKPT", "results/g2x_v2_embed_ckpt.npz")
REF = "results/g2x_v2_refstain.json"
OUT = "results/g2x_external_eval_v2.json"
Image.MAX_IMAGE_PIXELS = None
MEAN = np.array([0.485, 0.456, 0.406], np.float32)
STD = np.array([0.229, 0.224, 0.225], np.float32)
torch.set_num_threads(int(os.environ.get("P2_THREADS", "4")))

# --- reference LAB stats from internal BC23901_C2 tiles (p2 geometry)
if os.path.exists(REF):
    rs = json.load(open(REF))
else:
    MAXDIM = 4500
    import pandas as pd
    spots = pd.read_csv(os.path.join(HER2ST, "spots_BT23901_C2.csv.gz"), index_col=0)
    img = Image.open(os.path.join(HER2ST, "HE_BT23901_C2.jpg"))
    orig_w, orig_h = img.size
    img.draft("RGB", (MAXDIM, MAXDIM)); img = img.convert("RGB")
    if max(img.size) > MAXDIM:
        r = MAXDIM / max(img.size)
        img = img.resize((int(img.size[0] * r), int(img.size[1] * r)), Image.LANCZOS)
    arr = np.asarray(img)
    sx = spots["X"].values * (arr.shape[1] / orig_w)
    sy = spots["Y"].values * (arr.shape[0] / orig_h)
    d = np.sqrt((sx[:, None] - sx[None, :]) ** 2 + (sy[:, None] - sy[None, :]) ** 2)
    np.fill_diagonal(d, np.inf)
    pitch = float(np.median(d.min(axis=1)))
    half = max(24, int(round(0.375 * pitch)))
    s1 = np.zeros(3); s2 = np.zeros(3); npix = 0
    for x, y in zip(sx, sy):
        x0, x1 = int(round(x - half)), int(round(x + half))
        y0, y1 = int(round(y - half)), int(round(y + half))
        x0, y0 = max(0, x0), max(0, y0)
        x1, y1 = min(arr.shape[1], x1), min(arr.shape[0], y1)
        if x1 <= x0 or y1 <= y0: continue
        t = np.asarray(Image.fromarray(arr[y0:y1, x0:x1]).resize((224, 224), Image.BILINEAR), np.float32) / 255.0
        lab = rgb2lab(t).reshape(-1, 3)
        s1 += lab.sum(axis=0); s2 += (lab ** 2).sum(axis=0); npix += lab.shape[0]
    mu_lab = s1 / npix
    sd_lab = np.sqrt(np.maximum(s2 / npix - mu_lab ** 2, 1e-6))
    rs = {"reference_section": "BC23901_C2", "n_tiles": int(npix // (224 * 224)),
          "tile_px": 224, "half_px": half, "pitch_px": pitch,
          "lab_mean": mu_lab.tolist(), "lab_std": sd_lab.tolist(),
          "space": "CIELAB via skimage.color.rgb2lab (L in [0,100])",
          "note": "Reinhard normalization reference; computed from internal tiles only, never external"}
    json.dump(rs, open(REF, "w"), indent=1)
    print("reference stats:", rs["n_tiles"], "tiles", flush=True)
REF_MU = np.array(rs["lab_mean"], np.float32); REF_SD = np.array(rs["lab_std"], np.float32)

def reinhard(tile01):
    lab = rgb2lab(tile01).reshape(-1, 3)
    m = lab.mean(axis=0); s = lab.std(axis=0); s[s < 1e-3] = 1e-3
    out = (lab - m) / s * REF_SD + REF_MU
    out[:, 0] = np.clip(out[:, 0], 0, 100); out[:, 1:] = np.clip(out[:, 1:], -128, 127)
    return np.clip(lab2rgb(out.reshape(224, 224, 3)), 0, 1).astype(np.float32)

# --- frozen model (loaded first; fail fast if artifact missing)
fz = np.load("results/g2x_frozen_model.npz", allow_pickle=True)
bank, Pnames = fz["bank"], list(fz["Pnames"])
universe = [str(g) for g in fz["universe"]]
sc_mu, sc_sd, mu, sd = fz["sc_mean"], fz["sc_std"], fz["mu"], fz["sd"]
counts = fz["pw_idx_counts"]; flat = fz["pw_idx_flat"]
pw_idx, o = {}, 0
for n, c in zip(Pnames, counts):
    pw_idx[n] = list(flat[o:o + c]); o += c

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
print(f"in-tissue spots: {len(bcs)} | pitch {pitch:.1f}px | half {half:.1f}px", flush=True)

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

start = 0
if os.path.exists(CKPT):
    z = np.load(CKPT); Xs, done_bcs = list(z["Xs"]), list(z["bcs"]); start = len(done_bcs)
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
    im = reinhard(im)  # THE single v2 change (prespec def3b27)
    t = torch.from_numpy(((im - MEAN) / STD).transpose(2, 0, 1)[None]).float()
    with torch.no_grad():
        ff = model.forward_features(t)
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
N = N[[bpos[b] for b in done_bcs]]

pw_true = {}
for n, idx in pw_idx.items():
    Z = N[:, idx].toarray(); Z = (Z - mu[idx]) / sd[idx]
    pw_true[n] = Z.mean(axis=1)
Z = (X - sc_mu) / sc_sd
Pp = models.predict_bank(Z, bank)
per_pw = {}
for k, n in enumerate(Pnames):
    t, p = pw_true[n], Pp[:, k]
    per_pw[n] = float(np.corrcoef(p, t)[0, 1]) if t.std() > 0 and p.std() > 0 else float("nan")
vals = np.array(list(per_pw.values()), float)
med = float(np.nanmedian(vals))
verdict = "PASS" if med >= 0.10 else "FAIL"
res = {"prespec": "results/g2x_v2_prespec.json (def3b27); lineage: v1 FAIL 34afa2f stays on record",
       "single_change": "Reinhard CIELAB stain normalization of external tiles to internal BC23901_C2 reference (results/g2x_v2_refstain.json)",
       "gate": "median pathway Pearson >= 0.10 (identical to v1)",
       "n_spots": int(len(done_bcs)), "n_pathways": len(Pnames),
       "median_pathway_r": med, "per_pathway_r": per_pw,
       "frac_pathways_r_pos": float(np.mean(vals > 0)),
       "verdict": verdict,
       "v1_comparison": {"median_pathway_r": -0.1547, "frac_pos": 0.16},
       "declared_gaps": "single external section; one normalization method; v2 designed after v1 FAIL - both reported with lineage"}
json.dump(res, open(OUT, "w"), indent=1)
print(f"G2-X v2 median pathway r = {med:.4f} -> {verdict} (v1 -0.1547) (wall {time.time()-t0:.0f}s)", flush=True)
