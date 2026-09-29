"""Tier-2 #7 technical-artifact control: Reinhard-normalized feature extraction.
PRE-DECLARED here before results: identical geometry/labels/counts as
scripts/extract_features.py (features_cache2), single change = each spot patch
is Reinhard-normalized in CIELAB to the internal BC23901_C2 reference
(results/g2x_v2_refstain.json, computed from internal tiles only) before the
handcrafted feature vector is computed. Output features_cache3/. Diagnostic
control, no gate; comparison lives in scripts/p1_norm_control.py.
Run: python3 scripts/extract_features_norm.py"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from skimage.color import rgb2lab, lab2rgb
from src.morphoscan import dataio, features

DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
OUT = "features_cache3"
os.makedirs(OUT, exist_ok=True)
MAXDIM = 4500
Image.MAX_IMAGE_PIXELS = None
rs = json.load(open("results/g2x_v2_refstain.json"))
REF_MU = np.array(rs["lab_mean"], np.float32); REF_SD = np.array(rs["lab_std"], np.float32)

def reinhard(patch):
    t = np.asarray(patch, np.float32) / 255.0
    lab = rgb2lab(t).reshape(-1, 3)
    m = lab.mean(axis=0); s = lab.std(axis=0); s[s < 1e-3] = 1e-3
    out = (lab - m) / s * REF_SD + REF_MU
    out[:, 0] = np.clip(out[:, 0], 0, 100); out[:, 1:] = np.clip(out[:, 1:], -128, 127)
    return (np.clip(lab2rgb(out.reshape(patch.shape)), 0, 1) * 255).astype(np.uint8)

meta = dataio.read_metadata(os.path.join(DATA, "metadata.csv"))
log = open("/tmp/extract3.log", "a")
for _, row in meta.iterrows():
    sec = row["count_matrix"].replace("_stdata.tsv.gz", "")
    t0 = time.time()
    try:
        outp = os.path.join(OUT, f"{sec}.npz")
        if os.path.exists(outp):
            continue
        spots = dataio.read_spot_positions(os.path.join(DATA, row["spot_coordinates"]))
        spot_ids, gene_ids, counts = dataio.read_stdata(os.path.join(DATA, row["count_matrix"]))
        spots = spots[spots.index.astype(str).isin(set(spot_ids))]
        tum = dataio.read_tumor_coords(os.path.join(DATA, row["tumor_annotation"]))
        tum["ax"] = tum["xcoord"].round().astype(int)
        tum["ay"] = tum["ycoord"].round().astype(int)
        lab_map = {f"{r.ax}x{r.ay}": r.tumor for r in tum.itertuples()}
        img = Image.open(os.path.join(DATA, row["histology_image"]))
        orig_w, orig_h = img.size
        img.draft("RGB", (MAXDIM, MAXDIM))
        img = img.convert("RGB")
        if max(img.size) > MAXDIM:
            r = MAXDIM / max(img.size)
            img = img.resize((int(img.size[0] * r), int(img.size[1] * r)), Image.LANCZOS)
        arr = np.asarray(img)
        xcol = "X" if "X" in spots.columns else "x"
        ycol = "Y" if "Y" in spots.columns else "y"
        sx = (spots[xcol].values * (arr.shape[1] / orig_w)).astype(float)
        sy = (spots[ycol].values * (arr.shape[0] / orig_h)).astype(float)
        d = np.sqrt((sx[:, None] - sx[None, :]) ** 2 + (sy[:, None] - sy[None, :]) ** 2)
        np.fill_diagonal(d, np.inf)
        pitch = np.median(d.min(axis=1))
        half = max(24, int(round(0.375 * pitch)))
        keep_idx = [spot_ids.index(s) for s in spots.index.astype(str)]
        X = np.zeros((len(spots), len(features.feature_names())), np.float32)
        for k, (cx, cy) in enumerate(zip(sx, sy)):
            patch = dataio.extract_patch(arr, cx, cy, half)
            if patch.shape[0] < 8 or patch.shape[1] < 8:
                X[k, :] = np.nan
                continue
            X[k, :] = features.patch_vector(reinhard(patch))
        C = counts[keep_idx, :]
        labels = np.array([1 if lab_map.get(s) == "tumor" else 0 for s in spots.index.astype(str)], np.int8)
        np.savez_compressed(outp, X=X, counts_sparse=sparse.csr_matrix(C),
                            spot_ids=np.array(spots.index.astype(str)),
                            px=np.stack([sx, sy], axis=1), labels=labels,
                            genes=np.array(gene_ids), pitch=np.array([pitch]))
        log.write(f"OK {sec} spots={len(spots)} {time.time()-t0:.0f}s\n"); log.flush()
        del img, arr, X, C
        import gc; gc.collect()
    except Exception as e:
        log.write(f"FAIL {sec} {e!r}\n"); log.flush()
log.write("DONE\n"); log.close()
