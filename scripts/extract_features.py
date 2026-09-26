"""Extract per-spot morphology features + expression + labels for all 68 sections.
Caches per-section npz into features_cache/ (gitignored). Images are decoded at
<=MAXDIM on the long edge via JPEG draft mode; spot pixel coords are rescaled
accordingly. Patch half-width = 0.375 x median nearest-neighbor pitch so a patch
covers roughly the 100um spot (ST-Net used a 224px window at this resolution)."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__) + "/..")
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from src.morphoscan import dataio, features

DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
OUT = "features_cache2"  # corrected-geometry rerun; features_cache/ (buggy) kept untouched for provenance
os.makedirs(OUT, exist_ok=True)
MAXDIM = 4500
Image.MAX_IMAGE_PIXELS = None

meta = dataio.read_metadata(os.path.join(DATA, "metadata.csv"))
log = open("/tmp/extract2.log", "a")

for _, row in meta.iterrows():
    sec = row["count_matrix"].replace("_stdata.tsv.gz", "")
    t0 = time.time()
    try:
        outp = os.path.join(OUT, f"{sec}.npz")
        if os.path.exists(outp):
            continue
        spots = dataio.read_spot_positions(os.path.join(DATA, row["spot_coordinates"]))
        spot_ids, gene_ids, counts = dataio.read_stdata(os.path.join(DATA, row["count_matrix"]))
        counts_df_ids = set(spot_ids)
        spots = spots[spots.index.astype(str).isin(counts_df_ids)]
        # tumor labels: Coords file xcoord/ycoord ~ array coords + noise
        tum = dataio.read_tumor_coords(os.path.join(DATA, row["tumor_annotation"]))
        tum["ax"] = tum["xcoord"].round().astype(int)
        tum["ay"] = tum["ycoord"].round().astype(int)
        lab_map = {f"{r.ax}x{r.ay}": r.tumor for r in tum.itertuples()}
        img = Image.open(os.path.join(DATA, row["histology_image"]))
        orig_w, orig_h = img.size  # header size BEFORE draft decode
        img.draft("RGB", (MAXDIM, MAXDIM))
        img = img.convert("RGB")
        if max(img.size) > MAXDIM:
            r = MAXDIM / max(img.size)
            img = img.resize((int(img.size[0] * r), int(img.size[1] * r)), Image.LANCZOS)
        arr = np.asarray(img)
        # PHASE-1 BUG FIX (PHASE1-GEOMETRY-BUG.md): exact header-derived scale.
        # Old code used arr.width/(max(spotX)*1.02), misplacing the grid and
        # NaN-ing ~82% of spots on affected sections.
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
            X[k, :] = features.patch_vector(patch)
        C = counts[keep_idx, :]
        labels = np.array([1 if lab_map.get(s) == "tumor" else 0 for s in spots.index.astype(str)], np.int8)
        lab_found = sum(1 for s in spots.index.astype(str) if s in lab_map)
        np.savez_compressed(outp, X=X, counts_sparse=sparse.csr_matrix(C),
                            spot_ids=np.array(spots.index.astype(str)),
                            px=np.stack([sx, sy], axis=1), labels=labels,
                            genes=np.array(gene_ids), pitch=np.array([pitch]))
        log.write(f"OK {sec} spots={len(spots)} pitch={pitch:.0f} half={half} labmatch={lab_found}/{len(spots)} "
                  f"tum_frac={labels.mean():.2f} {time.time()-t0:.0f}s\n"); log.flush()
        del img, arr, X, C
        import gc; gc.collect()
    except Exception as e:
        log.write(f"FAIL {sec} {e!r}\n"); log.flush()
log.write("DONE\n"); log.close()
