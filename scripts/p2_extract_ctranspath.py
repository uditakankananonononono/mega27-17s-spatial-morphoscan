"""Phase 2: extract frozen CTransPath embeddings per spot, IDENTICAL patch
geometry to Phase 1 (extract_features.py): same image decode (MAXDIM 4500 JPEG
draft), same coord scaling, same half = 0.375 x median NN pitch. Only the
representation changes (handcrafted 41 -> frozen CTransPath 768-d).
Checkpointed per section; safe to re-run. Usage: python3 scripts/p2_extract_ctranspath.py [N sections max]"""
import os, sys, time, re
sys.path.insert(0, os.path.dirname(__file__) + "/..")
sys.path.insert(0, "/home/sandbox/ctranspath_src")
import numpy as np
from PIL import Image
import torch, torch.nn as nn
from ctran import ConvStem
import timm
from src.morphoscan import dataio

DATA = "/home/sandbox/mega27-17s-spatial-morphoscan/data/her2st"
OUT = "/home/sandbox/mega27-17s-spatial-morphoscan/embeddings_cache"
os.makedirs(OUT, exist_ok=True)
MAXDIM = 4500
Image.MAX_IMAGE_PIXELS = None
MEAN = np.array([0.485, 0.456, 0.406], np.float32)
STD = np.array([0.229, 0.224, 0.225], np.float32)

torch.set_num_threads(int(os.environ.get("P2_THREADS", "2")))
model = timm.create_model("swin_tiny_patch4_window7_224", pretrained=False)
model.patch_embed = ConvStem(embed_dim=96, norm_layer=torch.nn.LayerNorm)
# Post-wipe reconstruction (2026-09-27): timm 0.5.4 + jamesdolezal/CTransPath checkpoint.
# ConvStem built directly (timm 0.5.4 create_model ignores embed_layer kwarg);
# embed_dim=96 + LayerNorm matches checkpoint patch_embed keys exactly; NO downsample
# remap needed with this timm/checkpoint pair (the old +1 remap mismatched here).
# Verified: missing == head.* only, zero unexpected keys, forward pass OK.
# mmap=False + clone: default torch.load mmaps the checkpoint, leaving weights
# file-backed; under memory pressure they get evicted and re-faulted from the
# REMOTE home mount every forward pass (observed 10x inference slowdown).
sd = torch.load("/home/sandbox/models/ctranspath.pth", map_location="cpu", weights_only=False, mmap=False)
sd = sd["model"] if "model" in sd else sd
sd = {k: v.clone() for k, v in sd.items()}
missing, unexpected = model.load_state_dict(sd, strict=False)
del sd
assert all(m.startswith("head.") for m in missing) and not unexpected, (missing, unexpected)
model.eval()

meta = dataio.read_metadata(os.path.join(DATA, "metadata.csv"))
max_sec = int(sys.argv[1]) if len(sys.argv) > 1 else 10**9
done = 0
log = open("/tmp/p2_extract.log", "a")

for _, row in meta.iterrows():
    sec = row["count_matrix"].replace("_stdata.tsv.gz", "")
    outp = os.path.join(OUT, f"{sec}.npz")
    if os.path.exists(outp) or done >= max_sec:
        continue
    t0 = time.time()
    try:
        spots = dataio.read_spot_positions(os.path.join(DATA, row["spot_coordinates"]))
        spot_ids, gene_ids, counts = dataio.read_stdata(os.path.join(DATA, row["count_matrix"]))
        spots = spots[spots.index.astype(str).isin(set(spot_ids))]
        xcol = "X" if "X" in spots.columns else "x"
        ycol = "Y" if "Y" in spots.columns else "y"
        img = Image.open(os.path.join(DATA, row["histology_image"]))
        orig_w, orig_h = img.size  # header size BEFORE draft decode
        img.draft("RGB", (MAXDIM, MAXDIM))
        img = img.convert("RGB")
        if max(img.size) > MAXDIM:
            r = MAXDIM / max(img.size)
            img = img.resize((int(img.size[0] * r), int(img.size[1] * r)), Image.LANCZOS)
        arr = np.asarray(img)
        t_img = time.time()
        # PHASE-1 BUG FIX: exact header-derived scale factors (spot pixel coords
        # are in original-image space). Phase 1 used arr.width/(max(spotX)*1.02),
        # which misplaces the grid and drops ~82% of spots (210/256 on
        # BC23287_C1, verified against Phase-1's own features_cache NaN rows).
        sx = spots[xcol].values * (arr.shape[1] / orig_w)
        sy = spots[ycol].values * (arr.shape[0] / orig_h)
        d = np.sqrt((sx[:, None] - sx[None, :]) ** 2 + (sy[:, None] - sy[None, :]) ** 2)
        np.fill_diagonal(d, np.inf)
        pitch = np.median(d.min(axis=1))
        half = max(24, int(round(0.375 * pitch)))
        keep = [spot_ids.index(s) for s in spots.index.astype(str)]
        X = np.zeros((len(spots), 768), np.float16)
        batch, bidx = [], []
        def flush():
            if not batch: return
            t = torch.from_numpy(np.stack(batch)).permute(0, 3, 1, 2)
            with torch.no_grad():
                f = model.forward_features(t).mean(dim=(1, 2))
            for j, k in enumerate(bidx):
                X[k] = f[j].numpy().astype(np.float16)
            batch.clear(); bidx.clear()
        for k, (cx, cy) in enumerate(zip(sx, sy)):
            p = dataio.extract_patch(arr, cx, cy, half)
            if p.shape[0] < 8 or p.shape[1] < 8:
                X[k] = np.nan; continue
            p = np.asarray(Image.fromarray(p).resize((224, 224), Image.BILINEAR), np.float32) / 255.0
            batch.append((p - MEAN) / STD); bidx.append(k)
            if len(batch) == 8: flush()
        flush()
        np.savez_compressed(outp, spot_ids=np.array(spots.index.astype(str)), X=X,
                            sx=sx, sy=sy, keep_idx=np.array(keep))
        done += 1
        msg = f"{sec}: {len(spots)} spots half={half} {time.time()-t0:.0f}s"
        msg = msg + f" [img {t_img - t0:.0f}s infer {time.time() - t_img:.0f}s]"
        print(msg, flush=True); log.write(msg + "\n"); log.flush()
    except Exception as e:
        log.write(f"{sec}: ERROR {e}\n"); log.flush()
print(f"sections done this run: {done}", flush=True)
