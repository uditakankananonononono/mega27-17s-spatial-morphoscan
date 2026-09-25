"""Loaders for Mendeley 29ntw7sh4r (breast cancer ST cohort)."""
from __future__ import annotations
import gzip, io, os
import numpy as np
import pandas as pd


def read_metadata(path):
    return pd.read_csv(path)


def read_stdata(path):
    """Returns (spot_ids list, gene_ids list, counts float32 matrix)."""
    df = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    return df.index.tolist(), df.columns.tolist(), df.values.astype(np.float32)


def read_spot_positions(path):
    """spots_BT*.csv.gz: pixel positions of spots in the histology image.
    Returns DataFrame with whatever columns the file carries."""
    return pd.read_csv(path, compression="gzip")


def read_tumor_coords(path):
    """*_Coords.tsv.gz: single-line-per-spot annotation with xcoord/ycoord/lab/tumor."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["spot_key", "xcoord", "ycoord", "lab", "tumor"])
    return df


def extract_patch(img_arr, cx, cy, half):
    """Crop a (2*half) square patch centered at (cx, cy), edge-clamped."""
    h, w = img_arr.shape[:2]
    x0, x1 = max(0, int(cx - half)), min(w, int(cx + half))
    y0, y1 = max(0, int(cy - half)), min(h, int(cy + half))
    return img_arr[y0:y1, x0:x1]
