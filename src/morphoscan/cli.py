"""MorphoScan command-line tool.
  morphoscan scan --image IMG --coords CSV --model MODEL.npz   tumor-region scan
  morphoscan predict --image IMG --coords CSV --bank BANK.npz  expression predict
  morphoscan features --image IMG --x X --y Y --half H         one-patch features
"""
from __future__ import annotations
import argparse
import numpy as np
from PIL import Image
from . import features, metrics, dataio


def _load_image(path, max_dim=8000):
    img = Image.open(path)
    img.draft("RGB", (max_dim, max_dim))
    return np.asarray(img.convert("RGB"))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="morphoscan")
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("features")
    f.add_argument("--image", required=True)
    f.add_argument("--x", type=float, required=True)
    f.add_argument("--y", type=float, required=True)
    f.add_argument("--half", type=int, default=112)

    s = sub.add_parser("scan")
    s.add_argument("--image", required=True)
    s.add_argument("--coords", required=True, help="CSV with x,y columns (pixels)")
    s.add_argument("--model", required=True, help="npz with coef,intercept,mu,sd")
    s.add_argument("--half", type=int, default=112)
    s.add_argument("--out", required=True)

    args = ap.parse_args(argv)
    if args.cmd == "features":
        img = _load_image(args.image)
        v = features.patch_vector(dataio.extract_patch(img, args.x, args.y, args.half))
        for n, val in zip(features.feature_names(), v):
            print(f"{n}\t{val:.6g}")
    elif args.cmd == "scan":
        img = _load_image(args.image)
        m = np.load(args.model)
        import pandas as pd
        df = pd.read_csv(args.coords)
        df.columns = [c.strip().lower() for c in df.columns]
        X = np.stack([features.patch_vector(dataio.extract_patch(img, r.x, r.y, args.half))
                      for r in df.itertuples()])
        Z = (X - m["mu"]) / m["sd"]
        logit = Z @ m["coef"] + m["intercept"]
        prob = 1.0 / (1.0 + np.exp(-logit))
        out = df.copy(); out["tumor_probability"] = prob
        out.to_csv(args.out, index=False)
        print(f"wrote {args.out}: {len(out)} spots, mean p={prob.mean():.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
