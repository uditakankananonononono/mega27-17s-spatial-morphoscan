import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd
from PIL import Image
from src.morphoscan import cli, features as F


def _img(tmp_path):
    rng = np.random.default_rng(0)
    p = tmp_path / "t.png"
    Image.fromarray(rng.integers(0, 255, (120, 120, 3), dtype=np.uint8)).save(p)
    return str(p)

def test_features_cmd_prints_all_names(tmp_path, capsys):
    assert cli.main(["features", "--image", _img(tmp_path), "--x", "60", "--y", "60", "--half", "20"]) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert [l.split("\t")[0] for l in lines] == F.feature_names()

def test_scan_matches_manual_logistic(tmp_path):
    img = _img(tmp_path)
    n = len(F.feature_names())
    rng = np.random.default_rng(1)
    coef = rng.normal(0, 0.1, n)
    mp = tmp_path / "m.npz"
    np.savez(mp, coef=coef, intercept=np.array(0.2), mu=np.zeros(n), sd=np.ones(n))
    cp = tmp_path / "c.csv"
    pd.DataFrame({" X ": [40, 60], "Y": [50, 70]}).to_csv(cp, index=False)
    out = tmp_path / "o.csv"
    assert cli.main(["scan", "--image", img, "--coords", str(cp), "--model", str(mp),
                     "--half", "20", "--out", str(out)]) == 0
    df = pd.read_csv(out)
    from src.morphoscan import dataio
    arr = np.asarray(Image.open(img).convert("RGB"))
    for r in df.itertuples():
        v = F.patch_vector(dataio.extract_patch(arr, r.x, r.y, 20))
        p = 1 / (1 + np.exp(-(v @ coef + 0.2)))
        assert abs(r.tumor_probability - p) < 1e-9
