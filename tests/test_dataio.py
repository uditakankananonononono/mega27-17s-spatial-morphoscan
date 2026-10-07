import sys, os, gzip
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from src.morphoscan import dataio as D


def test_read_stdata_roundtrip(tmp_path):
    p = tmp_path / "s.tsv.gz"
    with gzip.open(p, "wt") as f:
        f.write("\tg1\tg2\n1x1\t1\t2\n1x2\t3\t4\n")
    spots, genes, m = D.read_stdata(str(p))
    assert spots == ["1x1", "1x2"] and genes == ["g1", "g2"]
    assert m.dtype == np.float32 and m.tolist() == [[1, 2], [3, 4]]

def test_read_spot_positions(tmp_path):
    p = tmp_path / "p.csv.gz"
    with gzip.open(p, "wt") as f:
        f.write("id,x,y\ns1,10,20\n")
    df = D.read_spot_positions(str(p))
    assert df.loc["s1", "x"] == 10

def test_read_metadata(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("a,b\n1,2\n")
    assert D.read_metadata(str(p)).shape == (1, 2)
