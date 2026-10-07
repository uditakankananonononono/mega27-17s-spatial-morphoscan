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


def test_read_tumor_coords_cr_terminated(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_bytes(b"hdr\tx\ty\tl\tt\rs1\t1.5\t2.5\tA\ttumor\rs2\t3\t4\tB\tnon\r")
    df = D.read_tumor_coords(str(p))
    assert list(df.columns) == ["spot_key", "xcoord", "ycoord", "lab", "tumor"]
    assert len(df) == 2 and df.loc[0, "spot_key"] == "s1" and df.loc[1, "tumor"] == "non"
    assert df.loc[0, "xcoord"] == 1.5
