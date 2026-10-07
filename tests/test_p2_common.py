import sys, os
ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts"))
import numpy as np
from scipy import sparse
import p2_common as C


def test_pathway_idx_min5_and_symbol_map():
    uni = [f"E{i}" for i in range(8)]
    sym = {f"E{i}": f"S{i}" for i in range(8)}
    hall = {"big": ["S0", "S1", "S2", "S3", "S4", "ZZ"], "small": ["S5", "S6"]}
    pw = C.pathway_idx(uni, sym, hall)
    assert list(pw) == ["big"] and pw["big"] == [0, 1, 2, 3, 4]

def test_pathway_scores_known_mean_z():
    N = sparse.csr_matrix(np.array([[1.0, 3.0, 9.0], [3.0, 5.0, 9.0]]))
    mu = np.array([2.0, 4.0, 9.0]); sd = np.array([1.0, 1.0, 1.0])
    s = C.pathway_scores(N, mu, sd, {"p": [0, 1]})["p"]
    assert np.allclose(s, [-1.0, 1.0])

def test_build_universe_intersection(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "FC", str(tmp_path))
    np.savez(tmp_path / "a.npz", genes=np.array(["g1", "g2", "g3"]))
    np.savez(tmp_path / "b.npz", genes=np.array(["g2", "g3", "g4"]))
    assert C.build_universe(["a", "b"]) == ["g2", "g3"]
