import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from src.morphoscan import features as F


def white(v=255):
    return np.full((32, 32, 3), v, dtype=np.uint8)

def test_od_white_is_zero():
    od = F.optical_density(white(255), eps=0.0)
    assert np.allclose(od, 0.0)

def test_od_dark_positive():
    assert F.optical_density(white(10)).mean() > 3.0

def test_deconv_shape():
    assert F.color_deconvolution(white()).shape == (32 * 32, 3)

def test_deconv_white_near_zero_concentration():
    C = F.color_deconvolution(white(255))
    assert abs(C.mean()) < 0.02  # eps=1.0 makes white slightly negative OD

def test_glcm_constant_image_energy_one():
    m = F.glcm(np.zeros((16, 16)), levels=4)[0]
    _, _, energy, _ = F.glcm_props(m)
    assert abs(energy - 1.0) < 1e-9

def test_glcm_contrast_zero_constant():
    m = F.glcm(np.full((16, 16), 7), levels=8)[0]
    contrast, homogeneity, _, _ = F.glcm_props(m)
    assert contrast == 0.0 and homogeneity > 0.99

def test_entropy_uniform_max():
    assert abs(F.shannon_entropy(np.ones(8) / 8) - 3.0) < 1e-9

def test_entropy_point_mass_zero():
    p = np.zeros(8); p[0] = 1.0
    assert F.shannon_entropy(p) == 0.0

def test_feature_vector_length_matches_names():
    img = np.random.default_rng(0).integers(0, 255, (48, 48, 3), dtype=np.uint8)
    v = F.patch_vector(img)
    assert v.shape[0] == len(F.feature_names())
    assert np.all(np.isfinite(v))

def test_feature_names_unique():
    n = F.feature_names()
    assert len(n) == len(set(n))

def test_textured_vs_flat_features_differ():
    flat = np.full((48, 48, 3), 200, dtype=np.uint8)
    rng = np.random.default_rng(3)
    tex = rng.integers(0, 255, (48, 48, 3), dtype=np.uint8)
    vf, vt = F.patch_vector(flat), F.patch_vector(tex)
    assert not np.allclose(vf, vt)
