# PHASE-1 GEOMETRY BUG - documented audit finding (2026-09-26 revival)

## What happened
Phase 1 (scripts/extract_features.py) mapped spot pixel coordinates to the
decoded histology image with `scale = arr.width / (max(spot_X) * 1.02)`. Spot
coordinates are in original-image pixel space; normalizing by the maximum spot
X misplaces the entire grid whenever the tissue does not span the full slide
width (the common case).

## Measured impact (verified 2026-09-26)
- BC23287_C1: 210 of 256 spots (82%) fall outside the decoded image and are
  NaN in Phase 1's OWN committed cache (features_cache/BC23287_C1.npz,
  nan rows = 210/256). The surviving 46 patches are at wrong positions.
- Phase-1 Track A/B gate outcomes (G-A1..G-A5, G-B1..G-B2) were therefore
  measured on a small, mis-positioned spot subset per section. The gates stand
  in PREREG.md; the measurements are superseded by corrected-geometry reruns
  (this is a bug fix with full disclosure, not re-fishing).

## The fix
Exact header-derived scale factors: read the JPEG header size BEFORE draft
decode, scale spot coords by decoded/original per axis. Verified: 100% of
spots in-bounds on all 68 sections, and a rendered spot-grid overlay sits
exactly on the tissue (visual check, spot grid covers tissue blob with no
spots on empty glass). Corrected pitch on BC23287_C1 = 139.2 px (Phase-1
cache carried 485.6 px from a different decode).

## Status
- Phase-2 pipeline (scripts/p2_extract_ctranspath.py) uses the corrected
  geometry from the start.
- Phase-1 corrected-geometry reruns are DONE (updated 2026-10-07; this line previously said queued): Track A 123b892 (results/trackA_results.json), Track B v2 fe5b87f (results/trackB_results_v2.json). Old results are preserved in results/phase1_buggy_geometry/. Gate verdicts: docs/GATE_LEDGER.md.
