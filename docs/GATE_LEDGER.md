# Phase-1 gate ledger (written 2026-10-07)

Locked gates are in PREREG.md (last modified 2026-09-25, never amended). PREREG-2.md states that the Phase-1 FAILs stand.

| Gate | Locked verdict | Corrected-geometry rerun (descriptive) |
|---|---|---|
| G-A2 (M2-M1 >= +0.005) | FAIL | Track A rerun (123b892): M2-M1 about +0.0025. Still FAIL. |
| G-A3 (match ST-Net median r ~0.19) | FAIL | Rerun per-gene median M1 0.0350, M2 0.0375. Still FAIL. |
| G-B1 (patient-held-out AUC >= 0.80) | FAIL (0.642 LR / 0.676 GB) | Track B v2 (fe5b87f): median AUC 0.892 LR / 0.877 GB. POST-HOC. See below. |

## G-B1 lineage (checked 2026-10-07)
- The locked gate was measured first (2026-09-25) and failed. The geometry bug was found and documented on 2026-09-26 (09ff7b2). The v2 script (4dd1bf3, identical recipe, features_cache2) was committed before its run on 2026-09-29, but no PREREG amendment names a verdict rule for v2, and PREREG.md has no later commit.
- Ruling recorded: v2 ran after the locked gate had already failed and under no pre-registered amendment, so it is post-hoc relative to G-B1. It cannot pass the locked gate. G-B1 stays FAIL. The v2 numbers are descriptive only.
- The scanner claim that stands on a pre-registered gate is PREREG-2 G2-S (median AUC 0.926, pass), which was prespecified before its run.

## Known open item
The paper's Track A/B numbers have not been re-audited against these corrected files in this pass.
