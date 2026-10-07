# MorphoScan — MEGA27 item 17-spatial

Pre-registered prediction of spatial gene expression from H&E morphology, plus a
tumor-region scanner, in breast cancer — with an honest benchmark against
published deep-learning baselines.

MorphoScan asks a narrow, falsifiable question: how much spatial expression
signal is recoverable from hematoxylin-and-eosin morphology alone, with a fully
interpretable feature set, and where does that honest number sit relative to
published deep-learning systems? Every gate was locked in `PREREG.md` before
any result existed. All negatives are preserved and reported in the same voice
as passes.

## Cohort

HER2-ST breast-cancer spatial transcriptomics cohort (Mendeley Data
`29ntw7sh4r` v5, doi:10.17632/29ntw7sh4r.5): 68 sections / 23 patients, all
files SHA-256 verified. Evaluation cohort: 10 HER2+ patients / 29 sections
(leave-one-patient-out); 13 non-HER2 patients / 39 sections frozen for the
cross-subtype transport test. Gene universe: 11,787-gene intersection.

## Locked gates — outcomes exactly as measured

**Superseded-measurement notice (2026-10-07 audit).** The table below is the Phase-1 record measured with a spot-to-image geometry mapping that placed most spots outside the decoded image (e.g. 210 of 256 spots NaN on BC23287_C1; see `PHASE1-GEOMETRY-BUG.md`). Those exact table values (Track A M1 0.014 over 280 genes; scanner AUC 0.642/0.676) match the buggy-geometry files in `results/phase1_buggy_geometry/` and are kept as the audit trail. They are NOT the current results. The locked gates in `PREREG.md` are unchanged.

Corrected-geometry reruns exist (the "rerun queued" line in `PHASE1-GEOMETRY-BUG.md` is stale):
- Track A rerun (commit 123b892, `results/trackA_results.json`): median held-out Pearson M1 0.0382 / M2 0.0408 across 278 evaluated genes, permutation p = 0.0196. The M2 - M1 gap is about +0.0025, still under the +0.005 G-A2 threshold. This file, not the 0.014 value, is the final Phase-1 Track A record.
- Track B v2 (commit fe5b87f, `results/trackB_results_v2.json`): median AUC 0.892 (LR) / 0.877 (GB), sensitivity 0.807, specificity 0.774, Brier 0.152. A formal re-verdict of G-B1 on v2 is not recorded in this README.
- Phase-2 pipeline (`results/p2_gate_verdicts.json`): median per-gene Pearson 0.0521 over 275 genes (G2-R and G2-P not passed; gap to ST-Net 0.19 is -0.138, characterization only) and scanner median AUC 0.926 (G2-S pass).

Read the table as the superseded Phase-1 record.

| Gate | Threshold (locked) | Measured | Outcome |
|---|---|---|---|
| G-A1 sanity | M1 > M0, permutation p < 0.05 | M0 0.000, M1 0.014, p = 0.020 | **PASS** |
| G-A2 method advance | M2 − M1 ≥ +0.005, Wilcoxon p < 0.05 | +0.0029, p = 3.03e-07 | **FAIL** (honest negative; M2 kept) |
| G-A3 benchmark | match/beat ST-Net median r ≈ 0.19 | 0.014 (HisToGene-class ≈ 0.21) | **FAIL** (gap logged, not hidden) |
| G-A4 transport | log cross-subtype degradation, no re-fishing | LumA 0.036 / LumB 0.015 / TNBC 0.008 | logged |
| G-A5 coherence | pathway enrichment of top-25 predictable genes | antigen processing (KEGG p = 9.2e-4); GRB7/ERBB2 signaling (REAC p = 6.7e-3) | **PASS** |
| G-B1 scanner | patient-held-out AUC ≥ 0.80 | 0.642 (LR) / 0.676 (GB) | **FAIL** (honest negative) |
| G-B2 calibration | report Brier + sens/spec | Brier 0.189, sens 0.719, spec 0.725 | **PASS** |

Track A runs 10/10 LOPO folds: median held-out per-gene Pearson M0 0.000 /
M1 0.014 / M2 0.016 across 280 evaluated genes. The comparison to deep
baselines is reported as measured: 41 handcrafted morphology features and a
linear ridge bank vs a fine-tuned DenseNet-121 — the probe is linear and
trains in seconds on a CPU; where it lands relative to those numbers is in
the paper.

## Audit gates

- **Tools:** 51 certified external tools/services (evidence artifacts in
  `results/tools/`; failed attempts logged as attempted-failed, never counted).
- **Datasets:** 629 accession-level records (`results/accession_ledger.csv`),
  level-of-use disclosed per row (69 analyzed-primary + 560 metadata-mined).
- **Tests:** 50 unit tests (`pytest tests/`; re-run 2026-10-07: 50 passed). The earlier "32" figure was stale.

## Paper

Research paper (currently 51 pages per pdfinfo on `paper/main.pdf`, CI-rendered with xelatex per `paper/RENDER_STATUS.md`; the earlier 21-page / lualatex description is stale): `paper/main.pdf` in this
repo, and the Drive copy:
https://drive.google.com/file/d/1U6AETHusYxMZxdN4loNa8dpNA5z-u90N/view

## CLI

```bash
PYTHONPATH=src python3 -c "from morphoscan import cli; cli.main(['scan', \
  '--image', 'HE_BT23287_C1.jpg', \
  '--coords', 'spots_BT23287_C1.csv.gz', \
  '--model', 'results/scanner_model.npz', \
  '--out', '/tmp/scan.csv'])"
# -> writes per-spot tumor probabilities (smoke-tested: 256 spots, mean p=0.164)
```

`morphoscan features --image IMG --x X --y Y` prints the 41-feature vector for
one patch. The CLI runs on any H&E JPEG plus a spot-coordinate CSV (x/y pixel
columns, any case).

## Layout

`PREREG.md` locked gates · `src/morphoscan/` package (metrics, features,
models, dataio, cli) · `scripts/` download, harvest, evidence, extraction,
Track A/B, ledgers, figures, paper tables · `tests/` unit tests ·
`results/` ledgers, per-gene table, scanner model, evidence artifacts ·
`figures/` · `paper/`. Randomness seeded (20260925).
