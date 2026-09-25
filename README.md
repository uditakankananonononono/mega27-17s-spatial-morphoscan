# MorphoScan (MEGA27 item 17-spatial)

Spatial-transcriptomics + medical-image scanning toolkit: predicts spot-level
spatial gene expression from H&E morphology and scans breast cancer sections
for tumor regions. Fully pre-registered (PREREG.md), all gates locked before
results, negatives preserved.

## Cohort
Mendeley Data 29ntw7sh4r v5 - 68 sections / 23 breast cancer patients
(LumA/LumB/HER2+/TNBC), ST-array spot counts + H&E images + tumor annotations.
HER2+ patients: leave-one-patient-out benchmark vs published ST-Net/HisToGene.
Non-HER2 patients: frozen external transport cohort.

## Layout
- `PREREG.md` - locked gates (Track A: G-A1..G-A5; Track B: G-B1, G-B2)
- `src/morphoscan/` - package: metrics, features (color deconvolution, GLCM,
  entropy), models (vectorized ridge bank, spatial smoothing), dataio, CLI
- `scripts/` - harvest_accessions, tools_evidence(1,2), extract_features,
  run_trackA, run_trackB, build_ledgers, make_figures
- `tests/` - 32 unit tests (metrics + features against known values)
- `results/` - trackA/trackB results JSON, per-gene table, accession_ledger.csv,
  tools_ledger.csv, tools/ evidence artifacts, scanner_model.npz
- `paper/` - Times New Roman research paper (lualatex + fontspec)
- `figures/` - PDF figures used in the paper

## CLI
```
python3 -m src.morphoscan.cli features --image HE.jpg --x 5200 --y 1160 --half 112
python3 -m src.morphoscan.cli scan --image HE.jpg --coords spots.csv --model results/scanner_model.npz --out scan.csv
```

## Honesty conventions
- accession_ledger level_of_use: analyzed-primary / metadata-mined (disclosed)
- tools_ledger status: certified (committed evidence artifact) vs attempted-failed
- No benchmark claim beyond measured evidence; gaps vs published baselines logged.
