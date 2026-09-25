# PREREG - MEGA27 item 17-spatial: MorphoScan
Locked 2026-09-25 before any model result is seen. Gates are unchangeable; negatives are preserved.

## Study
MorphoScan: a spatial-transcriptomics + medical-image scanning toolkit.
Track A: predict spot-level spatial gene expression from H&E morphology.
Track B: scan H&E sections and flag tumor regions from morphology alone.

## Data
Primary: Mendeley Data 29ntw7sh4r v5, "Human breast cancer in situ capturing
transcriptomics" (68 sections, 23 patients, LumA/LumB/HER2+/TNBC, ST array,
H&E images + spot count matrices + coordinates + tumor annotations).
HER2+ subset used for benchmark comparability with ST-Net (He et al. 2020,
Nat Biomed Eng; published median per-gene Pearson ~0.19, top-250 genes,
leave-one-patient-out). Non-HER2 patients = frozen external transport cohort.

## Locked protocol (Track A)
- Genes: top 250 by mean normalized count, selected on TRAIN folds only.
- Normalization: library-size to 1e4 + log1p, parameters from train only.
- Features: handcrafted morphology (color deconvolution stats, GLCM texture,
  entropy, color moments) from a fixed spot-centered patch. No deep net.
- Split: leave-one-patient-out on HER2+ patients; no spot from a test patient
  in any fit (scaler, genes, ridge lambda by inner CV on train patients).
- Models: M0 train-mean; M1 ridge probe; M2 = M1 + Gaussian-kernel spatial
  graph smoothing (MorphoScan advance); M3 MLP (nonlinear check).

## Locked gates (Track A)
- G-A1 sanity: M1 median per-gene Pearson on held-out patients > M0 with
  permutation p < 0.05.
- G-A2 method advance: M2 beats M1 median per-gene Pearson by >= +0.005 with
  paired Wilcoxon p < 0.05 across the 250 genes. Else honest negative, M2 kept.
- G-A3 benchmark: report median Pearson vs ST-Net published ~0.19 (and
  HisToGene-class ~0.21). PASS = match or beat 0.19; gap logged honestly.
- G-A4 transport: train HER2+, test LumA/LumB/TNBC patients. Log degradation
  honestly; expected negative is documented, not re-fished.
- G-A5 coherence: top-25 predictable genes checked for pathway enrichment
  (g:Profiler/Reactome). Reported either way.

## Locked gates (Track B)
- G-B1: tumor-vs-other patch classifier, patient-held-out AUC >= 0.80.
  Else gap logged.
- G-B2: report calibration (Brier) and thresholded sensitivity/specificity.

## Honesty rules
No page padding, no fabricated counts, no unearned benchmark claims.
Every ledger row carries level-of-use (analyzed / metadata-mined / attempted).
Tool certifications require a committed evidence artifact.
