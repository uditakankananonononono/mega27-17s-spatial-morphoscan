# PREREG-2 - MorphoScan Phase 2: foundation embeddings + pathway-level prediction

Locked 2026-09-26 16:37 IST BEFORE any Phase-2 model result exists. PREREG.md is
untouched: Phase-1 FAILs (G-A2, G-A3, G-B1) stand as documented negatives and are
not re-fished. Phase 2 is new locked work under the 2026-09-26 standing rules
(WhatsApp 4:11:18 / 4:11:49 / 4:12:25 / 4:14:37): benchmark beat required, new
discovery required, negatives never terminal, redirection via ChatGPT
(JUDGE_ROUNDS.md round 1, verbatim).

## What changes vs Phase 1

- Representation: FROZEN CTransPath embeddings (Swin-Tiny, 768-dim; weights
  publicly released, github.com/Xiyue-Wang/TransPath, Medical Image Analysis
  2022) on spot-centered tiles resized to 224x224. No fine-tuning. LOCKED
  CONDITIONAL PIVOT: if G2-R FAILs, the representation (and only the
  representation) may switch to UNI or CONCH (HF-gated; needs account approval)
  before any other change.
- Targets: (a) MSigDB Hallmark-50 pathway scores per spot = mean of z-scored
  member-gene log-normalized expression, z-parameters from TRAIN folds only;
  (b) the same top-250 per-gene task as Phase 1, for direct comparability.
- Identical cohort, LOPO split, normalization, and no-leak discipline as
  PREREG.md. External cohort: a public no-auth 10x Visium human breast cancer
  dataset; the final all-HER2+ model is frozen before touching it.

## Locked gates (Phase 2)

- G2-R representation: CTransPath+ridge median held-out per-gene Pearson
  >= 0.08 AND beats Phase-1 handcrafted M1 in >= 80% of the 250 genes
  (paired Wilcoxon p < 0.05).
- G2-P pathway: median held-out pathway Pearson >= 0.20 across Hallmark-50,
  permutation p < 0.01; immune-program AUROC (interferon-alpha/gamma,
  inflammatory response, complement, allograft rejection; high vs low
  activity) >= 0.75.
- G2-B benchmark: median per-gene Pearson vs ST-Net published ~0.19 on the
  same cohort; PASS = match/beat 0.19, gap logged honestly either way.
  Locked pivot ladder on FAIL: (1) spatial-graph aggregation over embeddings,
  (2) multi-scale tiles, (3) UNI/CONCH representation.
- G2-X external: frozen model on the external Visium breast cohort, no
  retraining: median pathway Pearson >= 0.10.
- G2-A ablations: shuffled spot-coordinate control at chance for any spatial
  variant; shuffled-embedding control at chance for all predictors.
- G2-S scanner: CTransPath-embedding tumor-vs-other classifier,
  patient-held-out AUC >= 0.85 (Phase-1 G-B1 failed at 0.642-0.676; this is
  the new scanner gate, not a re-fish).

## Locked discovery targets (FDR < 0.05, train-only selection)

- D1: a visually-predictable-gene set (>= 50 held-out genes at Pearson
  >= 0.15) whose pathway enrichment replicates Phase-1's antigen-processing
  signal on held-out patients.
- D2: a spatial pathway-activity architecture (e.g., immune/stromal hotspot
  colocalization with tumor annotation) that is non-random against the
  shuffled-coordinate null.

## Honesty rules

Same as PREREG.md: no padding, no fabricated counts, no unearned benchmark
claims, every ledger row carries level-of-use, negatives preserved in the same
voice as passes.
