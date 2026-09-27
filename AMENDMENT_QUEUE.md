# AMENDMENT_QUEUE.md - Lane B (mega27-17s), LOCKED from ChatGPT judge verdict 2026-09-27

Locked from JUDGE_VERDICT_CHATGPT_20260927.md (conversation 6ab8c1ed-f308-83e8-b0fd-5e9144ebbd9f,
gate MET entry in JUDGE_ROUNDS.md, commit e5b1456). Protocol: queue locked BEFORE execution;
each item closes only with committed, live-verified evidence (script + output + hash).
Sub-bullets truncated below are verbatim in the archived verdict file.

Status: DONE = evidence committed. PARTIAL = exists, gap named. TODO = not started.
BLOCKED marks items needing inputs destroyed in the 2026-09-27 sandbox wipe
(features_cache/, embeddings_cache/, data/her2st/, CTransPath weights - re-acquisition tracked separately).

## Tier 1 - cheap, high rigor value (no wiped inputs needed)
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 2 | Permutation resolution 50 -> 10k (or Monte Carlo stopping rule) | TODO | Rerun G2-P permutation at 10k; report resolution-corrected p. |
| 1 | Bootstrap CIs for all headline metrics (median r, fold medians, M2 delta, scanner AUC) | TODO | 10k bootstrap over committed per-gene/per-fold results. |
| 5 | Per-gene confidence intervals | TODO | Ships with #1. |
| 13 | Label-permutation scanner control (Track B AUC -> ~0.5) | TODO | Permute tumor labels within patients, confirm collapse. |
| 17 | Scanner feature-importance stability (bootstrap coefficients) | TODO | With #1 bootstrap machinery. |
| 9 | Multi-task gene prediction: elastic-net multi-output, PLS, reduced-rank vs independent ridge | TODO | sklearn over committed features where available; blocked portions noted. |
| 11 | Nested evaluation of the 250-gene selection rule | TODO | Sensitivity of results to the selection rule. |
| 12 | Negative-control gene experiment (gene-level) | PARTIAL | PancreasBeta negative control exists at niche level (D2). Gap: gene-level negative controls' correlation distribution. |
| 3 | Complete spatial null model for M2 (randomized graphs, shuffled coords, distance-preserving perms) | TODO | Uses committed geometry; no image inputs needed. |
| 4 | Simple spatial baselines: NN prediction, Gaussian smoothing, section-mean interpolation, Moran-based predictor | TODO | As #3. |

## Tier 2 - model/spec completeness (paper + sklearn-level compute)
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 8 | Full MLP specification + result table (layers, units, optimizer, epochs, params) | TODO | Paper-level; extract from committed M3 code. |
| 14 | Scanner subgroup analyses (AUC by subgroups per verdict) | TODO | From committed scanner outputs. |
| 15 | Scanner calibration diagnostics | TODO | Reliability curve, Brier, ECE. |
| 18 | Full pipeline reproducibility artifacts | PARTIAL | Repo scripts/results committed. Gap: per verdict sub-items (env lock, run manifest). |
| 19 | Computational cost benchmarking | TODO | Rerun timing on committed pipelines; paper table. |

## Tier 3 - BLOCKED on sandbox-wipe re-acquisition
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 6 | Morphology ablations (stain/texture/color/entropy/nuclear-density only) | TODO-BLOCKED | Needs features_cache (wiped). |
| 7 | Technical artifact controls (Macenko/Reinhard normalization) | TODO-BLOCKED | Needs raw HER2ST images (wiped). |
| 10 | Learned image baseline, same patient-held-out split | TODO-BLOCKED | Needs CTransPath weights + embeddings (wiped). |
| 16 | External cohort validation, frozen model, no retraining | TODO-BLOCKED | = PREREG-2 G2-X gate. Data prep spec done pre-wipe (10x Parent_Visium_Human_BreastCancer URLs live-verified); needs re-download + embeddings rebuild. |
| 20 | Challenge-style hidden evaluation protocol | TODO-CONSTRAINED | Extends #16; freeze code+model, evaluate on unseen public sections. |

Execution order: 2 -> 1+5 -> 13 -> 17 -> 3 -> 4 -> 9 -> 11 -> 12 -> 14,15,8,18,19 -> (rebuild) 6,7,10 -> 16 -> 20.
Every closure lands as: script + results JSON + tests where fitting + paper section increment (50+pp directive, user 12:02:56), committed with live-verified numbers.
