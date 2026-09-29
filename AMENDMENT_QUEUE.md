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
| 2 | Permutation resolution 50 -> 10k (or Monte Carlo stopping rule) | DONE | results/p2_perm10k_verdict.json: 10,000 perms (200 committed + 9,800 appended, Cholesky-shared fits verified identical estimator at 1.8e-05). Pathway real median r 0.3261, perm_p_10k = 0.0011 - G2-P survives at 10k resolution. Replay of committed reps showed rng stream offset (documented); null-distribution equivalence of committed vs appended verified (means/sds/quantiles match). |
| 1 | Bootstrap CIs for all headline metrics (median r, fold medians, M2 delta, scanner AUC) | DONE | results/p2_bootstrap_ci.json (seed 20260927, 10k patient-level boots): pathway r 0.326 [0.104, 0.571]; immune AUC 0.633 [0.575, 0.743]; m1 0.046 [0.036, 0.072]; m2 0.052 [0.039, 0.080]; M2-M1 delta 0.0061 [0.0046, 0.0116] (excludes 0); scanner AUC 0.926 [0.887, 0.954]. |
| 5 | Per-gene confidence intervals | DONE | results/p2_bootstrap_ci.json per_gene_r_m2: 266 genes, 2k boots over held-out sections; 171/266 with CI95_lo > 0. |
| 13 | Label-permutation scanner control (Track B AUC -> ~0.5) | DONE | results/p2_scanner_labelperm.json: exact replication of p2_model scanner (LogReg C=1.0), 20 label perms/section, 27 sections. Perm AUC median 0.506, p95 0.577 < 0.6 - signal is real, not artifact. True AUCs match fidelity rerun exactly; vs pre-rebuild CSV dev <=2.97e-3 (LogReg sensitivity to 1e-8 embedding deltas; scanner CSV refreshed to canonical rebuilt output). |
| 17 | Scanner feature-importance stability (bootstrap coefficients) | DONE | results/p2_scanner_stability.json: 10 LOPO LogReg coefficient vectors, pairwise Spearman median 0.857 [0.839, 0.874 fold-bootstrap 95%], min 0.753; top-50 features 62% Jaccard overlap, 100% sign consistency - reproducible morphological signature, not patient-specific noise. |
| 9 | Multi-task gene prediction: elastic-net multi-output, PLS, reduced-rank vs independent ridge | PARTIAL (positive) | results/p2_multitask.json: PLS(25) beats committed ridge in 10/10 folds (median 0.0673 vs 0.0462), RRR(25) in 7/10 (median 0.0596); ridge-refold sanity max abs diff 0.0103 vs committed. MultiTaskElasticNet leg pending: 3 configs never finished fold 1 (5000-iter >35 CPU-min, 1000-iter cyclic >15 CPU-min, 500-iter random 26+ CPU-min under ~4% VM duty); reordered to per-fold-checkpointed pass 2, running. |
| 11 | Nested evaluation of the 250-gene selection rule | DONE | results/p2_selection_sensitivity.json (pre-declared rules, ea70157): per-fold ridge medians - top-100 mean 0.0571, top-250 mean (committed reference) 0.0502, top-500 mean 0.0472, top-250 variance 0.0546 (mean/var panel overlap only 171.5/250). Headline median moves +-0.005 across rules; reference rule mid-pack; every rule's median stays below the locked G2-R 0.08, so the gate conclusion is selection-rule-robust. Coheres with #12: predictability not concentrated in top-expressed genes. |
| 12 | Negative-control gene experiment (gene-level) | DONE (expectation falsified, informative) | results/p2_gene_negcontrol.json: pre-declared design (5 draws x 250 genes from train-expression ranks 251-2250, seed 20260929, same ridge harness). Pre-declared ~0 expectation FALSIFIED: control median 0.0342 [0.0087, 0.0876] - per-gene predictability is broad across expressed genes (cell-type-composition-driven), not panel-specific. Selected-250 ridge beats controls 9/10 folds (0.0462 vs 0.0342, Wilcoxon p=0.0488); PLS beats controls 10/10 (0.0673, p=0.0020). Named anomaly BC23287 (controls 0.0814 > selected 0.0588) reported. Niche-level PancreasBeta control (D2) already existed. |
| 3 | Complete spatial null model for M2 (randomized graphs, shuffled coords, distance-preserving perms) | DONE | results/p2_spatial_null.json (pre-declared 3-null design, f8625b7; B=20, committed alpha + W recipe): real smoothing gain median +0.0073 (8/10 folds positive; BC23450 -0.0001, BC23901 -0.0041 named). Value-shuffle null -0.0503 (10/10 below real), graph-shuffle -0.0011 (9/10) - gain needs true neighborhoods. Distance-shell null retains +0.0039 (9/10 below real) - about half the gain is generic local denoising, half requires true geometry; both legs verified. |
| 4 | Simple spatial baselines: NN prediction, Gaussian smoothing, section-mean interpolation, Moran-based predictor | TODO | As #3. |

## Tier 2 - model/spec completeness (paper + sklearn-level compute)
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 8 | Full MLP specification + result table (layers, units, optimizer, epochs, params) | TODO | Paper-level; extract from committed M3 code. |
| 14 | Scanner subgroup analyses (AUC by subgroups per verdict) | DONE | results/p2_scanner_subgroups.json (29 sections): subtype-robust - HER2_luminal AUC 0.933 (n=14) vs HER2_non_luminal 0.926 (n=15); tumor-fraction tertiles 0.885-0.948; size tertiles ~0.93; hardest difficulty tertile (near-balanced labels) 0.847; one anti-predictive outlier (BC24105_D1, 0.065) sits in extreme-fraction stratum - named, not hidden. |
| 15 | Scanner calibration diagnostics | TODO | Reliability curve, Brier, ECE. |
| 18 | Full pipeline reproducibility artifacts | PARTIAL | Repo scripts/results committed. Gap: per verdict sub-items (env lock, run manifest). |
| 19 | Computational cost benchmarking | TODO | Rerun timing on committed pipelines; paper table. |

## Tier 3 - BLOCKED on sandbox-wipe re-acquisition
| # | Addition | Status | Evidence / Action |
|---|----------|--------|-------------------|
| 6 | Morphology ablations (stain/texture/color/entropy/nuclear-density only) | REBUILD-IN-PROGRESS | HER2ST re-acquired 2026-09-27 16:27 (273/273 files, sha256-verified vs Mendeley API); features_cache2 re-extraction launching. |
| 7 | Technical artifact controls (Macenko/Reinhard normalization) | REBUILD-IN-PROGRESS | Images re-acquired with HER2ST (68 JPGs, sha256-verified). |
| 10 | Learned image baseline, same patient-held-out split | REBUILD-IN-PROGRESS | CTransPath weights re-downloaded + loader verified (missing=head.* only, 0 unexpected); embeddings re-extraction launching. |
| 16 | External cohort validation, frozen model, no retraining | REBUILD-IN-PROGRESS | = PREREG-2 G2-X gate. External 10x data re-download owed (URLs live-verified pre-wipe); HER2ST side rebuilding now. |
| 20 | Challenge-style hidden evaluation protocol | TODO-CONSTRAINED | Extends #16; freeze code+model, evaluate on unseen public sections. |

Execution order: 2 -> 1+5 -> 13 -> 17 -> 3 -> 4 -> 9 -> 11 -> 12 -> 14,15,8,18,19 -> (rebuild) 6,7,10 -> 16 -> 20.
Every closure lands as: script + results JSON + tests where fitting + paper section increment (50+pp directive, user 12:02:56), committed with live-verified numbers.
