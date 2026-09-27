## Weaknesses (20 computational-only)

1. Primary endpoint uses a median across genes that may obscure unstable fold behavior (Methods 3.4, Results 4.2).
   The headline M1 performance (median per-gene held-out Pearson r = 0.014) is dominated by the center of a distribution where most genes are near zero. A median can become statistically significant with many weakly shifted genes while having limited predictive utility. The paper does not report confidence intervals for the median, bootstrap uncertainty, or effect-size distributions.

2. The permutation test design does not match the model-selection pipeline (Methods 3.4, Equation 6).
   The gene-axis permutation preserves expression structure but destroys gene identity. However, gene selection, feature scaling, and ridge penalty selection occur inside folds. The permutation procedure does not clearly rerun the complete gene-selection pipeline if any preprocessing decisions depend on the observed gene targets. This may underestimate null variance.

3. Only 50 permutations are used for a p-value threshold near 0.05 (Methods 3.4).
   With 50 permutations, the minimum attainable two-sided resolution is coarse. A reported p = 0.020 is based on very few null samples and has substantial Monte Carlo uncertainty.

4. The benchmark comparison against ST-Net is not a fully matched computational reproduction (Introduction 1.1, Results 4.3).
   The paper compares patient-held-out MorphoScan results with published ST-Net numbers but does not reproduce ST-Net's preprocessing, patch extraction, architecture, training procedure, or exact evaluation code. The claim of “apples-to-apples” is therefore limited to cohort identity, not pipeline identity.

5. The MLP nonlinear check is underspecified (Methods 3.3 M3).
   The paper states that an MLP is used as a sanity check but provides no architecture, parameter count, training procedure, regularization, stopping criteria, or performance results. A negative or positive MLP result cannot be interpreted.

6. Spatial smoothing risks inflating correlation through spatial autocorrelation without a spatial null model (Methods 3.3 M2, Results 4.2).
   M2 improves correlation by averaging neighboring predictions. Because expression itself is spatially autocorrelated, improvement may reflect generic spatial smoothness rather than morphology-derived signal. Moran's I is measured but no spatial permutation control is reported.

7. The smoothing comparison lacks comparison against simpler spatial baselines (Methods 3.3 M2).
   The claimed “method advance” compares M2 only against M1. It does not compare against predicting neighborhood averages, Gaussian smoothing of the true training expression, nearest-neighbor interpolation, or random spatial kernels.

8. The top-250 gene selection procedure is insufficiently described (Methods 3.1).
   The paper states that the top-250 gene list is selected inside folds but does not specify whether ranking occurs by variance, expression abundance, prior ST-Net gene list overlap, or predictive criteria. Different choices could materially alter performance.

9. Gene-level multiple testing is unclear despite BH mention (Methods 3.4, Results 4.5).
   The paper discusses Benjamini–Hochberg but does not report adjusted p-values for the 250 gene correlations or enrichment analyses. It is unclear where multiple testing correction is applied.

10. The enrichment analysis is vulnerable to post hoc selection bias (Results 4.5, Appendix H).
   G-A5 evaluates the “top predictable genes” after observing held-out performance. Even though prediction evaluation is separated, pathway interpretation of the highest-ranked genes is not corrected for ranking bias.

11. The scanner evaluation aggregates heterogeneous patient distributions without subgroup analysis (Methods 3.5, Results 5).
   The AUC values combine sections with potentially extreme tumor fractions. The paper does not evaluate performance stratified by tumor fraction, subtype, section size, or scanner difficulty.

12. The scanner threshold analysis is insufficiently calibrated (Results 5).
   Sensitivity/specificity at threshold 0.5 is reported, but there is no patient-level threshold selection analysis, calibration curve, expected calibration error, or comparison with prevalence-based baselines.

13. Feature engineering may unintentionally encode technical artifacts (Methods 3.2).
   Color statistics, RGB moments, and stain features can capture scanner-specific variation. No computational stain normalization ablation or slide-level batch-effect analysis is provided.

14. Patient-held-out evaluation may still contain section-level leakage through preprocessing (Methods 3.1).
   The paper states patient-level separation but does not explicitly audit whether image normalization parameters, stain matrices, patch dimensions, or coordinate-derived statistics are computed only from training patients.

15. The transport experiment is not statistically quantified (Results 4.4).
   Transport medians are reported for Luminal A, Luminal B, and TNBC, but no confidence intervals, uncertainty estimates, or comparison against within-subtype null expectations are provided.

16. The dataset is small relative to model claims (Data 2.1, Results 4.2).
   Ten HER2-positive patients are used for the main benchmark. Even with patient-held-out evaluation, only ten independent test units exist, limiting confidence in generalization estimates.

17. The regression formulation ignores possible gene-gene structure (Methods 3.3 M1).
   Independent ridge regression predicts each gene separately. This discards biological covariance among genes and may underperform multi-task approaches.

18. The feature space is not compared against simpler dimensional baselines (Methods 3.2).
   A 41-feature morphology model is presented, but there is no comparison against RGB-only features, stain-only features, texture-only features, or random-feature controls.

19. The reproducibility claim relies heavily on unit tests rather than scientific rerun validation (Section 9).
   Thirty-two unit tests verify functions, but there is no reported full clean-environment rerun, container hash, dependency lock validation, or independent reproduction result.

20. The accession/tool ledger adds volume but little computational validation (Sections 6.1–6.2).
   The 51 tools and 629 records demonstrate documentation effort, but they do not directly validate the prediction pipeline. The manuscript spends substantial space on audit counts without measuring computational robustness.

## Additions (20 computational-only)

1. Add bootstrap confidence intervals for all headline metrics (Results 4.2).
   Report bootstrap intervals for median per-gene r, fold-level medians, M2 improvement, and scanner AUC.

2. Increase permutation testing resolution (Methods 3.4).
   Replace 50 permutations with at least 1,000–10,000 permutations or use an analytically justified Monte Carlo stopping rule.

3. Add a complete spatial null model for M2 (Methods 3.3).
   Test smoothing using randomized spatial graphs, shuffled coordinates, or distance-preserving permutation graphs to determine whether gains exceed generic autocorrelation.

4. Benchmark against simple spatial baselines (Methods M2).
   Add:
   
   nearest-neighbor prediction,
   
   Gaussian smoothing of training-fold expression,
   
   section-level mean interpolation,
   
   Moran-based spatial predictor.

5. Report per-gene confidence intervals (Results 4.2).
   Provide uncertainty around each gene correlation instead of only ranked point estimates.

6. Add morphology ablation experiments (Methods 3.2).
   Train separate models using:
   
   stain features only,
   
   texture only,
   
   color only,
   
   entropy only,
   
   nuclear-density proxy only.

7. Add technical artifact controls (Methods 3.2).
   Include slide normalization experiments:
   
   Macenko normalization,
   
   Reinhard normalization,
   
   raw-image comparison.

8. Provide a full MLP specification and result table (Methods 3.3 M3).
   Report layers, hidden units, optimizer, epochs, regularization, parameter count, and final performance.

9. Add multi-task gene prediction models (Methods 3.3).
   Compare independent ridge against:
   
   elastic-net multi-output regression,
   
   partial least squares,
   
   reduced-rank regression.

10. Add a learned image baseline trained under the same split (Results 4.3).
   Implement a lightweight CNN or pretrained feature extractor using identical patient-held-out folds to contextualize the handcrafted model.

11. Perform nested evaluation of the 250-gene selection rule (Methods 3.1).
   Report whether results change when selecting genes by:

12. Add a negative-control gene experiment (Results 4.5).
   Evaluate genes expected not to be morphology-linked and show their correlation distribution.

13. Add a label-permutation scanner control (Track B).
   Randomly permute tumor labels within patients and confirm AUC collapses toward 0.5.

14. Add scanner subgroup analyses (Results 5).
   Report AUC by:

15. Add calibration diagnostics for the scanner (Results 5).
   Include:

16. Add external cohort computational validation if available (Transport 4.4).
   Evaluate the frozen model on another public breast spatial dataset without retraining.

17. Add feature importance stability analysis (Results 5, Appendix D).
   Run bootstrap resampling to determine whether scanner coefficients remain stable.

18. Add full pipeline reproducibility artifacts (Section 9).
   Include:

19. Add computational cost benchmarking (Discussion).
   Report:

20. Add a prospective “challenge-style” hidden evaluation protocol (Conclusion).
   Freeze the code and model, then evaluate on an unseen public section set or held-out accession to demonstrate that the reported patient-held-out performance is not specific to the current cohort.
