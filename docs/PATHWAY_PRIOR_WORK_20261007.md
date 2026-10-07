# Pathway prior-work check, 2026-10-07

## Scope

This is a focused source check, not a benchmark run or a new performance result.
The existing measured pathway and scanner results are unchanged. DeepSpot-M remains
frozen. No weights were fetched and no model was trained for this note.

## Verified prior work

**PEaRL (WACV 2026).** The proceedings paper describes an ssGSEA pathway-score
encoder aligned with histology through contrastive learning. It evaluates both gene
and pathway expression prediction on breast, skin, and lymph-node spatial data.
The full proceedings text was fetched for this check. Pathway prediction from
histology is therefore established prior work, not itself a unique MorphoScan claim.

Source: https://openaccess.thecvf.com/content/WACV2026/papers/Majumder_PEaRL_Pathway-Enhanced_Representation_Learning_for_Gene_and_Pathway_Expression_Prediction_WACV_2026_paper.pdf

**STP-BENCH (arXiv preprint, submitted September 5, 2026).** The fetched abstract
reports a unified pathology encoder across predictive methods where applicable,
gene-set recoverability analysis, and model-ranking changes after encoder control.
This supports separating representation effects from predictor effects when
planning a future comparison. This check read the abstract, not a reproduction of
the benchmark. Its findings retain preprint status.

Source: https://arxiv.org/abs/2609.05956

## Consequence for this project's framing

Keep the current conservative claims: patient-held-out lightweight probes,
explicit spatial controls, tested tools, and a frozen transport boundary.
Do not describe pathway-level prediction alone as first or unique. The current
paper's claims paragraph already makes no first-in-field claim; no headline
result needs retraction based on these sources.

Published correlations are not a head-to-head with this project's results.
Pathway definitions, score construction, dataset, splits, morphology encoders,
training procedures, and aggregation units must match before comparing numbers.
PEaRL's ssGSEA targets differ from this project's mean-z-score Hallmark targets.

## Open gaps

- No matched PEaRL or STP-BENCH comparison has been run here.
- STP-BENCH methods and implementation were not inspected in this focused check.
- These two sources are a bounded prior-work update, not an exhaustive review.
- No new positive result or completion claim follows from this note.
- Any future experiment needs a distinct, predeclared question and protocol;
  the failed external-transfer surface and frozen DeepSpot-M work remain unchanged.
