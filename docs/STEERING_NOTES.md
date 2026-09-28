# STEERING_NOTES.md - user ground rules 2026-09-28 1:17 PM IST + literature consultation

User rules (verbatim via main agent): (1) negatives stay a MINORITY per project - achieved by
steering to promising directions, never by softening measured numbers; (2) when a direction
stalls, ask ChatGPT for the next steering direction and keep learning from newly published
literature - document what was consulted and taken; (3) creativity and novelty lead;
discovery from existing data counts.

## Literature consultation 2026-09-28 (web search, 2026-06+ publications)
Consulted recent histology-to-spatial-expression prediction literature. What we take:

1. **DeepSpot-M** (medRxiv 2026-06): multimodal foundation model for transcriptome-wide
   "virtual spatial transcriptomics" from histology - the current frontier baseline class
   for our exact task. Steering: position our CTransPath-ridge as the interpretable,
   lightweight alternative; if DeepSpot-M weights/code are public, a head-to-head on
   HER2ST is a novelty-forward benchmark (tool-beats-benchmark element).
2. **HistoGPA** (arXiv 2607.24364): gene-prior attention framework for histology-based
   spatial gene prediction - explicit gene-prior conditioning contrasts with our
   train-selected top-250; a gene-prior augmentation is a candidate addition.
3. **VOICE** (alphaXiv 2608.08366): vision-omics foundation model, direct + retrieval
   prediction of in-situ single-cell expression - retrieval-based prediction is an
   alternative null/comparison.
4. **R4ST** (reference-guided graph-generative reconstruction) and foundation-model
   contrastive histology-to-transcriptomics translation (Comput. Biol. 2026-09):
   generative/contrastive direction gaining ground.

Steering decision (rule 1/3): our VERIFIED positives lead the paper - pathway-level
prediction (median r 0.326, 95% CI [0.104, 0.571], permutation p resolution now 10k),
scanner tumor-vs-other AUC 0.926 [0.887, 0.954] with label-permutation control PASSED,
spatial smoothing gain (M2-M1 delta CI excludes 0). Novelty path: benchmark against a
2026 foundation baseline (DeepSpot-M preferred if artifacts public) rather than only
HGB; pathway-level prediction as our distinct axis (gene-level r is modest everywhere
in the literature - pathway programs are where morphology carries signal). Negatives
stay compact in limitations.

## 2026-09-28 follow-up: DeepSpot-M artifacts CONFIRMED public (live check)

- Code: https://github.com/ratschlab/DeepSpotM (PolyForm Noncommercial 1.0.0 - non-commercial OK)
- Weights: https://huggingface.co/ratschlab/DeepSpotM (CC-BY-NC-SA-4.0 - non-commercial + attribution OK)
- `DeepSpotM.from_pretrained("ratschlab/DeepSpotM", source="scgpt")` predicts a ~19k-gene panel
  from 224x224 H&E tiles; `predict_genes` allows a targeted gene subset (cheaper).
- Feasibility for HER2ST head-to-head: our spots are 100um diameter (approx 224px at 20x);
  image tiles are already extracted for CTransPath - same tiles can feed DeepSpot-M. Then
  compare per-spot pathway-level r on the same 10 pathways and gene-level r on the shared
  gene panel. Compute: single forward pass per tile per fold is cheap; full WSI example needs
  pyvips (installable). License: our use is non-commercial academic benchmarking - within terms.
- NEXT ACTION (after #9 multitask): download weights (safetensors ~size TBD), run tile-level
  inference on held-out fold patients, pathway-score both prediction sets, report r comparison.
