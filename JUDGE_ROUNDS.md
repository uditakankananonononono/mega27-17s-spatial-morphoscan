# JUDGE_ROUNDS - MorphoScan (MEGA27 item 17-spatial)

Standing rules (user, WhatsApp, verbatim in channel history):
- 2026-09-26 4:11:18: "...ask CHATGPT about more ideas... complete each project... never count a negative as a result of a research project... minimum 10 judging rounds about weaknesses in the project and what more to include to make it better."
- 2026-09-26 4:11:49: "make sure each project comes out beating benchmarks (improve until it does) and actually takes out a new discovery in each case"
- 2026-09-26 4:12:25: "ask chatgpt how to redirected when you see a negative and failure is not moving forward"
- 2026-09-26 4:14:37: "take inspiration from previosu isef winners like natasha's project"

Every round logs the verbatim prompt and verbatim response. ChatGPT output is
external advice, not empirical evidence: every resource it names is
independently verified before use, and all new gates are locked by us before
outcomes are seen. Locked PREREG.md gates are unchangeable; redirections
become NEW locked phases.

---

## Round 1 - 2026-09-26 16:29-16:31 IST - REDIRECTION on the 4 failures (rule 6)

Surface: ChatGPT (user's own account, Free tier), cloud browser.
Conversation URL: https://chatgpt.com/c/6ab7a5c8-8f20-83ee-9eb7-52b4c77b87cf

### Verbatim prompt

You are a senior science-fair judge and computational-biology mentor. A pre-registered study (gates locked before outcomes; locked gates cannot change, so fixes become NEW locked phases, not re-fished old ones) measured these honest results. Give ranked redirection options per failure - concrete, free/open-data, executable by one student. PROJECT: MorphoScan - predict spot-level spatial gene expression from H&E morphology in breast cancer (HER2-ST cohort: 68 sections, 23 patients; leave-one-patient-out on 10 HER2+ patients, 29 sections, 250 genes; 41 handcrafted morphology features + ridge; M2 adds Gaussian spatial-graph smoothing). MEASURED: G-A1 PASS ridge beats train-mean, median per-gene Pearson 0.014 vs 0.000, p=0.020. G-A2 FAIL M2-M1 = +0.0029 needed +0.005 (Wilcoxon p=3e-7, right direction, too small). G-A3 FAIL 0.014 vs published ST-Net ~0.19 and HisToGene ~0.21 on the same cohort. G-A4 transport LumA 0.036, LumB 0.015, TNBC 0.008. G-A5 PASS top-25 predictable genes enrich antigen processing p=9.2e-4. G-B1 FAIL tumor-region scanner patient-held-out AUC 0.642-0.676 needed 0.80. G-B2 PASS Brier 0.189. FAILURES TO REDIRECT: (1) A2 method advance too small; (2) A3 benchmark gap vs deep nets; (3) B1 scanner AUC; (4) weak cross-subtype transport. For EACH give 2-3 ranked pivots with: the exact free resource (e.g. open H&E foundation-model embeddings like UNI or CONCH, public spatial-transcriptomics datasets), the new lockable gate with a numeric threshold, and why a judge would accept it. End with the single strongest overall redirection.

### Verbatim response

Senior ISEF computational-biology judge assessment

The pre-registration discipline matters here: A2 and A3 are not “tuning failures”; they are failed hypotheses. The correct scientific move is not to keep adding features until the original gates pass. The strongest redirections should create a new locked phase with a different biological question or a stronger measurement regime.

Overall diagnosis:

Morphology → expression is real but weak (A1, A5 pass).

Naive spatial smoothing is not enough (A2).

Handcrafted morphology is the bottleneck (A3).

Scanner/tumor classification is underpowered because morphology alone is noisy and patient-shifted (B1).

Subtype transport is biologically plausible but unstable (A4).

Below are ranked pivots.

Failure 1: A2 — Spatial graph smoothing gives only +0.0029 Pearson
Original hypothesis failed

“Adding spatial neighborhood information substantially improves spot-level expression prediction.”

The result says:

Direction correct.

Effect statistically significant.

But biologically/practically tiny.

A judge will not accept “try more smoothing parameters.” That is hyperparameter fishing.

Pivot A2-1 (strongest): Replace smoothing with learned spatial representation
New hypothesis

Spatial context is useful, but fixed Gaussian smoothing cannot learn tumor architecture.

Free resources

Public spatial transcriptomics:

HER2-ST cohort already available.

10x Visium breast cancer datasets from:

Spatial Transcriptomics repository

10x Genomics public datasets

Human Tumor Atlas Network (HTAN)

Model:

PyTorch Geometric

Graph Attention Network (GAT)

GraphSAGE

New locked gate

Before touching test patients:

LOPO patient-held-out median gene Pearson:

Baseline ridge: 0.014

New locked target:

GAT must achieve median Pearson ≥0.035 and improve ≥70% of evaluated genes over ridge.

Secondary:

permutation spatial-edge control must remove the gain:

real graph > shuffled graph by ≥0.015 Pearson.

Why judges accept

Because this changes the scientific question:

Not:

“Can smoothing help?”

But:

“Can learned tissue topology extract spatial biological programs?”

That is a genuine method advance.

Pivot A2-2: Predict biological programs instead of individual genes
Hypothesis

Individual genes are noisy; morphology predicts pathway-level states.

Free resources

MSigDB Hallmark gene sets

KEGG pathways

GO biological processes

Convert 250 genes → pathway scores.

New locked gate

Predict 20–50 pathway activities:

Accept if:

median pathway Pearson ≥0.20

permutation p <0.01

immune pathway AUROC ≥0.75

Why judges accept

Spatial transcriptomics is noisy. Showing morphology recovers immune/stromal/tumor programs is more biologically meaningful than chasing every gene.

Pivot A2-3: Morphology uncertainty map
Hypothesis

The failure is not prediction failure; it is identifying where morphology contains signal.

Gate

Model must identify:

top 10% highest-confidence spots contain ≥2× enrichment of predictable genes versus random spots.

Why accepted

Turns a weak predictor into a biological discovery tool.

Failure 2: A3 — Huge gap vs ST-Net/HisToGene
Original failure

Handcrafted 41 morphology features:

Pearson:

MorphoScan: 0.014

ST-Net: ~0.19

HisToGene: ~0.21

The judge interpretation:

The representation is the bottleneck, not the regression model.

Pivot A3-1 (strongest): Replace handcrafted features with pathology foundation embeddings
Hypothesis

Modern histology foundation models extract morphology unavailable to handcrafted statistics.

Free resources

Open models:

UNI (Virchow/Pathology Foundation Model family)

CONCH vision-language pathology model

CTransPath

Available through:

Hugging Face

PyTorch

Pipeline:

H&E tile → embedding vector → ridge/MLP/GNN → expression.

New locked gate

Patient-held-out LOPO:

Must achieve:

median gene Pearson ≥0.08

AND

beat handcrafted ridge in ≥80% of genes.

Why judges accept

This is the correct modernization path.

You are no longer competing with 2018-era morphology features against 2024 deep models.

Pivot A3-2: Narrow to predictable gene discovery
Hypothesis

Not all genes are visually encoded.

Instead of predicting all 250:

Find the visually predictable subset.

Gate

Discovery phase:

Select genes only using training patients.

Locked test:

At least:

50 genes with Pearson ≥0.15

FDR <0.05

AND pathway enrichment.

Why judges accept

A failed universal predictor can become a biomarker discovery system.

Pivot A3-3: Cross-dataset validation
Hypothesis

Current score is limited by cohort size.

Free resources

Combine:

HER2-ST

breast cancer Visium datasets

TCGA H&E + RNA-seq

Gate

External dataset:

same frozen model

Pearson improvement over mean baseline

no retraining

Threshold:

≥0.05 median Pearson.

Why accepted

External validation often matters more than squeezing the original cohort.

Failure 3: B1 — Tumor-region scanner AUC only 0.64–0.68
Interpretation

The classifier learned weak morphology correlates.

A judge will ask:

“Why does this matter if humans already know tumor regions?”

Need either stronger representation or a different biological endpoint.

Pivot B1-1 (strongest): Foundation-model tumor segmentation
Hypothesis

The failure is feature quality.

Free resources

Use:

UNI embeddings

CONCH embeddings

Pathology segmentation datasets:

PANDA

TCGA breast WSIs

Classifier:

embedding → logistic regression / small MLP.

New gate

Patient-held-out:

AUC ≥0.85

AND:

Dice score ≥0.75 on public annotated regions.

Why accepted

A pathology AI judge recognizes this as a legitimate computational pathology benchmark.

Pivot B1-2: Change task from tumor detection to biological niche classification
Hypothesis

Tumor/non-tumor is too crude.

Predict:

immune-rich

stromal

necrotic

proliferative regions

Resources

TCGA pathology annotations

TIL datasets

Gate

Macro AUROC:

≥0.80 across ≥3 tissue compartments.

Why accepted

More biologically useful than binary tumor scanning.

Pivot B1-3: Use weak supervision
Hypothesis

Spot labels are insufficient.

Resource

CLAM multiple-instance learning

TCGA WSIs

Gate

Slide-level cancer subtype prediction:

AUC ≥0.85.

Why accepted

Weak supervision is standard in computational pathology.

Failure 4: Weak subtype transport

Results:

Subtype	Pearson
LumA	0.036
LumB	0.015
TNBC	0.008

Interpretation:

Morphology-expression relationship is subtype-dependent.

Pivot A4-1 (strongest): Domain adaptation across subtypes
Hypothesis

The model learned HER2-specific morphology.

Free resources

Methods:

CORAL domain adaptation

adversarial domain adaptation

scikit-learn/PyTorch

Datasets:

TCGA BRCA

CPTAC breast cancer

Gate

Train:

HER2 + LumA

Test:

TNBC

Accept:

TNBC median Pearson ≥0.03

AND

not worse than HER2 baseline by >20%.

Why accepted

This directly addresses biological generalization.

Pivot A4-2: Build subtype-specific models
Hypothesis

One universal predictor is unrealistic.

Gate

For each subtype:

minimum:

Pearson ≥0.05

≥25 predictable genes

Why accepted

Precision medicine often requires subtype-specific models.

Pivot A4-3: Predict conserved pathways instead of genes
Gate

Across all four subtypes:

pathway Pearson ≥0.20

immune/proliferation pathways replicated in ≥3 subtypes.

Why accepted

Biological programs are more conserved than individual transcripts.

Ranking all pivots by ISEF strength
Rank	Pivot	Feasibility	Novelty	Judge appeal
1	H&E foundation embeddings + pathway prediction	★★★★★	★★★★	★★★★★
2	Learned spatial graph model	★★★★	★★★★★	★★★★★
3	Predict morphology-linked gene programs	★★★★★	★★★★	★★★★★
4	Subtype domain adaptation	★★★	★★★★★	★★★★
5	Foundation segmentation for tumor regions	★★★★	★★★	★★★★
6	Subtype-specific predictors	★★★★★	★★★	★★★
Single strongest overall redirection
MorphoScan 2.0: “Discovering morphology-encoded spatial biological programs using pathology foundation embeddings and patient-held-out validation.”

New locked hypothesis:

Frozen H&E foundation representations contain spatially localized information that predicts biologically meaningful transcriptomic programs better than handcrafted morphology.

Locked gates:

Representation gate

UNI/CONCH embeddings outperform handcrafted features:

median gene Pearson ≥0.08

Biology gate

pathway prediction:

median pathway Pearson ≥0.20

Generalization gate

external breast cancer dataset:

median pathway Pearson ≥0.10

Ablation gate

shuffled tissue coordinates reduce performance by ≥50%.

Why this survives judging:

It does not pretend the original method succeeded.

It converts a weak gene predictor into a mechanistically interpretable spatial biology discovery platform.

It uses free public data.

It is executable by one student.

Every claim can be locked before seeing outcomes.

The original MorphoScan idea is not dead; the failed gates reveal that the scientific unit should move upward from “predict every gene” to “identify which biological programs are visually encoded in tissue architecture.” That is a much stronger ISEF-level question.

### Independent assessment (ours, not ChatGPT's)

Adopt as the Phase-2 spine, subject to verification before any gate is locked:
1. A3-1 + A2-2 fused = the strongest overall redirection (ChatGPT's own #1 and
   #2/#3 ranks agree): frozen pathology foundation embeddings (UNI or CONCH or
   CTransPath - verify license + CPU-feasible inference) + pathway-level
   prediction target. New locked phase, not a re-fish of G-A2/G-A3.
2. B1-1 (foundation-embedding tumor segmentation, patient-held-out AUC >= 0.85)
   as the scanner redirection; verify a free annotated-region source for Dice.
3. A4 pivots deferred until Phase-2 representation exists; A4-3 (conserved
   pathways across subtypes) folds naturally into the pathway-level target.
4. A2-1 (GAT on learned spatial graphs) kept as Phase-3 candidate only if the
   Phase-2 representation gate passes; graph attention on 29 sections is
   overfit-prone at this cohort size.

Verification queue BEFORE locking PREREG-2 (no gate locks on ChatGPT's say-so):
- UNI / CONCH / CTransPath: weights availability, license, embedding dim,
  CPU inference cost for ~100k spots-tiles.
- ST-Net published 0.19: exact gene set, split, and normalization so any
  Phase-2 comparison is same-cohort same-split fair.
- Public breast Visium / HTAN cohorts: accession, size, subtype labels for
  the external-validation gate.
- PANDA / TCGA-BRCA annotations usable for scanner Dice gate.

Original FAIL gates G-A2, G-A3, G-B1 remain documented negatives in PREREG.md;
they are not edited, hidden, or re-fished. Phase 2 is new locked work.
