# Literature calibration: raw cross-site transfer of histology->expression models
Collected 2026-09-29 after G2-X FAIL (34afa2f: frozen internal model, external 10x Visium breast, median pathway r = -0.1547 vs gate >= 0.10). Purpose: establish whether raw-transfer failure is field-typical before designing v2.

## DeepSpot-M (medRxiv 2026.06.19.26356060, ratschlab; PyPI v1.0.0, HF ratschlab/DeepSpotM)
- Trained on 730,000 paired histology-transcriptomic profiles, 14 cancer types, multiple institutions.
- ZERO-SHOT on external Xenium samples, no single-cell training: mean Pearson 0.19 across 15 genes. Finetuning on ONE external sample doubles it to 0.38 (test-time adaptation needed for usable transfer).
- Stromal genes transfer best zero-shot: COL1A2 0.47, LUM 0.41 (their Fig 3C,D).
- Source: https://www.medrxiv.org/content/10.64898/2026.06.19.26356060v1.full-text

## Concordance with our G2-X FAIL
- Our only positive-transferring programs on the external section were stromal: Wnt-beta-catenin 0.235, Angiogenesis 0.222 (EMT 0.095). DeepSpot-M shows the identical stromal-transfer bias (COL1A2/LUM). Independent replication of the same direction, across different models and cohorts.
- Even a pan-cancer foundation model gets only ~0.19 mean Pearson raw zero-shot on external assays and needs per-slide adaptation for 0.38. Our internal model (29 sections, one study) dropping below zero raw on a new site's section is consistent with the field, not an anomaly.

## Pathology foundation models under site/stain shift (general)
- "Early convolutional networks achieved strong in-domain performance but failed across centres ... PFMs remain fragile under realistic perturbations and unseen domains" - The Good, the Bad, and the Brittle: Benchmarking Robustness and Generalisation of Histopathology Foundation Models (Yajnik & Minhas, Warwick), https://arxiv.org/html/2607.04401 (2026-07-05).
- Stain/scanner variation measurably degrades pathology foundation models: Impact of tissue staining and scanner variation on the performance of pathology foundation models, https://www.biorxiv.org/content/10.1101/2025.08.18.670932v2 (also PMC12932120).

## Implication for v2 design (feeds G2-X v2 prespec + DeepSpot-M head-to-head)
1. Raw-transfer FAIL is field-typical; the paper's external-validity contribution is quantifying the gap + the stromal concordance, not passing the raw gate.
2. DeepSpot-M head-to-head on the SAME external section: if a 730k-profile foundation model also degrades raw on this exact tissue, the calibration is direct, not by analogy. Their published zero-shot mean (0.19, Xenium genes) gives the reference band.
3. Stain normalization (Macenko/Vahadane) is the literature-standard first lever; v2 prespec must keep the frozen model + same gate, normalize tiles only, with full lineage disclosure.
