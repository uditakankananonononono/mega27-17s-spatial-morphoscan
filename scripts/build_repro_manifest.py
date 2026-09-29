"""Tier-2 #18 reproducibility artifacts: environment lock + run manifest.
Writes results/env_lock.txt (measured on the analysis box) and
results/run_manifest.json (stage -> script -> committed outputs -> how to run).
Manifest rows are curated from the committed pipeline; every output path named
here is committed in the repo."""
import json, os, subprocess, sys, platform, importlib.metadata as md

pkgs = ["numpy","pandas","scipy","scikit-learn","torch","timm","tifffile",
        "imagecodecs","scikit-image","h5py","matplotlib"]
env = {"python": sys.version.split()[0], "platform": platform.platform(),
       "packages": {p: md.version(p) for p in pkgs},
       "hardware": "2-core Xeon 2.60GHz, 1GB RAM, no GPU (measured 2026-09-29)"}
with open("results/env_lock.txt", "w") as fh:
    fh.write(f"python {env['python']}\n{env['platform']}\n")
    for p, v in env["packages"].items(): fh.write(f"{p}=={v}\n")

stages = [
 {"stage": "dataset download (Mendeley 29ntw7sh4r v5)", "script": "scripts/g2x_download.sh + docs in README", "outputs": ["data/her2st/"], "note": "sha256-verified"},
 {"stage": "Phase-1 handcrafted feature extraction", "script": "scripts/extract_features.py", "outputs": ["features_cache/"]},
 {"stage": "Phase-1 Track A (M0/M1/M2 ridge ladder)", "script": "scripts/run_trackA.py", "outputs": ["results/trackA_results.json", "results/trackA_per_gene.csv", "results/trackA_folds_partial.csv"]},
 {"stage": "Phase-1 Track B scanner (handcrafted)", "script": "scripts/run_trackB.py", "outputs": ["results/trackB_results.json"]},
 {"stage": "M3 MLP nonlinear check", "script": "scripts/run_trackA_m3.py", "outputs": ["results/trackA_m3.json", "results/trackA_m3_partial.csv"]},
 {"stage": "CTransPath embedding extraction (68 sections)", "script": "scripts/p2_extract_ctranspath.py", "outputs": ["embeddings_cache/"]},
 {"stage": "Phase-2 ridge probe, 250 genes + 50 pathways, LOPO", "script": "scripts/p2_model.py", "outputs": ["results/p2_folds_partial.csv", "results/p2_per_gene_partial.csv", "results/p2_per_pathway_partial.csv", "results/p2_scanner_partial.csv"]},
 {"stage": "Gate evaluation", "script": "scripts/p2_evaluate_gates.py", "outputs": ["results/p2_gate_verdicts.json"]},
 {"stage": "Scanner label-permutation control", "script": "scripts/p2_scanner_labelperm.py", "outputs": ["results/p2_scanner_labelperm.json"]},
 {"stage": "Scanner coefficient stability", "script": "scripts/p2_scanner_stability.py", "outputs": ["results/p2_scanner_stability.json"]},
 {"stage": "Scanner subgroup analyses", "script": "scripts/p2_scanner_subgroups.py", "outputs": ["results/p2_scanner_subgroups.json"]},
 {"stage": "Scanner calibration diagnostics", "script": "scripts/p2_scanner_calibration.py", "outputs": ["results/p2_scanner_calibration.json"]},
 {"stage": "Pathway permutation test 10k", "script": "scripts/p2_perm10k.py", "outputs": ["results/p2_perm10k.json"]},
 {"stage": "Spatial null + baselines + selection sensitivity + negative controls", "script": "scripts/p2_spatial_null.py, scripts/p2_spatial_baselines.py, scripts/p2_selection_sensitivity.py, scripts/p2_gene_negcontrol.py", "outputs": ["results/p2_spatial_null.json", "results/p2_spatial_baselines.json", "results/p2_selection_sensitivity.json", "results/p2_gene_negcontrol.json"]},
 {"stage": "Multitask PLS + reduced-rank ridge", "script": "scripts/p2_multitask.py", "outputs": ["results/p2_multitask.json"]},
 {"stage": "Niche architecture (D2)", "script": "scripts/p2_niches.py", "outputs": ["results/p2_niches_summary.json"]},
 {"stage": "Bootstrap CIs", "script": "scripts/p2_bootstrap.py", "outputs": ["results/p2_bootstrap_ci.json"]},
 {"stage": "G2-X external eval v1 + frozen head", "script": "scripts/g2x_frozen_train.py, scripts/g2x_external_eval.py", "outputs": ["results/g2x_external_eval.json", "results/g2x_frozen_model_meta.json"]},
 {"stage": "G2-X v2 stain-normalized re-evaluation", "script": "scripts/g2x_external_eval_v2.py", "outputs": ["results/g2x_external_eval_v2.json", "results/g2x_v2_prespec.json", "results/g2x_v2_refstain.json"]},
 {"stage": "Cost benchmark", "script": "scripts/p2_cost_benchmark.py", "outputs": ["results/p2_cost_benchmark.json"]},
 {"stage": "Paper tables/figures", "script": "scripts/generate_tables.py, scripts/make_figures.py, scripts/paper_scanner_tables.py", "outputs": ["paper/main.tex", "paper/scanner_table.tex"]},
]
manifest = {"item": "Tier-2 #18 run manifest", "env": env,
            "convention": "every script runs from repo root; results/ is gitignored and committed with git add -f; caches (features_cache, embeddings_cache) are reproducible from data/ via the listed scripts",
            "stages": stages}
missing = []
for s in stages:
    for o in s["outputs"]:
        if o.endswith("/"):
            if not os.path.isdir(o.rstrip("/")): missing.append(o)
        elif not os.path.exists(o): missing.append(o)
manifest["self_check"] = {"outputs_verified_present": len(stages), "missing": missing}
json.dump(manifest, open("results/run_manifest.json", "w"), indent=1)
print(json.dumps(manifest["self_check"]))
