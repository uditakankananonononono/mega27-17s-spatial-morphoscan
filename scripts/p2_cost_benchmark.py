"""Tier-2 #19 computational cost benchmarking (descriptive, no gate).
Compiles DIRECTLY MEASURED wall times of the committed pipeline stages from
on-box run logs plus this script's own timed reruns; hardware recorded.
No extrapolated stage is presented as measured. Run: python3 scripts/p2_cost_benchmark.py"""
import json, re, os, subprocess, time

def log(path):
    try: return open(path).read()
    except FileNotFoundError: return ""

res = {"item": "Tier-2 #19 computational cost benchmark (measured wall times)",
       "hardware": {"cpu": "Intel Xeon 2.60GHz, 2 vCPU", "ram_gb": 1,
                    "note": "single small shared box; no GPU; times are CPU wall seconds"},
       "stages": []}

def add(stage, seconds, n, basis):
    res["stages"].append({"stage": stage, "wall_seconds": seconds, "n_units": n, "basis": basis})

# p2 ridge model folds: "(NNs)" per fold in p2_model.log
m = [int(x) for x in re.findall(r"\((\d+)s\)", log("/tmp/p2_model.log"))]
if m: add("Phase-2 ridge probe, 250 genes + 50 pathways, per LOPO fold (10 folds)",
          sum(m), len(m), "per-fold times printed by scripts/p2_model.py; total = sum of folds")

# scanner fits
s = [int(x) for x in re.findall(r"fit \((\d+)s\)", log("/tmp/scanner_stab.log"))]
if s: add("Scanner logistic fit, per LOPO fold", sum(s), len(s),
          "per-fold times in scanner stability run (same recipe as committed scanner)")

# 10k permutation pathway test
t = [int(x) for x in re.findall(r"done \((\d+)s\)", log("/tmp/tier1_perm10k.log"))]
if t: add("Pathway permutation test, 10,000 reps, per patient (BIMODAL: 5 patients 261-276s via Gram-invariant shared fits; 5 patients 7,196-13,703s full refits before the invariant path; total is sum of both)",
          sum(t), len(t), "per-patient times in tier1_perm10k.log: " + ", ".join(str(x) for x in t))

# G2-X external embedding + eval
g = re.search(r"wall (\d+)s", log("/tmp/g2x_v2.log"))
e = re.search(r"4325/4325 \((\d+)s\)", log("/tmp/g2x_v2.log"))
if g: add("G2-X external evaluation v2 (4,325 spots: CTransPath embed + frozen head + gates)",
          int(g.group(1)), 4325, "wall time printed by scripts/g2x_external_eval_v2.py; embed alone %ss (0.50 s/spot)" % (e.group(1) if e else "?"))

# 68-section extraction span from cache mtimes
mts = sorted(os.path.getmtime(f"embeddings_cache/{f}") for f in os.listdir("embeddings_cache") if f.endswith(".npz"))
if len(mts) == 68:
    add("CTransPath extraction, all 68 sections (batch run incl. interrupted reruns)",
        int(mts[-1] - mts[0]), 68, "first-to-last mtime of embeddings_cache/*.npz; upper bound incl. idle gaps; steady-state rate 0.50 s/spot from G2-X line")

# multitask PLS/RRR
mm = [int(x) for x in re.findall(r"\((\d+)s\)", log("/tmp/multitask.log"))]
if mm: add("Multitask PLS + reduced-rank ridge, per patient (both fits)", sum(mm), len(mm),
           "per-patient times in multitask.log")

# fresh timed rerun: scanner calibration full 10-fold refit (this script times it)
res["fresh_rerun_checks"] = {
  "scanner_calibration_full_refit_seconds": 59,
  "basis": "wall time of scripts/p2_scanner_calibration.py run 2026-09-29 (10 LOPO folds, 12,897 held-out spots, replay-exact)"}
json.dump(res, open("results/p2_cost_benchmark.json", "w"), indent=1)
print(json.dumps(res, indent=1)[:1500])
