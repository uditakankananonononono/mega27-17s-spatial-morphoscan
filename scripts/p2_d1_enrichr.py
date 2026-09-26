"""D1 enrichment (PREREG-2 round-2 adoption): split confirmatory / exploratory.
Input set: held-out genes with median r_m1 >= 0.15 (the locked D1 count set).
Confirmatory (single pre-registered test, carried from Phase-1): antigen
processing & presentation enriched in GO_Biological_Process_2023.
Exploratory: full Enrichr libraries, BH-adjusted, reported with adjusted p.
Also runs the r_m2>=0.15 set as a sensitivity analysis."""
import os, sys, json, time
import numpy as np, pandas as pd
import urllib.request, urllib.parse

R = "results"
pg = list(__import__("csv").DictReader(open(os.path.join(R, "p2_per_gene_partial.csv"))))
def f(x):
    try: return float(x)
    except: return np.nan
gene_r1, gene_r2 = {}, {}
for r in pg:
    gene_r1.setdefault(r["gene"], []).append(f(r["r_m1"]))
    gene_r2.setdefault(r["gene"], []).append(f(r["r_m2"]))
med1 = {g: float(np.nanmedian(v)) for g, v in gene_r1.items()}
med2 = {g: float(np.nanmedian(v)) for g, v in gene_r2.items()}
sym = json.load(open("data/ensg_symbol_map.json"))
set1 = sorted([sym.get(g, g) for g, v in med1.items() if v >= 0.15])
set2 = sorted([sym.get(g, g) for g, v in med2.items() if v >= 0.15])
set1 = [g for g in set1 if not g.startswith("ENSG")]
set2 = [g for g in set2 if not g.startswith("ENSG")]
print(f"confirmatory set (r_m1>=0.15): {len(set1)} genes")
print(f"sensitivity set (r_m2>=0.15): {len(set2)} genes")

BASE = "https://maayanlab.cloud/Enrichr"
def enrich(gene_list, library):
    boundary = "----mega27"
    body = "\r\n".join([
        f"--{boundary}",
        'Content-Disposition: form-data; name="list"', "", "\n".join(gene_list),
        f"--{boundary}",
        'Content-Disposition: form-data; name="description"', "", "mega27 D1",
        f"--{boundary}--", ""]).encode()
    req = urllib.request.Request(BASE + "/addList", data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}",
                 "User-Agent": "mega27-research/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        uid = json.loads(r.read())["userListId"]
    time.sleep(0.5)
    url = f"{BASE}/enrich?userListId={uid}&backgroundType={urllib.parse.quote(library)}"
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read())[library]

LIBS = ["GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020", "KEGG_2021_Human"]
ANTIGEN = "Antigen Processing And Presentation"

out = {"set_r_m1": set1, "set_r_m2_sensitivity": set2, "results": {}}
for tag, genes in [("primary", set1), ("sensitivity_r_m2", set2)]:
    if len(genes) < 5:
        out["results"][tag] = {"error": "too few genes"}; continue
    res = {}
    for lib in LIBS:
        try:
            rows = enrich(genes, lib)
        except Exception as e:
            res[lib] = {"error": str(e)}; continue
        res[lib] = [dict(term=r[1], p=r[2], adj_p=r[6], overlap=f"{r[3] and len(r[5]) or 0}/{r[4] if len(r)>4 else '?'}",
                         genes=r[5]) for r in rows[:15]]
        time.sleep(0.5)
    out["results"][tag] = res

# confirmatory verdict: antigen-processing term in GO BP, primary set
conf = {"hypothesis": "antigen processing & presentation enrichment (GO_Biological_Process_2023)",
        "alpha": 0.05, "one_sided": True}
go = out["results"].get("primary", {}).get("GO_Biological_Process_2023", [])
hit = [t for t in go if isinstance(t, dict) and ANTIGEN.lower() in t["term"].lower()]
conf["hit_terms"] = hit
conf["pass"] = bool(hit and hit[0]["p"] < 0.05)
out["confirmatory"] = conf
json.dump(out, open(os.path.join(R, "p2_d1_enrichr.json"), "w"), indent=1)
print("CONFIRMATORY:", json.dumps({k: v for k, v in conf.items() if k != "hit_terms"}))
for t in hit[:3]:
    print(" ", t["term"], "p=", t["p"], "adj_p=", t["adj_p"], t["genes"])
print("TOP EXPLORATORY (GO BP, primary):")
for t in go[:8]:
    if isinstance(t, dict): print(" ", t["term"], "adj_p=", t["adj_p"])
print("DONE")
