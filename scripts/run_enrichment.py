"""G-A5: g:Profiler enrichment of top-25 predictable genes (by median held-out r)."""
import json, urllib.request
import pandas as pd

G = pd.read_csv("results/trackA_per_gene.csv")
top = (G.groupby("gene")["r_m1"].median().sort_values(ascending=False).head(25))
genes = [g for g in top.index if g.startswith("ENSG")]
body = json.dumps({"organism": "hsapiens", "query": genes,
                   "sources": ["GO:BP", "KEGG", "REAC"]}).encode()
req = urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
                             data=body, headers={"Content-Type": "application/json"})
res = json.loads(urllib.request.urlopen(req, timeout=60).read())
terms = [{"source": t["source"], "term": t["name"], "p": t["p_value"],
          "intersection": len(t["intersections"])} for t in res["result"][:15]]
out = {"genes": genes, "n_genes": len(genes), "terms": terms}
json.dump(out, open("results/enrichment_top25.json", "w"), indent=1)
print(json.dumps(out, indent=1)[:1500])
