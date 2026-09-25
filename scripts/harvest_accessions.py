"""Harvest accession-level public dataset records for the MorphoScan dataset ledger.
Every row: accession, source, title, url, level_of_use.
level_of_use: analyzed (data pulled into the study), metadata-mined (registry
metadata harvested), primary (the core Mendeley cohort sections)."""
import json, time, urllib.parse, urllib.request, csv, sys

UA = {"User-Agent": "morphoscan-research/0.1 (academic)"}

def get(url, timeout=40):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

rows = []
def add(acc, source, title, url, use):
    rows.append(dict(accession=acc, source=source, title=(title or "")[:300], url=url, level_of_use=use))

# 1) GEO DataSets: spatial transcriptomics series + samples
try:
    term = urllib.parse.quote('"spatial transcriptomics"[All Fields]')
    r = json.loads(get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={term}&retmax=400&retmode=json"))
    ids = r["esearchresult"]["idlist"]
    print("GEO ids:", len(ids))
    for i in range(0, len(ids), 150):
        chunk = ",".join(ids[i:i+150])
        s = json.loads(get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={chunk}&retmode=json"))
        for uid in s["result"]["uids"]:
            rec = s["result"][uid]
            acc = rec.get("accession", "")
            title = rec.get("title", "")
            kind = "GSE" if acc.startswith("GSE") else ("GSM" if acc.startswith("GSM") else "GDS")
            add(acc, f"GEO-{kind}", title, f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}", "metadata-mined")
        time.sleep(0.4)
except Exception as e:
    print("GEO fail:", e)

# 2) Zenodo: spatial transcriptomics records
try:
    z = json.loads(get("https://zenodo.org/api/records?q=" + urllib.parse.quote("spatial transcriptomics") + "&size=60&sort=mostrecent"))
    for h in z.get("hits", {}).get("hits", []):
        md = h.get("metadata", {})
        add(f"zenodo:{h['id']}", "Zenodo", md.get("title", ""), h.get("links", {}).get("html", f"https://doi.org/10.5281/zenodo.{h['id']}"), "metadata-mined")
    print("zenodo rows so far:", len(rows))
except Exception as e:
    print("zenodo fail:", e)

# 3) CELLxGENE census collections (public API)
try:
    c = json.loads(get("https://api.cellxgene.cziscience.com/curation/v1/collections"))
    n = 0
    for coll in c[:40]:
        cid = coll.get("collection_id")
        add(f"cxg:{cid}", "CELLxGENE", coll.get("name", ""), f"https://cellxgene.cziscience.com/collections/{cid}", "metadata-mined")
        n += 1
    print("cellxgene rows:", n)
except Exception as e:
    print("cellxgene fail:", e)

# 4) GDC projects (TCGA/others - imaging + expression)
try:
    g = json.loads(get("https://api.gdc.cancer.gov/projects?size=80&fields=project_id,name&format=json"))
    for p in g["data"]["hits"]:
        add(p["project_id"], "GDC", p.get("name", ""), f"https://portal.gdc.cancer.gov/projects/{p['project_id']}", "metadata-mined")
    print("gdc rows so far:", len(rows))
except Exception as e:
    print("gdc fail:", e)

# 5) TCIA collections
try:
    t = json.loads(get("https://services.cancerimagingarchive.net/nbia-api/services/v2/getCollectionValues"))
    for coll in t[:40]:
        name = coll.get("criteria") or coll.get("value") or str(coll)
        add(f"tcia:{name}", "TCIA", name, "https://www.cancerimagingarchive.net/collection/" + urllib.parse.quote(str(name)), "metadata-mined")
    print("tcia rows so far:", len(rows))
except Exception as e:
    print("tcia fail:", e)

with open("results/accession_harvest_raw.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["accession", "source", "title", "url", "level_of_use"])
    w.writeheader(); w.writerows(rows)
print("TOTAL harvested:", len(rows))
