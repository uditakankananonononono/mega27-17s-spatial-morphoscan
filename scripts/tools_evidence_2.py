"""Supplemental JSON-POST tools (need proper JSON bodies)."""
import json, time, urllib.parse, urllib.request, csv, os
OUT = "results/tools"
UA = {"User-Agent": "morphoscan-research/0.1 (academic)", "Content-Type": "application/json"}
GENES = ["ERBB2", "ESR1", "MKI67", "PGR", "GATA3", "KRT19", "KRT5", "PTK6", "GRB7", "STARD3"]
rows = []

def call(idx, name, category, purpose, url, payload, extract=None, timeout=60):
    try:
        req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers=UA, method="POST")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            raw = r.read()
        path = os.path.join(OUT, f"{idx:02d}_{name.replace(' ','_').replace('/','-')}.json")
        open(path, "wb").write(raw[:400000])
        n = None
        try:
            if extract: n = extract(json.loads(raw))
        except Exception: pass
        rows.append(dict(tool=name, category=category, purpose=purpose, status="certified", detail=f"evidence {path}", records=n if n is not None else ""))
        print("OK", name, n)
    except Exception as e:
        rows.append(dict(tool=name, category=category, purpose=purpose, status="attempted-failed", detail=repr(e)[:180], records=""))
        print("FAIL", name, repr(e)[:100])
    time.sleep(0.4)

call(60, "g:Profiler g:GOSt", "enrichment", "GO/KEGG/Reactome enrichment of 10 key breast genes",
     "https://biit.cs.ut.ee/gprofiler/api/gost/profile",
     {"organism": "hsapiens", "query": GENES, "sources": ["GO:BP", "KEGG", "REAC"], "no_evidences": True},
     lambda j: len(j.get("result", [])))
call(61, "Open Targets GraphQL", "drug-db", "approved drugs targeting ERBB2",
     "https://api.platform.opentargets.org/api/v4/graphql",
     {"query": "{ target(ensemblId: \"ENSG00000141736\") { approvedSymbol knownDrugs { count rows { drug { name drugType maximumClinicalTrialPhase } } } } }"},
     lambda j: j["data"]["target"]["knownDrugs"]["count"])
call(62, "Figshare API", "data-repository", "figshare articles on spatial transcriptomics",
     "https://api.figshare.com/v2/articles/search",
     {"search_for": "spatial transcriptomics", "limit": 15},
     lambda j: len(j))
call(63, "Enrichr addList+enrich", "enrichment", "Enrichr GO BP enrichment of key genes",
     "https://maayanlab.cloud/Enrichr/addList",
     None, timeout=45) if False else None
# Enrichr needs multipart; use GET-based userListId flow instead
try:
    import urllib.parse as up
    data = up.urlencode({"list": "\n".join(GENES), "description": "morphoscan key genes"}).encode()
    req = urllib.request.Request("https://maayanlab.cloud/Enrichr/addList", data=data,
                                 headers={"User-Agent": UA["User-Agent"]})
    j = json.loads(urllib.request.urlopen(req, timeout=45).read())
    uid = j["userListId"]
    j2 = json.loads(urllib.request.urlopen(f"https://maayanlab.cloud/Enrichr/enrich?userListId={uid}&backgroundType=GO_Biological_Process_2023", timeout=45).read())
    open(os.path.join(OUT, "63_Enrichr_GO_BP.json"), "w").write(json.dumps(j2)[:400000])
    n = len(j2.get("GO_Biological_Process_2023", []))
    rows.append(dict(tool="Enrichr enrichment", category="enrichment", purpose="GO BP enrichment of key genes", status="certified", detail="evidence results/tools/63_Enrichr_GO_BP.json", records=n))
    print("OK Enrichr enrich", n)
except Exception as e:
    rows.append(dict(tool="Enrichr enrichment", category="enrichment", purpose="GO BP enrichment of key genes", status="attempted-failed", detail=repr(e)[:180], records=""))
    print("FAIL Enrichr", repr(e)[:100])

def getj(url, timeout=45):
    req = urllib.request.Request(url, headers={"User-Agent": "morphoscan-research/0.1"})
    return json.loads(urllib.request.urlopen(req, timeout=timeout).read())

for idx, name, cat, purp, url, ext in [
    (64, "bioRxiv details API", "literature", "recent biorxiv submissions scanned for spatial transcriptomics",
     "https://api.biorxiv.org/details/biorxiv/2026-09-01/2026-09-25/0", lambda j: len(j.get("collection", []))),
    (65, "OSF Preprints", "literature", "preprints on spatial transcriptomics",
     "https://api.osf.io/v2/preprints/?filter[q]=spatial%20transcriptomics&page[size]=15", lambda j: len(j.get("data", []))),
    (66, "EBI EB-eye search", "search-index", "cross-resource search: breast cancer expression atlas entries",
     "https://www.ebi.ac.uk/ebisearch/ws/rest/gx?query=breast%20cancer&size=15&format=json", lambda j: j.get("hitCount", 0)),
]:
    try:
        j = getj(url)
        open(os.path.join(OUT, f"{idx:02d}_{name.replace(' ','_')}.json"), "w").write(json.dumps(j)[:400000])
        rows.append(dict(tool=name, category=cat, purpose=purp, status="certified", detail=f"evidence results/tools/{idx:02d}_{name.replace(' ','_')}.json", records=ext(j)))
        print("OK", name)
    except Exception as e:
        rows.append(dict(tool=name, category=cat, purpose=purp, status="attempted-failed", detail=repr(e)[:180], records=""))
        print("FAIL", name, repr(e)[:100])
    time.sleep(0.4)

with open("results/tools_ledger_raw2.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tool", "category", "purpose", "status", "detail", "records"])
    w.writeheader(); w.writerows(rows)
print("supplement certified:", sum(1 for r in rows if r["status"] == "certified"), "/", len(rows))
