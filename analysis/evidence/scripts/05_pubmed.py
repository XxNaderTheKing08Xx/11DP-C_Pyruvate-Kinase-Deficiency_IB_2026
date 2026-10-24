"""PubMed searches (NCBI E-utilities) and record retrieval for the R486W literature.

Usage:
  python scripts/05_pubmed.py search .        # run the logged searches
  python scripts/05_pubmed.py fetch . PMID ... # save MEDLINE records + abstracts
Run from the analysis/evidence folder. Outputs go to raw/pubmed_*.
"""
import json
import sys
import time
from pathlib import Path

import requests

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
mode, out = sys.argv[1], Path(sys.argv[2]) / "raw"

SEARCHES = {
    "S1_R486W": '(PKLR OR "pyruvate kinase") AND (R486W OR "Arg486Trp" OR "1456C>T" OR "1456T" OR "1456 C>T" OR "1456C-->T")',
    "S2_cohorts": '"pyruvate kinase deficiency"[tiab] AND (mutation*[tiab] OR variant*[tiab]) AND (Italy OR Italian OR Spain OR Spanish OR Portugal OR Portuguese OR France OR French OR "southern Europe")',
    "S3_mitapivat_genotype": '(mitapivat OR AG-348 OR AG348) AND "pyruvate kinase" AND (genotype OR missense OR PKLR)',
    "S4_gene_therapy": '"pyruvate kinase deficiency" AND ("gene therapy" OR lentiviral OR RP-L301)',
    "S5_prevalence_gene_frequency": '"pyruvate kinase deficiency" AND prevalence AND (gene frequency OR allele frequency)',
}

if mode == "search":
    log = {}
    for tag, term in SEARCHES.items():
        r = requests.get(EU + "esearch.fcgi", params={"db": "pubmed", "term": term, "retmode": "json", "retmax": 200})
        r.raise_for_status()
        res = r.json()["esearchresult"]
        log[tag] = {"term": term, "count": res["count"], "pmids": res["idlist"], "querytranslation": res.get("querytranslation")}
        time.sleep(0.4)
        if res["idlist"]:
            s = requests.get(EU + "esummary.fcgi", params={"db": "pubmed", "id": ",".join(res["idlist"]), "retmode": "json"})
            s.raise_for_status()
            (out / f"pubmed_esummary_{tag}.json").write_text(s.text)
            summ = s.json()["result"]
            lines = []
            for p in res["idlist"]:
                x = summ[p]
                au = x["authors"][0]["name"] if x["authors"] else ""
                lines.append(f"{p}\t{x.get('pubdate','')[:4]}\t{au}\t{x.get('source','')}\t{x.get('title','')}")
            (out / f"pubmed_hits_{tag}.tsv").write_text("\n".join(lines) + "\n")
            time.sleep(0.4)
    (out / "pubmed_search_log.json").write_text(json.dumps(log, indent=1))
    for k, v in log.items():
        print(k, v["count"])
elif mode == "fetch":
    ids = sys.argv[3:]
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        r = requests.get(EU + "efetch.fcgi", params={"db": "pubmed", "id": ",".join(chunk), "rettype": "medline", "retmode": "text"})
        r.raise_for_status()
        for rec in r.text.strip().split("\n\n"):
            pmid = [l for l in rec.splitlines() if l.startswith("PMID-")]
            if pmid:
                p = pmid[0].split("-", 1)[1].strip()
                (out / f"pubmed_medline_{p}.txt").write_text(rec + "\n")
                print("saved", p)
        time.sleep(0.4)
