"""Task 3.x: ClinVar records for PKLR codon 486 (p.Arg486Trp and p.Arg486Gln).

Queries NCBI E-utilities (esearch + esummary, db=clinvar), saves raw JSON in raw/clinvar/
and writes results/02_clinvar.txt.
Run: python scripts/02_clinvar.py   (from the analysis/sequence folder)
"""
import json
import time
from pathlib import Path
import requests

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "raw/clinvar"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


queries = {
    "Arg486Trp": 'PKLR[gene] AND "p.Arg486Trp"',
    "Arg486Gln": 'PKLR[gene] AND "p.Arg486Gln"',
    "c.1456C>T": '"NM_000298.6(PKLR):c.1456C>T"',
    "c.1457G>A": '"NM_000298.6(PKLR):c.1457G>A"',
    "codon486_any": 'PKLR[gene] AND (Arg486 OR R486)',
}
all_ids = set()
for name, term in queries.items():
    r = requests.get(f"{EU}/esearch.fcgi", params={"db": "clinvar", "term": term, "retmode": "json", "retmax": 100})
    r.raise_for_status()
    d = r.json()
    (RAW / f"esearch_{name.replace('>', '_')}.json").write_text(json.dumps(d, indent=1))
    ids = d["esearchresult"]["idlist"]
    log(f"esearch [{term}] -> count {d['esearchresult']['count']}, ids {ids}")
    all_ids.update(ids)
    time.sleep(0.5)

if all_ids:
    r = requests.get(f"{EU}/esummary.fcgi", params={"db": "clinvar", "id": ",".join(sorted(all_ids)), "retmode": "json"})
    r.raise_for_status()
    d = r.json()
    (RAW / "esummary_codon486.json").write_text(json.dumps(d, indent=1))
    for uid in d["result"]["uids"]:
        rec = d["result"][uid]
        gc = rec.get("germline_classification", {}) or {}
        log("---")
        log("VariationID:", uid, "| accession:", rec.get("accession"), rec.get("accession_version"))
        log("title:", rec.get("title"))
        log("protein change:", rec.get("protein_change"))
        log("germline classification:", gc.get("description"), "| review status:", gc.get("review_status"),
            "| last evaluated:", gc.get("last_evaluated"))
        traits = [t.get("trait_name") for t in gc.get("trait_set", [])]
        log("conditions:", traits)
        for v in rec.get("variation_set", []):
            log("variation:", v.get("variation_name"), "| cdna:", v.get("cdna_change"),
                "| xrefs:", [(x.get("db_source"), x.get("db_id")) for x in v.get("variation_xrefs", [])])
            for loc in v.get("variation_loc", []):
                if loc.get("assembly_name") == "GRCh38":
                    log("GRCh38 location:", loc.get("chr"), loc.get("start"), loc.get("stop"))
(BASE / "results/02_clinvar.txt").write_text("\n".join(out) + "\n")
