"""Retrieve ClinVar records for PKLR p.Arg486Trp and p.Arg486Gln via NCBI E-utilities.

Usage: python scripts/01_clinvar.py .
Run from the analysis/evidence folder. Writes raw responses to raw/clinvar_*.
"""
import json
import sys
import time
from pathlib import Path

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
out = Path(sys.argv[1]) / "raw"
out.mkdir(parents=True, exist_ok=True)

queries = {
    "R486W": 'PKLR[gene] AND ("p.Arg486Trp" OR "c.1456C>T")',
    "R486Q": 'PKLR[gene] AND ("p.Arg486Gln" OR "c.1457G>A")',
}
for tag, term in queries.items():
    r = requests.get(EUTILS + "esearch.fcgi",
                     params={"db": "clinvar", "term": term, "retmode": "json", "retmax": 50})
    r.raise_for_status()
    (out / f"clinvar_esearch_{tag}.json").write_text(r.text)
    ids = r.json()["esearchresult"]["idlist"]
    print(tag, term, "->", ids)
    time.sleep(0.5)
    if not ids:
        continue
    s = requests.get(EUTILS + "esummary.fcgi",
                     params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"})
    s.raise_for_status()
    (out / f"clinvar_esummary_{tag}.json").write_text(s.text)
    time.sleep(0.5)
    for vid in ids:
        f = requests.get(EUTILS + "efetch.fcgi",
                         params={"db": "clinvar", "id": vid, "rettype": "vcv", "is_variationid": "true"})
        f.raise_for_status()
        (out / f"clinvar_efetch_vcv_{tag}_{vid}.xml").write_text(f.text)
        time.sleep(0.5)
