"""Europe PMC metadata (open-access status, PMCID, licence, DOI) and full-text XML for selected PMIDs.

Usage:
  python scripts/06_europepmc.py meta . PMID ...      -> raw/europepmc_meta.tsv (+ json per PMID)
  python scripts/06_europepmc.py fulltext . PMCID ... -> raw/fulltext_<PMCID>.xml
"""
import json
import sys
import time
from pathlib import Path

import requests

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
mode, out = sys.argv[1], Path(sys.argv[2]) / "raw"
if mode == "meta":
    rows = []
    for p in sys.argv[3:]:
        r = requests.get(BASE + "search", params={"query": f"EXT_ID:{p} AND SRC:MED", "resultType": "core", "format": "json"})
        r.raise_for_status()
        res = r.json()["resultList"]["result"]
        (out / f"europepmc_meta_{p}.json").write_text(json.dumps(res, indent=1))
        if not res:
            rows.append(f"{p}\tNOT FOUND")
            continue
        x = res[0]
        rows.append("\t".join(str(v) for v in [
            p, x.get("pmcid", ""), x.get("doi", ""), x.get("isOpenAccess", ""), x.get("license", ""),
            x.get("inEPMC", ""), x.get("hasPDF", ""), x.get("authorString", "")[:60], x.get("pubYear", ""),
            x.get("journalInfo", {}).get("journal", {}).get("isoabbreviation", ""),
            x.get("journalInfo", {}).get("volume", ""), x.get("pageInfo", ""), x.get("title", "")]))
        time.sleep(0.3)
    hdr = "pmid\tpmcid\tdoi\tisOpenAccess\tlicense\tinEPMC\thasPDF\tauthors\tyear\tjournal\tvolume\tpages\ttitle"
    path = out / "europepmc_meta.tsv"
    existing = path.read_text().splitlines()[1:] if path.exists() else []
    keep = [l for l in existing if l.split("\t")[0] not in sys.argv[3:]]
    path.write_text("\n".join([hdr] + keep + rows) + "\n")
    print("\n".join(rows))
elif mode == "fulltext":
    for pmc in sys.argv[3:]:
        r = requests.get(BASE + f"{pmc}/fullTextXML")
        print(pmc, r.status_code, len(r.content))
        if r.ok:
            (out / f"fulltext_{pmc}.xml").write_bytes(r.content)
        time.sleep(0.3)
