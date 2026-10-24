"""Task 3.1/3.2: multiple sequence alignments with Clustal Omega (EMBL-EBI Job Dispatcher REST API).

Submits results/core_panel.fasta and results/broad_panel.fasta to
https://www.ebi.ac.uk/Tools/services/rest/clustalo with default Clustal Omega parameters
(stype=protein, guidetreeout=true, dismax/mbed defaults, iterations=0) and outfmt=clustal_num,
polls until finished, and saves every result type to raw/clustalo/<panel>/<jobid>.<type>.
The job IDs are appended to raw/clustalo/jobs.tsv.
The e-mail field is mandatory for the service; a non-personal placeholder is used.
Run: python scripts/04_clustalo.py
"""
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import requests

BASE = Path(__file__).resolve().parent.parent
URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
EMAIL = "noreply@example.org"

panels = sys.argv[1:] or ["core_panel", "broad_panel"]
for panel in panels:
    seq = (BASE / f"results/{panel}.fasta").read_text()
    params = {"email": EMAIL, "title": f"PKLR_{panel}", "sequence": seq, "stype": "protein",
              "outfmt": "clustal_num", "guidetreeout": "true", "dismatout": "true", "order": "aligned"}
    r = requests.post(f"{URL}/run", data=params)
    r.raise_for_status()
    job = r.text.strip()
    print(panel, "job", job)
    while True:
        st = requests.get(f"{URL}/status/{job}").text.strip()
        if st in ("FINISHED", "ERROR", "FAILURE", "NOT_FOUND"):
            break
        time.sleep(5)
    print(panel, "status", st)
    outdir = BASE / f"raw/clustalo/{panel}"
    outdir.mkdir(parents=True, exist_ok=True)
    types = requests.get(f"{URL}/resulttypes/{job}", headers={"Accept": "application/json"}).json()
    for t in types["types"]:
        ident = t["identifier"]
        res = requests.get(f"{URL}/result/{job}/{ident}")
        if res.ok:
            (outdir / f"{job}.{ident}.txt").write_text(res.text)
            print("  saved", ident)
    with open(BASE / "raw/clustalo/jobs.tsv", "a") as fh:
        fh.write(f"{panel}\t{job}\t{st}\t{datetime.now(timezone.utc).isoformat(timespec='seconds')}\t"
                 f"stype=protein; outfmt=clustal_num; guidetreeout=true; dismatout=true; order=aligned; other parameters default\n")
