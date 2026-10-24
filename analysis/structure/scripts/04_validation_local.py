"""Per-residue wwPDB validation metrics (RSCC, RSRZ, rotamer, Ramachandran, clashes) around residue 486.
Usage: python -I scripts/04_validation_local.py 2vgb 2vgi"""
import sys, gzip, csv
import xml.etree.ElementTree as ET
RES = [361, 362, 363, 364, 441, 445, 453, 482, 483, 484, 485, 486, 487, 488]
rows = []
for pid in sys.argv[1:]:
    root = ET.parse(gzip.open(f"raw/{pid}_validation.xml.gz")).getroot()
    e = root.find("Entry")
    print(pid.upper(), "clashscore", e.get("clashscore"), "rama_outliers%", e.get("percent-rama-outliers"),
          "rota_outliers%", e.get("percent-rota-outliers"), "RSRZ_outliers%", e.get("percent-RSRZ-outliers"),
          "DCC_Rfree", e.get("DCC_Rfree"), "EDS_R", e.get("EDS_R"), "Wilson B", e.get("WilsonBestimate"))
    partners = {}
    for g in root.findall("ModelledSubgroup"):
        for c in g.findall("clash"):
            partners.setdefault(c.get("cid"), []).append(f"{g.get('chain')}:{g.get('resname')}{g.get('resnum')}:{c.get('atom')}")
    for g in root.findall("ModelledSubgroup"):
        try: n = int(g.get("resnum"))
        except: continue
        if n in RES and g.get("model") == "1":
            me = f"{g.get('chain')}:{g.get('resname')}{g.get('resnum')}"
            clashes = [f"{c.get('atom')}--{'/'.join(x for x in partners[c.get('cid')] if not x.startswith(me + ':' + c.get('atom')))} d={c.get('dist')} overlap={c.get('clashmag')}" for c in g.findall("clash")]
            rows.append([pid.upper(), g.get("chain"), g.get("resname"), n, g.get("rscc"), g.get("rsrz"),
                         g.get("rota"), g.get("rama"), g.get("avgoccu"), ";".join(clashes)])
with open("data/validation_local.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["pdb", "chain", "resname", "resnum", "RSCC", "RSRZ", "rotamer", "ramachandran", "avg_occupancy", "clashes(atom--partner d=distance_A overlap=A)"])
    w.writerows(rows)
for r in rows:
    if r[3] == 486 or r[9]: print("\t".join(map(str, r)))
