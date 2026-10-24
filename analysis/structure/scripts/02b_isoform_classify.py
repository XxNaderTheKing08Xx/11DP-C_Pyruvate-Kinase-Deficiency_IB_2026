"""Classify each RCSB entry of the survey by the isoform-specific N-terminus present in its construct sequence:
R-type-specific residues 1-33 (MSIQENISSLQLRSWVSKSQRDLAKSILIGAPG) vs L-type-specific start (MEGPAG...).
Constructs starting after residue 33 contain only the region shared by both isoforms. Output: data/rcsb_survey_isoform.tsv"""
import json
d = json.load(open("data/rcsb_survey.json"))
rows = []
for e in d:
    for p in e["polymer_entities"]:
        s = p["entity_poly"]["pdbx_seq_one_letter_code_can"]
        if "SIQENISSLQ" in s or "SVQENTLPQQ" in s: iso = "R-type N-terminus present"
        elif "MEGPAGYLRR" in s or "EGPAGYLRR" in s: iso = "L-type N-terminus present"
        else: iso = "shared region only (construct lacks isoform-specific N-terminus)"
        rows.append((e["rcsb_id"], p["rcsb_polymer_entity"]["pdbx_description"], len(s), s[:20], iso, e["exptl"][0]["method"],
                     (e["rcsb_entry_info"]["resolution_combined"] or [None])[0], e["struct"]["title"]))
with open("data/rcsb_survey_isoform.tsv", "w") as f:
    f.write("pdb\tentity_description\tseq_length\tseq_start\tisoform_evidence\tmethod\tresolution_A\ttitle\n")
    for r in sorted(rows): f.write("\t".join(map(str, r)) + "\n")
from collections import Counter
print(Counter((r[4], r[5]) for r in rows))
for r in sorted(rows):
    if not r[4].startswith("L-type"): print(r[0], r[2], r[3], r[4], r[5], r[6], r[7][:80])
