"""Summarize gnomAD JSON responses (from 03_gnomad.py) as Markdown tables.

Usage: python scripts/04_gnomad_table.py raw/gnomad_gnomad_r4_R486W_1-155291918-G-A.json
Prints, for exome, genome, and joint (exome + genome) data: AC, AN, AF, homozygotes,
overall and per genetic ancestry group (sex-specific subsets omitted).
"""
import json
import sys

LABELS = {
    "afr": "African/African American", "amr": "Admixed American", "asj": "Ashkenazi Jewish",
    "eas": "East Asian", "fin": "Finnish", "nfe": "Non-Finnish European", "mid": "Middle Eastern",
    "sas": "South Asian", "ami": "Amish", "remaining": "Remaining", "oth": "Other",
}
d = json.load(open(sys.argv[1]))["data"]["variant"]
print(f"Variant {d['variant_id']} ({d['reference_genome']}), rsIDs {d.get('rsids')}")
for tc in d.get("transcript_consequences") or []:
    if tc.get("is_mane_select") or tc.get("is_canonical"):
        print("  ", tc)
for p in d.get("in_silico_predictors") or []:
    print("  predictor", p)
for part in ("exome", "genome", "joint"):
    x = d.get(part)
    if not x:
        print(f"\n{part}: no data")
        continue
    af = x["ac"] / x["an"] if x["an"] else float("nan")
    print(f"\n### {part}: AC {x['ac']}, AN {x['an']}, AF {af:.3e}, homozygotes {x['homozygote_count']}, filters {x.get('filters')}")
    print("| Group | AC | AN | AF | Homozygotes |\n|---|---|---|---|---|")
    for p in x["populations"]:
        if "_" in p["id"] or p["id"] in ("XX", "XY"):
            continue
        a = p["ac"] / p["an"] if p["an"] else float("nan")
        print(f"| {LABELS.get(p['id'], p['id'])} ({p['id']}) | {p['ac']} | {p['an']} | {a:.3e} | {p['homozygote_count']} |")
