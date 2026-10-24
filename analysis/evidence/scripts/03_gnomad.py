"""Query the gnomAD GraphQL API for PKLR c.1456C>T (R486W) and c.1457G>A (R486Q).

Usage: python scripts/03_gnomad.py .
Variant IDs are GRCh38 (chrom-pos-ref-alt), taken from the ClinVar VCV records
(VCV000001513: NC_000001.11:g.155291918G>A; VCV000873616: NC_000001.11:g.155291917C>T).
Datasets: gnomad_r4 (v4.1, GRCh38) and gnomad_r2_1 (v2.1.1, GRCh37) for comparison.
"""
import json
import sys
import time
from pathlib import Path

import requests

API = "https://gnomad.broadinstitute.org/api"
out = Path(sys.argv[1]) / "raw"

QUERY = """
query V($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id reference_genome chrom pos ref alt rsids caid
    exome { ac an af homozygote_count hemizygote_count filters
            populations { id ac an homozygote_count hemizygote_count } }
    genome { ac an af homozygote_count hemizygote_count filters
             populations { id ac an homozygote_count hemizygote_count } }
    joint { ac an homozygote_count hemizygote_count filters
            populations { id ac an homozygote_count hemizygote_count } }
    transcript_consequences { gene_symbol transcript_id transcript_version hgvsc hgvsp
                              consequence_terms is_canonical is_mane_select }
    in_silico_predictors { id value flags }
  }
  meta { clinvar_release_date }
}
"""
QUERY_R2 = """
query V($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id reference_genome chrom pos ref alt rsids
    exome { ac an af homozygote_count filters populations { id ac an homozygote_count } }
    genome { ac an af homozygote_count filters populations { id ac an homozygote_count } }
  }
}
"""
jobs = [
    ("R486W", "gnomad_r4", "1-155291918-G-A", QUERY),
    ("R486Q", "gnomad_r4", "1-155291917-C-T", QUERY),
    ("R486W", "gnomad_r2_1", "1-155261709-G-A", QUERY_R2),
]
for tag, ds, vid, q in jobs:
    r = requests.post(API, json={"query": q, "variables": {"variantId": vid, "dataset": ds}},
                      headers={"Content-Type": "application/json"}, timeout=60)
    print(tag, ds, vid, r.status_code)
    (out / f"gnomad_{ds}_{tag}_{vid}.json").write_text(json.dumps(r.json(), indent=1))
    time.sleep(2)
