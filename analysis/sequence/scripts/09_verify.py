"""Task 2.8/3.9: independent recomputation of the key numbers by different methods.

1. Variant: second detection of the WT/mutant difference by difflib (string matcher), independent
   of the positional comparison in 01_characterize.py.
2. Pairwise identities: (a) Biopython global aligner (05_pairwise.py), (b) EMBOSS needle (05),
   (c) manual position-by-position count (valid because all four R-type proteins are 574 aa and
   the global alignments contain no gaps; all four R-type CDSs are 1725 nt), (d) Clustal Omega
   percent-identity matrix of the core MSA (identical pairs / aligned non-gap pairs),
   (e) NCBI BLAST+ blastp/blastn bl2seq (local, BLAST+ 2.17.0) when the binary path is given.
3. Residue at human 486: already cross-checked MSA vs pairwise in 06_conservation.py and with
   RefSeq models for macaque/marmoset in 08; here re-read directly from the raw UniProt/RefSeq
   sequences by index (rat/mouse/dog R-type 486; rat L-type 455).
Run: python scripts/09_verify.py <BLAST+ bin dir or 'none'>
"""
import difflib
import re
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO

B = Path(__file__).resolve().parent.parent
RAW, RES = B / "raw", B / "results"
BL = sys.argv[1] if len(sys.argv) > 1 else "none"
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


def read1(p):
    return str(next(SeqIO.parse(p, "fasta")).seq).upper()


# 1. variant via difflib
w, m = read1(B / "supplied/PKD_CDS_WT.fasta"), read1(B / "supplied/PKD_CDS_mutant.fasta")
ops = [op for op in difflib.SequenceMatcher(None, w, m, autojunk=False).get_opcodes() if op[0] != "equal"]
log("difflib opcodes (CDS):", ops, "-> 1-based positions", [(o[0], o[1] + 1, w[o[1]:o[2]], m[o[3]:o[4]]) for o in ops])
pw, pm = read1(B / "supplied/PKD_protein_WT.fasta"), read1(B / "supplied/PKD_protein_mutant.fasta")
ops = [op for op in difflib.SequenceMatcher(None, pw, pm, autojunk=False).get_opcodes() if op[0] != "equal"]
log("difflib opcodes (protein):", [(o[0], o[1] + 1, pw[o[1]:o[2]], pm[o[3]:o[4]]) for o in ops])

# 2c. manual counts
prot = {"Human": RAW / "uniprot/P30613.fasta", "Rat": RAW / "uniprot/P12928.fasta",
        "Mouse": RAW / "uniprot/P53657.fasta", "Dog": RAW / "uniprot/Q29536.fasta"}
cds = {"Human": RAW / "ncbi/NM_000298.6_cds.fasta", "Rat": RAW / "ncbi/XM_006232593.5_cds.fasta",
       "Mouse": RAW / "ncbi/NM_013631.2_cds.fasta", "Dog": RAW / "ncbi/NM_001256262.2_cds.fasta"}
P = {k: read1(v) for k, v in prot.items()}
C = {k: read1(v) for k, v in cds.items()}
log("\nManual position-by-position counts (no gaps; equal lengths):")
manual = {}
for sp in ["Rat", "Mouse", "Dog"]:
    ip = sum(a == b for a, b in zip(P["Human"], P[sp]))
    ic = sum(a == b for a, b in zip(C["Human"], C[sp]))
    syn = sum(1 for i in range(574) if C["Human"][3*i:3*i+3] != C[sp][3*i:3*i+3] and P["Human"][i] == P[sp][i])
    manual[sp] = (ip, ic)
    log(f"Human vs {sp}: protein {ip}/574 identical ({100*ip/574:.1f}%), {574-ip} differences; "
        f"CDS {ic}/1725 identical ({100*ic/1725:.1f}%), {1725-ic} nucleotide differences; "
        f"codons that differ without amino-acid change (vs Swiss-Prot protein): {syn}")

# 2b/2a from files written by 05
txt = (RES / "05_pairwise.txt").read_text()
log("\nNeedle/Biopython (from results/05_pairwise.txt):")
for line in txt.splitlines():
    if line.startswith(("protein Human", "CDS Human")):
        log("  " + line)

# 2d. Clustal Omega distmat
dm = (RES / "msa_local/core_panel.distmat").read_text().split("\n")[1:5]
log("\nClustal Omega 1.2.4 percent identity (core MSA):")
for line in dm:
    f = line.split()
    log(f"  {f[0]}: " + " ".join(f[1:]))

# 2e. BLAST+ bl2seq
if BL != "none":
    log("\nBLAST+ bl2seq (blastp BLOSUM62 11/1 defaults; blastn defaults; best HSP):")
    for sp in ["Rat", "Mouse", "Dog"]:
        for prog, d in [("blastp", prot), ("blastn", cds)]:
            r = subprocess.run([f"{BL}/{prog}", "-query", str(d["Human"]), "-subject", str(d[sp]),
                                "-outfmt", "6 pident length mismatch gapopen qstart qend sstart send evalue bitscore nident"],
                               capture_output=True, text=True, check=True)
            first = r.stdout.strip().splitlines()[0].split("\t")
            log(f"  {prog} Human vs {sp}: pident {first[0]}%, HSP length {first[1]}, identical {first[10]}, "
                f"mismatches {first[2]}, gap openings {first[3]}, query {first[4]}-{first[5]}, subject {first[6]}-{first[7]}, E {first[8]}")

# 3. direct indexing
log("\nDirect indexing of raw sequences at 486 (R-type) / 455 (L-type):")
for sp, s in P.items():
    log(f"  {sp} Swiss-Prot residue 486 = {s[485]}; codon 486 in R-type CDS = {C[sp][1455:1458]}")
log(f"  Rat NP_036756.3 (L-type) residue 455 = {read1(RAW / 'ncbi/NP_036756.3.fasta')[454]}")
log(f"  Human NP_870986.1 (L-type) residue 455 = {read1(RAW / 'ncbi/NP_870986.1.fasta')[454]}")
log(f"  Mouse NP_038659.2 residue 486 = {read1(RAW / 'ncbi/NP_038659.2.fasta')[485]}; "
    f"Dog NP_001243191.1 residue 486 = {read1(RAW / 'ncbi/NP_001243191.1.fasta')[485]}")
(RES / "09_verify.txt").write_text("\n".join(out) + "\n")
