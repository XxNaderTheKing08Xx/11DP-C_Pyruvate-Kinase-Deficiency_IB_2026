"""Task 2.1/2.2: characterize the supplied PKLR sequences.

Translates both supplied CDSs (standard genetic code, NCBI table 1), checks the
translations against the supplied proteins, lists every WT-vs-mutant difference at
CDS and protein level, and compares the supplied WT sequences with RefSeq
NM_000298.6 (CDS) / NP_000289.1, Ensembl ENST00000342741.6 / ENSP00000339933.4,
UniProt P30613 and the L-type isoform NM_181871.4 / NP_870986.1.

Run:  python scripts/01_characterize.py   (from the analysis/sequence folder)
Output: results/01_characterize.txt
"""
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq3

BASE = Path(__file__).resolve().parent.parent
SUP, RAW = BASE / "supplied", BASE / "raw"
out_lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out_lines.append(s)


def read1(p):
    return str(next(SeqIO.parse(p, "fasta")).seq).upper()


cds_wt = read1(SUP / "PKD_CDS_WT.fasta")
cds_mu = read1(SUP / "PKD_CDS_mutant.fasta")
pr_wt = read1(SUP / "PKD_protein_WT.fasta")
pr_mu = read1(SUP / "PKD_protein_mutant.fasta")

log("== Supplied files")
for n, s in [("CDS WT", cds_wt), ("CDS mutant", cds_mu), ("protein WT", pr_wt), ("protein mutant", pr_mu)]:
    log(f"{n}: length {len(s)}")
log("CDS lengths divisible by 3:", len(cds_wt) % 3 == 0, len(cds_mu) % 3 == 0)
log("start codon WT/mutant:", cds_wt[:3], cds_mu[:3], "| stop codon:", cds_wt[-3:], cds_mu[-3:])

tr_wt = str(Seq(cds_wt).translate(table=1))
tr_mu = str(Seq(cds_mu).translate(table=1))
log("translation WT ends with stop '*':", tr_wt.endswith("*"), "| internal stops:", tr_wt[:-1].count("*"))
log("translation mutant ends with stop '*':", tr_mu.endswith("*"), "| internal stops:", tr_mu[:-1].count("*"))
log("translation(CDS WT) == supplied protein WT:", tr_wt.rstrip("*") == pr_wt)
log("translation(CDS mutant) == supplied protein mutant:", tr_mu.rstrip("*") == pr_mu)

log("\n== WT vs mutant differences")
nt_diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(cds_wt, cds_mu)) if a != b]
aa_diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(pr_wt, pr_mu)) if a != b]
log("number of nucleotide differences:", len(nt_diffs), nt_diffs)
log("number of amino-acid differences:", len(aa_diffs), aa_diffs)
for pos, a, b in nt_diffs:
    cod = (pos - 1) // 3 + 1
    phase = (pos - 1) % 3 + 1
    cw, cm = cds_wt[(cod - 1) * 3:cod * 3], cds_mu[(cod - 1) * 3:cod * 3]
    aw, am = str(Seq(cw).translate()), str(Seq(cm).translate())
    log(f"c.{pos}{a}>{b}: codon {cod}, codon position {phase}; {cw} -> {cm}; {aw} -> {am}")
    log(f"HGVS protein: p.{seq3(aw)}{cod}{seq3(am)} (1-letter p.{aw}{cod}{am})")
    ctx = 5
    log(f"CDS context c.{pos-ctx}..c.{pos+ctx}: WT {cds_wt[pos-1-ctx:pos-1]}[{cds_wt[pos-1]}]{cds_wt[pos:pos+ctx]}"
        f"  mutant {cds_mu[pos-1-ctx:pos-1]}[{cds_mu[pos-1]}]{cds_mu[pos:pos+ctx]}")
    log(f"dinucleotide c.{pos}-c.{pos+1} in WT: {cds_wt[pos-1:pos+1]} (CpG: {cds_wt[pos-1:pos+1] == 'CG'})")
    log(f"transition (C>T/G>A/A>G/T>C): {(a, b) in [('C','T'),('T','C'),('A','G'),('G','A')]}")

log("\n== Comparison with reference records (identity over full length)")


def ident(a, b):
    if len(a) != len(b):
        return f"lengths differ ({len(a)} vs {len(b)})"
    d = [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]
    return f"length {len(a)}; identical: {not d}; differences: {d}"


refs_nt = {
    "RefSeq NM_000298.6 CDS (44..1768)": read1(RAW / "ncbi/NM_000298.6_cds.fasta"),
    "Ensembl ENST00000342741.6 CDS": read1(RAW / "ensembl/ENST00000342741_cds.fasta"),
}
refs_aa = {
    "RefSeq NP_000289.1": read1(RAW / "ncbi/NP_000289.1.fasta"),
    "UniProt P30613 (canonical, SV=2)": read1(RAW / "uniprot/P30613.fasta"),
    "Ensembl ENSP00000339933.4": read1(RAW / "ensembl/ENSP00000339933.fasta"),
}
for k, v in refs_nt.items():
    log(f"supplied CDS WT vs {k}: {ident(cds_wt, v)}")
    log(f"supplied CDS mutant vs {k}: {ident(cds_mu, v)}")
for k, v in refs_aa.items():
    log(f"supplied protein WT vs {k}: {ident(pr_wt, v)}")

log("\n== R-type vs L-type")
ltype_cds = read1(RAW / "ncbi/NM_181871.4_cds.fasta")
ltype_pr = read1(RAW / "ncbi/NP_870986.1.fasta")
log("L-type NM_181871.4 CDS length:", len(ltype_cds), "| NP_870986.1 length:", len(ltype_pr))
log("R-type (supplied) length:", len(pr_wt))
off = len(pr_wt) - len(ltype_pr)
log("length difference R - L (aa):", off)
# find the longest common suffix between R and L proteins
k = 0
while k < min(len(pr_wt), len(ltype_pr)) and pr_wt[-1 - k] == ltype_pr[-1 - k]:
    k += 1
log(f"longest common C-terminal suffix: {k} aa -> R-type {len(pr_wt)-k+1}..{len(pr_wt)} = L-type {len(ltype_pr)-k+1}..{len(ltype_pr)}")
log(f"R-type-specific N-terminus: R 1..{len(pr_wt)-k} = {pr_wt[:len(pr_wt)-k]}")
log(f"L-type-specific N-terminus: L 1..{len(ltype_pr)-k} = {ltype_pr[:len(ltype_pr)-k]}")
log(f"R-type residue 486 = {pr_wt[485]}; L-type residue 455 = {ltype_pr[454]} (R486 corresponds to L455: offset 31)")
log("Supplied CDS equals R-type CDS:", cds_wt == refs_nt["RefSeq NM_000298.6 CDS (44..1768)"],
    "| equals L-type CDS:", cds_wt == ltype_cds)
log(f"supplied WT protein window 476..496: {pr_wt[475:496]}")

(BASE / "results/01_characterize.txt").write_text("\n".join(out_lines) + "\n")
