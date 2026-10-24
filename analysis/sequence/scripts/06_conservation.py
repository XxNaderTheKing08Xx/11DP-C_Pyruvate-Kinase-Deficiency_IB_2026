"""Task 3.2/3.9: residue at human PKLR position 486 (R-type numbering) in every panel species.

Method A (primary): Clustal Omega multiple alignments (results/msa_local/*.clustalo.aln, or the
EBI outputs when available) -> column containing human Arg486 -> residue and its own-sequence
number in every species; conservation of the window 476..496 (human numbering).
Method B (independent check): pairwise global alignment (Needle-equivalent, BLOSUM62 10/0.5)
of each sequence against human P30613 -> residue aligned to human 486.
Also reports the codon at the equivalent position in the R-type CDSs of rat, mouse and dog,
and the residue of human PKM (P14618 M2 and P14618-2 M1) aligned to PKLR 486 (paralog context).
Run: python scripts/06_conservation.py
"""
import json
import sys
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

BASE = Path(__file__).resolve().parent.parent
RAW, RES = BASE / "raw", BASE / "results"
HPOS, WIN = 486, 10
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


def read1(p):
    return str(next(SeqIO.parse(p, "fasta")).seq).upper()


meta = {}
for line in (RES / "panel_metadata.tsv").read_text().splitlines()[1:]:
    f = line.split("\t")
    meta[f"{f[0].replace(' ', '_')}_{f[4]}"] = dict(label=f[0], species=f[1], clade=f[2], acc=f[4], length=f[9])


def analyse(alnfile, tag):
    aln = AlignIO.read(alnfile, "clustal")
    ids = [r.id for r in aln]
    hi = next(i for i, r in enumerate(aln) if r.id.startswith("Human"))
    hseq = str(aln[hi].seq)
    # map ungapped human position -> column
    col_of, n = {}, 0
    for c, ch in enumerate(hseq):
        if ch != "-":
            n += 1
            col_of[n] = c
    res = {}
    for r in aln:
        s = str(r.seq)
        pos_at_col, n = {}, 0
        for c, ch in enumerate(s):
            if ch != "-":
                n += 1
                pos_at_col[c] = n
        res[r.id] = (s, pos_at_col)
    c486 = col_of[HPOS]
    log(f"\n== {tag}: {alnfile.name}; {len(aln)} sequences; alignment length {aln.get_alignment_length()}; "
        f"human 486 at column {c486 + 1}")
    rows = []
    for rid in ids:
        s, pc = res[rid]
        aa = s[c486]
        own = pc.get(c486, "gap")
        win = "".join(s[col_of[p]] for p in range(HPOS - WIN, HPOS + WIN + 1))
        m = meta.get(rid, {})
        rows.append(dict(id=rid, label=m.get("label", rid), aa=aa, own=own, win=win))
        log(f"{m.get('label', rid):15s} {m.get('acc', ''):20s} residue {aa} at own position {own}; "
            f"window (human {HPOS-WIN}-{HPOS+WIN} columns) {win}")
    # per-column conservation in window
    log("column conservation in window (human numbering: residue counts):")
    cons_line = []
    for p in range(HPOS - WIN, HPOS + WIN + 1):
        colres = [res[r][0][col_of[p]] for r in ids]
        counts = {}
        for x in colres:
            counts[x] = counts.get(x, 0) + 1
        top = max(counts.values())
        cons_line.append(f"{p}:{''.join(sorted(counts, key=lambda k: -counts[k]))}({top}/{len(ids)})")
    log("  " + "  ".join(cons_line))
    ident_cols = sum(1 for p in range(HPOS - WIN, HPOS + WIN + 1)
                     if len({res[r][0][col_of[p]] for r in ids}) == 1)
    log(f"fully identical columns in window: {ident_cols}/{2*WIN+1}")
    return rows


alns = {}
for panel in ["core_panel", "broad_panel"]:
    ebi = sorted((RAW / f"clustalo/{panel}").glob("*.aln-clustal_num.txt")) if (RAW / f"clustalo/{panel}").exists() else []
    f = ebi[-1] if ebi else RES / f"msa_local/{panel}.clustalo.aln"
    alns[panel] = analyse(f, panel)

# Method B: pairwise
al = PairwiseAligner()
al.mode = "global"
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.open_end_gap_score = al.extend_end_gap_score = 0
human = read1(RAW / "uniprot/P30613.fasta")


def residue_at_human(seq, hpos=HPOS):
    a = al.align(human, seq)[0]
    # a.coordinates-based mapping
    for (hs, he), (ss, se) in zip(*a.aligned):
        if hs <= hpos - 1 < he:
            j = ss + (hpos - 1 - hs)
            return seq[j], j + 1
    return "-", "gap"


log("\n== Method B (pairwise to human P30613), broad panel")
panel_seqs = {r.id: str(r.seq) for r in SeqIO.parse(RES / "broad_panel.fasta", "fasta")}
msa_b = {r["id"]: (r["aa"], r["own"]) for r in alns["broad_panel"]}
agree = True
for rid, s in panel_seqs.items():
    aa, own = residue_at_human(s)
    ok = (aa, own) == msa_b[rid]
    agree &= ok
    log(f"{meta[rid]['label']:15s} pairwise: {aa}{own}   MSA: {msa_b[rid][0]}{msa_b[rid][1]}   agree: {ok}")
log("Method A and B agree for all broad-panel sequences:", agree)

log("\n== Paralog context (not part of the panels)")
for lab, f in [("human PKM isoform M2 (P14618)", RAW / "uniprot/P14618.fasta"),
               ("human PKM isoform M1 (P14618-2)", RAW / "uniprot/P14618-2.fasta")]:
    s = read1(f)
    aa, own = residue_at_human(s)
    a = al.align(human, s)[0]
    log(f"{lab}: residue aligned to PKLR Arg486 = {aa}{own}")

log("\n== Codon at the position equivalent to human c.1456-1458 (R-type CDSs; all 1725 nt, 574 codons, no indels)")
for lab, f in [("Human NM_000298.6", RAW / "ncbi/NM_000298.6_cds.fasta"),
               ("Rat XM_006232593.5", RAW / "ncbi/XM_006232593.5_cds.fasta"),
               ("Mouse NM_013631.2", RAW / "ncbi/NM_013631.2_cds.fasta"),
               ("Dog NM_001256262.2", RAW / "ncbi/NM_001256262.2_cds.fasta")]:
    c = read1(f)
    cod = c[(HPOS - 1) * 3:HPOS * 3]
    log(f"{lab}: codon 486 = {cod} ({Seq(cod).translate()}); context codons 484-488: "
        f"{' '.join(c[(i-1)*3:i*3] for i in range(484, 489))}")
log("L-type numbering: human R-type 486 = L-type 455 (NP_870986.1); rat R-type 486 = L-type 455 (NP_036756.3)")
rat_l = read1(RAW / "ncbi/NP_036756.3.fasta")
hum_l = read1(RAW / "ncbi/NP_870986.1.fasta")
log(f"check: NP_870986.1 residue 455 = {hum_l[454]}; NP_036756.3 residue 455 = {rat_l[454]}")

json.dump(alns, open(RES / "residue486_by_species.json", "w"), indent=1)
(RES / "06_conservation.txt").write_text("\n".join(out) + "\n")
