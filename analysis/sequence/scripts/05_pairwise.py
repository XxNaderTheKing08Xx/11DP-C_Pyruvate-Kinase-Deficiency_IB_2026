"""Task 3.1/3.2 (+3.9): pairwise global alignments human vs rat, mouse and dog PKLR.

Primary method: Biopython PairwiseAligner, global (Needleman-Wunsch), parameters equal to
EMBOSS Needle defaults:
  protein: BLOSUM62, gap open 10, gap extend 0.5, end gaps not penalized
  CDS:     EDNAFULL (= NUC.4.4), gap open 10, gap extend 0.5, end gaps not penalized
Statistics follow the EMBOSS definitions: identity = identical aligned pairs / alignment
length (columns, gaps included); similarity = pairs with positive substitution score /
alignment length; gaps = gap columns / alignment length.
Independent check: EMBOSS needle 6.6.0 (same parameters) run on the same pairs; its report
files are saved in results/needle/ and parsed for comparison.
Sequences: Swiss-Prot R-type proteins (P30613, P12928, P53657, Q29536) and R-type CDSs
(NM_000298.6, XM_006232593.5 [rat; the only rat RefSeq R-type CDS], NM_013631.2, NM_001256262.2).
Run: python scripts/05_pairwise.py <path to EMBOSS bin dir or 'none'>
"""
import os
import re
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

BASE = Path(__file__).resolve().parent.parent
RAW, RES = BASE / "raw", BASE / "results"
EMBOSS = sys.argv[1] if len(sys.argv) > 1 else "none"
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


def read1(p):
    return str(next(SeqIO.parse(p, "fasta")).seq).upper()


def make_aligner(matrix):
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load(matrix)
    al.open_gap_score = -10
    al.extend_gap_score = -0.5
    al.open_end_gap_score = 0
    al.extend_end_gap_score = 0
    return al


def stats(al, a, b):
    aln = al.align(a, b)[0]
    s1, s2 = aln[0], aln[1]
    m = al.substitution_matrix
    L = len(s1)
    ident = sim = gaps = 0
    for x, y in zip(s1, s2):
        if x == "-" or y == "-":
            gaps += 1
        else:
            ident += x == y
            sim += m[x][y] > 0
    return dict(length=L, identity=ident, similarity=sim, gaps=gaps, score=aln.score,
                pid=100 * ident / L, psim=100 * sim / L, pgap=100 * gaps / L, s1=s1, s2=s2)


prot = {
    "Human P30613": read1(RAW / "uniprot/P30613.fasta"),
    "Rat P12928": read1(RAW / "uniprot/P12928.fasta"),
    "Mouse P53657": read1(RAW / "uniprot/P53657.fasta"),
    "Dog Q29536": read1(RAW / "uniprot/Q29536.fasta"),
}
cds = {
    "Human NM_000298.6": read1(RAW / "ncbi/NM_000298.6_cds.fasta"),
    "Rat XM_006232593.5": read1(RAW / "ncbi/XM_006232593.5_cds.fasta"),
    "Mouse NM_013631.2": read1(RAW / "ncbi/NM_013631.2_cds.fasta"),
    "Dog NM_001256262.2": read1(RAW / "ncbi/NM_001256262.2_cds.fasta"),
}
log("== CDS translation checks (standard code)")
for (k, v), (kp, p) in zip(cds.items(), prot.items()):
    t = str(Seq(v).translate(table=1)).rstrip("*")
    d = [(i + 1, x, y) for i, (x, y) in enumerate(zip(t, p)) if x != y] if len(t) == len(p) else "length differs"
    log(f"{k}: {len(v)} nt; translation {len(t)} aa; equals {kp}: {t == p}; differences (translation vs Swiss-Prot): {d}")

pa, na = make_aligner("BLOSUM62"), make_aligner("NUC.4.4")
rows = []
log("\n== Pairwise global alignments (Biopython PairwiseAligner, Needle-equivalent)")
hp, hc = "Human P30613", "Human NM_000298.6"
for (kp, p), (kc, c) in zip(list(prot.items())[1:], list(cds.items())[1:]):
    sp_ = stats(pa, prot[hp], p)
    sc_ = stats(na, cds[hc], c)
    log(f"protein {hp} vs {kp}: length {sp_['length']}, identity {sp_['identity']}/{sp_['length']} ({sp_['pid']:.1f}%), "
        f"similarity {sp_['similarity']}/{sp_['length']} ({sp_['psim']:.1f}%), gaps {sp_['gaps']}/{sp_['length']} ({sp_['pgap']:.1f}%), score {sp_['score']:.1f}")
    log(f"CDS {hc} vs {kc}: length {sc_['length']}, identity {sc_['identity']}/{sc_['length']} ({sc_['pid']:.1f}%), "
        f"gaps {sc_['gaps']}/{sc_['length']} ({sc_['pgap']:.1f}%), score {sc_['score']:.1f}")
    rows.append((kp, kc, sp_, sc_))
    (RES / f"pairwise_protein_human_vs_{kp.split()[0].lower()}.txt").write_text(
        f"# {hp} vs {kp}; BLOSUM62, open 10, extend 0.5, end gaps free\n{sp_['s1']}\n{sp_['s2']}\n")

# rat/mouse/dog among themselves (for completeness)
log("\n== Other protein pairs")
names = list(prot)
for i in range(1, 4):
    for j in range(i + 1, 4):
        s = stats(pa, prot[names[i]], prot[names[j]])
        log(f"{names[i]} vs {names[j]}: identity {s['pid']:.1f}%, similarity {s['psim']:.1f}%, gaps {s['gaps']}")

# EMBOSS needle independent check
if EMBOSS != "none":
    nd = RES / "needle"
    nd.mkdir(exist_ok=True)
    env = dict(os.environ, PATH=EMBOSS + ":" + os.environ["PATH"])
    tmp = nd / "seqs"
    tmp.mkdir(exist_ok=True)
    for k, v in list(prot.items()) + list(cds.items()):
        (tmp / f"{k.replace(' ', '_')}.fasta").write_text(f">{k.replace(' ', '_')}\n{v}\n")
    log("\n== EMBOSS needle 6.6.0 (gapopen 10, gapextend 0.5; protein EBLOSUM62, DNA EDNAFULL; endweight false)")
    agree = True
    for kp, kc, sp_, sc_ in rows:
        for (a, b, st, kind) in [(hp, kp, sp_, "protein"), (hc, kc, sc_, "cds")]:
            o = nd / f"needle_{kind}_{a.split()[0]}_vs_{b.split()[0]}.needle"
            rel = lambda p: os.path.relpath(p, BASE)  # relative paths, so the report headers contain no local paths
            cmd = ["needle", "-asequence", rel(tmp / f"{a.replace(' ', '_')}.fasta"),
                   "-bsequence", rel(tmp / f"{b.replace(' ', '_')}.fasta"),
                   "-gapopen", "10", "-gapextend", "0.5", "-outfile", rel(o), "-auto"]
            subprocess.run(cmd, env=env, check=True, cwd=BASE)
            txt = o.read_text()
            L = int(re.search(r"# Length: (\d+)", txt).group(1))
            I = int(re.search(r"# Identity:\s+(\d+)/", txt).group(1))
            S = int(re.search(r"# Similarity:\s+(\d+)/", txt).group(1))
            G = int(re.search(r"# Gaps:\s+(\d+)/", txt).group(1))
            sc = float(re.search(r"# Score: ([\d.]+)", txt).group(1))
            same = (L, I, G) == (st["length"], st["identity"], st["gaps"]) and (kind == "cds" or S == st["similarity"])
            agree &= same
            log(f"{kind} {a} vs {b}: needle length {L}, identity {I} ({100*I/L:.1f}%), similarity {S} ({100*S/L:.1f}%), "
                f"gaps {G}, score {sc}; Biopython identity {st['identity']}, gaps {st['gaps']}, score {st['score']:.1f} -> agree: {same}")
    log("All needle vs Biopython counts agree:", agree)

with open(RES / "pairwise_summary.tsv", "w") as fh:
    fh.write("pair_protein\tprot_len\tprot_identity\tprot_pid\tprot_similarity\tprot_psim\tprot_gaps\t"
             "pair_cds\tcds_len\tcds_identity\tcds_pid\tcds_gaps\n")
    for kp, kc, sp_, sc_ in rows:
        fh.write(f"Human P30613 vs {kp}\t{sp_['length']}\t{sp_['identity']}\t{sp_['pid']:.1f}\t{sp_['similarity']}\t"
                 f"{sp_['psim']:.1f}\t{sp_['gaps']}\tHuman NM_000298.6 vs {kc}\t{sc_['length']}\t{sc_['identity']}\t"
                 f"{sc_['pid']:.1f}\t{sc_['gaps']}\n")
(RES / "05_pairwise.txt").write_text("\n".join(out) + "\n")
