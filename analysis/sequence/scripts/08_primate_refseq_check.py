"""Independent check (RefSeq models, NCBI annotation) of the residue aligned to human PKLR Arg486
in macaque (XM_001112902.5 -> XP_001112902.2) and marmoset (XM_003735138.5, XM_008984562.5),
which carry Gln in the Ensembl 116 models. Reports codon and residue (pairwise Needle-equivalent)."""
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
B = Path(__file__).resolve().parent.parent
h = str(next(SeqIO.parse(B / "raw/uniprot/P30613.fasta", "fasta")).seq)
al = PairwiseAligner(); al.mode = "global"; al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5; al.open_end_gap_score = al.extend_end_gap_score = 0
lines = []
for acc in ["XM_001112902.5", "XM_003735138.5", "XM_008984562.5"]:
    rec = next(SeqIO.parse(B / f"raw/ncbi/{acc}_cds.fasta", "fasta"))
    c = str(rec.seq); p = str(Seq(c).translate()).rstrip("*")
    a = al.align(h, p)[0]
    for (x0, x1), (y0, y1) in zip(*a.aligned):
        if x0 <= 485 < x1:
            j = y0 + 485 - x0
    s = f"{acc} ({rec.description.split('[protein_id=')[1].split(']')[0]}): {len(p)} aa; residue aligned to human 486 = {p[j]}{j+1}; codon {c[3*j:3*j+3]}; window {p[j-10:j+11]}"
    print(s); lines.append(s)
(B / "results/08_primate_refseq_check.txt").write_text("\n".join(lines) + "\n")
