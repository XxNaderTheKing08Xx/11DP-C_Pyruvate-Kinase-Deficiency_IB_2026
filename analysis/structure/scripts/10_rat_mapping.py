"""Map human R-type residue 486 onto rat L-type PK (P12928-2) and the rat structures 6ECH/6ECK by pairwise alignment
(Biopython PairwiseAligner, global, BLOSUM62, gap open -10, gap extend -0.5). Output: data/rat_mapping.txt"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
al = Align.PairwiseAligner(); al.mode = "global"
al.substitution_matrix = substitution_matrices.load("BLOSUM62"); al.open_gap_score = -10; al.extend_gap_score = -0.5
hR = str(next(SeqIO.parse("../sequence/supplied/PKD_protein_WT.fasta", "fasta")).seq)
seqs = {"human_L_P30613-2": str(next(SeqIO.parse("raw/uniprot_P30613-2.fasta", "fasta")).seq),
        "rat_R_P12928-1": str(next(SeqIO.parse("raw/uniprot_P12928.fasta", "fasta")).seq),
        "rat_L_P12928-2": str(next(SeqIO.parse("raw/uniprot_P12928-2.fasta", "fasta")).seq)}
for pid in ("6ECH", "6ECK"):
    d = MMCIF2Dict(f"raw/{pid}.cif")
    seqs[f"{pid}_SEQRES"] = d["_entity_poly.pdbx_seq_one_letter_code_can"][0].replace("\n", "")
lines = [f"human R-type (supplied FASTA = P30613-1) length {len(hR)}; residue 486 = {hR[485]}; window 476-496: {hR[475:496]}"]
def map_pos(a, pos1):
    """position in target (1-based) aligned to query position pos1"""
    qi = pos1 - 1
    for (qs, qe), (ts, te) in zip(*a.aligned):
        if qs <= qi < qe:
            return ts + (qi - qs) + 1
    return None
for name, s in seqs.items():
    a = al.align(hR, s)[0]
    p = map_pos(a, 486)
    ident = sum(1 for (qs, qe), (ts, te) in zip(*a.aligned) for i in range(qe - qs) if hR[qs + i] == s[ts + i])
    alnlen = a.shape[1]
    lines.append(f"{name}: length {len(s)}; aligned position of human R486 = index {p} in this sequence, residue {s[p-1]}; "
                 f"identity {ident}/{sum(qe-qs for qs,qe in a.aligned[0])} aligned pairs; window: {s[p-11:p+10]}")
    if name.endswith("SEQRES"):
        lines.append(f"   note: SEQRES index {p}; author numbering = SEQRES index - 1 (from _struct_ref_seq: seq 2-544 -> auth 1-543) -> author residue {p-1}")
# rat window identity with human R around 486
s = seqs["rat_L_P12928-2"]; a = al.align(hR, s)[0]
lines.append("Alignment human R-type (top) vs rat L-type (bottom), residues around 486:")
txt = str(a).split("\n")
lines += [l for l in txt if l.strip()][:0]
p = map_pos(a, 486)
lines.append(f"human R 470-500: {hR[469:500]}")
lines.append(f"rat   L 439-469: {s[438:469]}")
lines.append("differences in window human R 456-516 vs rat L 425-485 (R-type numbering): " +
             ", ".join(f"{hR[i]}{i+1}{s[i-31]}" for i in range(455, 516) if hR[i] != s[i - 31]))
open("data/rat_mapping.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
