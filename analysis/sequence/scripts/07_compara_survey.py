"""Survey: residue aligned to human PKLR Arg486 in every Ensembl Compara (release 116) orthologue
of ENSG00000143627 (pairwise Needle-equivalent alignment to P30613, BLOSUM62, 10/0.5, free end gaps).
Output: results/compara_residue486_survey.tsv and a summary on stdout."""
import json
from pathlib import Path
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
BASE = Path(__file__).resolve().parent.parent
human = str(next(SeqIO.parse(BASE / "raw/uniprot/P30613.fasta", "fasta")).seq)
al = PairwiseAligner(); al.mode = "global"
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5
al.open_end_gap_score = al.extend_end_gap_score = 0
hs = json.loads((BASE / "raw/ensembl/compara_PKLR_orthologues.json").read_text())["data"][0]["homologies"]
rows = []
for h in hs:
    t = h["target"]; s = t["seq"].replace("-", "")
    a = al.align(human, s)[0]
    aa, pos, win = "-", "gap", ""
    for (x0, x1), (y0, y1) in zip(*a.aligned):
        if x0 <= 485 < x1:
            j = y0 + 485 - x0; aa, pos = s[j], j + 1
    # window 476-496 as aligned to human
    m = {}
    for (x0, x1), (y0, y1) in zip(*a.aligned):
        for k in range(x1 - x0):
            m[x0 + k] = s[y0 + k]
    win = "".join(m.get(i, "-") for i in range(475, 496))
    rows.append((h["taxonomy_level"], t["species"], h["type"], t["id"], t["protein_id"], len(s), round(t["perc_id"], 1), aa, pos, win))
with open(BASE / "results/compara_residue486_survey.tsv", "w") as fh:
    fh.write("taxonomy_level\tspecies\thomology_type\tgene\tprotein\tlength\tperc_id_target\tres486\town_pos\twindow476_496\n")
    for r in sorted(rows, key=lambda r: (r[0], r[1])):
        fh.write("\t".join(map(str, r)) + "\n")
print(Counter(r[7] for r in rows))
for r in sorted(rows):
    if r[7] != "R": print(r)
