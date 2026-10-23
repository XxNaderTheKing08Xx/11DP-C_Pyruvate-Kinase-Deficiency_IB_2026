# Independent verification of the sequence characterization (task 2.8)

The sequence characterization (task 2.1) and the pairwise comparisons (task 3.1) were repeated by a second member with different methods, so that no key number depends on a single tool or a single person. Script: `scripts/09_verify.py`; full output: `results/09_verify.txt`.

## What was checked, and how

| Claim | First method (tasks 2.1, 3.1) | Independent check (task 2.8) | Result |
|---|---|---|---|
| The wild-type and mutant CDSs differ at one position only | Position-by-position comparison (`01_characterize.py`) | Python `difflib` on the full CDSs | One substitution, c.1456C>T — agrees |
| The proteins differ at one residue only | Translation and comparison | `difflib` on the supplied protein sequences | One substitution, R486W — agrees |
| Human–rat, human–mouse, human–dog protein identity | Biopython global alignment (Needleman–Wunsch, BLOSUM62, gap opening 10, extension 0.5) | EMBOSS Needle; manual position-by-position count; Clustal Omega identity matrix; BLAST+ blastp | 523/574 (91.1%), 526/574 (91.6%), 531/574 (92.5%) with every method — agree |
| CDS identity | Position-by-position count over the 1725-nucleotide R-type CDSs | EMBOSS Needle; BLAST+ blastn | Rat 1491/1725 (86.4%), mouse 1490/1725 (86.4%), dog 1544/1725 (89.5%). Needle and blastn give 86.5–89.5%, because they insert a few one-nucleotide gaps; the gap-free count is the one reported |
| Residue at human position 486 in each species | Read from the alignment column | Direct indexing of the raw sequences (R-type numbering; L-type = R-type − 31) | Human Arg (CGG), rat Gln (CAA), mouse Arg (CGC), dog Arg (CGG); rat L-type NP_036756.3 Gln455 — agree |

## Additional observations from the check

- Synonymous codon differences (codons that differ without an amino-acid change): rat 152, mouse 155, dog 121. Most nucleotide differences between the orthologous CDSs therefore leave the protein unchanged.
- Clustal Omega also gives rat–mouse protein identity 96.3%, higher than either rodent with human, as expected for the two most closely related species of the panel.

## What a verification of this kind cannot detect

All the methods start from the same downloaded records. An error in the choice of a record (for example, a paralogue taken for an orthologue, or an L-type sequence taken for an R-type one) would be reproduced by every method. That risk is addressed separately, by the record-selection criteria in `README.md` (reviewed Swiss-Prot entries, Ensembl Compara one-to-one orthologues, comparison with human PKM).
