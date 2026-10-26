# Analysis procedures · step by step

Every analysis of the report, with the files we started from (`input/`), what we did on each website or program (`procedure/`: numbered screenshots and page prints), and what we obtained (`output/`, `records/`). The steps follow the course labs (P1 Databases, P4 PDB, P5 Sequence alignment, P6 Structural similarity) and our project plan of 9 October.

| Folder | Lab / task | What it contains | Used in the report |
|---|---|---|---|
| [`1_assigned_sequences/`](1_assigned_sequences/) | — | The four FASTA files supplied by the course | All Results |
| [`2_identification/`](2_identification/) | P5 §1.1–1.2; task 2.2 | BLAST (blastn, blastp) of the wild-type and mutant CDS and protein: forms, results, text reports, hit tables; NCBI records NM_000298.6, NP_000289.1, Gene 5313 | Results (identification), Methods |
| [`3_variant/`](3_variant/) | P5 §2.1–2.2; task 2.1 | BLAST "Align two or more sequences", wild type vs mutant (CDS and protein), "Query-anchored with dots for identities" view: c.1456C>T, p.Arg486Trp | Results, Table 1 |
| [`4_homologs/`](4_homologs/) | P5 §1.3–1.4; task 2.5 | BLAST restricted to each of ten species (RefSeq protein, RefSeq RNA): rat, mouse, dog, opossum, platypus, turtle, frog, coelacanth, zebrafish, chicken | Results (species), Table 2 |
| [`5_conservation/`](5_conservation/) | P5 §2.3–2.4; tasks 3.2, 3.9 | Clustal Omega on the EMBL-EBI website: the ten sequences, screenshots, alignment, percent identity matrix | Table 2, Figure 3 |
| [`6_clinical/`](6_clinical/) | P1; task 2.3 | gnomAD search path and variant page; UniProt P30613 variant annotation; ClinVar record summary (screenshots by hand: see `procedure/ClinVar/`) | Results, Table 1 |
| [`7_structure/`](7_structure/) | P4, P6; tasks 2.6, 3.3, 3.4 | RCSB PDB pages and files of 2VGB and 2VGI; the PyMOL command list (`structure_commands.pml`), a screenshot of the PyMOL window after each step, PyMOL's output log, and Figures 1 and 2 | Results (structure), Figures 1 and 2 |
| [`8_lab1_three_species/`](8_lab1_three_species/) | P5 §2.3–2.4 (Lab 1) | BLAST alignment of human with rat and mouse, CDS and protein | Lab 1 notebook |
| [`9_dog/`](9_dog/) | task 3.2 (outgroup) | TimeTree divergence times, dog CDS, PubMed abstracts on canine PK deficiency | Results and Discussion (dog) |

**Check values** (what each step must give): wild type vs mutant 1724/1725 nucleotides and 573/574 amino acids identical, no gaps; NM_000298.6 and NP_000289.1 at 100% identity; Clustal Omega identities with human from 67.8% (zebrafish) to 92.5% (dog); position 486 Arg in nine species, Gln in rat; PyMOL `RMSD = 0.587 (517 to 517 atoms)`, distances 26.1, 2.7, 3.0, and 3.4 Å.

BLAST results stay on the NCBI server for 36 hours; the request ID (RID) of each search is at the top of its text report.
