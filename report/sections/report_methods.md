## Methods

### Data sources

The wild-type and mutant coding and protein sequences were provided by the course (IB Teaching Staff 2026). We identified them with BLAST at NCBI (Altschul *et al.* 1990), using blastn for the coding sequences and blastp for the proteins against the default databases, and chose the best curated RefSeq hits (O'Leary *et al.* 2016): the *PKLR* transcript NM_000298.6 (NCBI 2026b) and its protein NP_000289.1, which is identical to UniProt P30613 (UniProt Consortium 2026). Variant information came from ClinVar (record VCV000001513; NCBI 2026a), gnomAD v4.1.2 (gnomAD 2026), and UniProt. The reviewed PKLR entries of rat (P12928), mouse (P53657), and dog (Q29536) came from UniProt (UniProt Consortium 2026); we checked them against BLAST searches restricted to each species (databases RefSeq protein and RefSeq RNA), in which we chose the PKLR entry covering the whole human query, that is, the red-cell form. For other vertebrates, we used the one-to-one orthologues of human *PKLR* (ENSG00000143627) listed by Ensembl (Dyer *et al.* 2025) and the corresponding UniProt entries where available (UniProt Consortium 2025). Divergence times between species were taken from TimeTree (Kumar *et al.* 2022). The crystal structures of the wild-type enzyme (PDB 2VGB) and of the R486W mutant (PDB 2VGI) (Valentini *et al.* 2002) were downloaded from the RCSB Protein Data Bank (RCSB PDB 2026). Copies of the database records are stored in the `references/` folder of our repository, and every analysis step (input files, screenshots, and output files) in the `analysis/procedures/` folder.

### Sequence analysis

Pairwise alignments were made with BLAST on the NCBI website, with the option "Align two or more sequences" (blastn for coding sequences, blastp for proteins) and the default settings; the variant was located in the "Query-anchored with dots for identities" view. The multiple alignment of the ten PKLR sequences was made with Clustal Omega (Sievers *et al.* 2011) on the EMBL-EBI website (Madeira *et al.* 2024), with the default settings; percent identities were taken from its percent identity matrix, and position 486 was read in the column aligned with human Arg486.

### Structure analysis

Structures were analysed with PyMOL (Schrödinger, LLC). For each structure we used chain A. Domains were coloured with the residue limits given by Valentini *et al.* (2002): N, 57–84; A, 85–159 and 263–431; B, 160–262; C, 432–574. The mutant was superposed on the wild type with the PyMOL command `align` on the Cα atoms of chain A, without outlier rejection (`cycles=0`), which gives the RMSD. Distances and hydrogen bonds were measured with the PyMOL commands `distance` and `get_distance`. Figures 1 and 2 were made with the same list of PyMOL commands, which is stored with the analysis.

### Statistical analysis

No statistical analysis was performed.
