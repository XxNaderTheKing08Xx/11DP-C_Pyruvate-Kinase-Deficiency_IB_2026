# Pyruvate Kinase Deficiency: from a *PKLR* variant to haemolytic anaemia


Aarón Godino Martínez, Bruno Carrasquilla, Ettore Colombo, Leo Nader Al Hourani, Pol Monllau

Group 11DP-C, BSc in Bioinformatics (UPC-UB-UAB-UPF), Barcelona, Spain



## Abstract

## Keywords

## Background

Pyruvate kinase (PK) catalyses the last, irreversible step of glycolysis, from phosphoenolpyruvate (PEP) and adenosine diphosphate (ADP) to pyruvate and adenosine triphosphate (ATP) (Bianchi and Fermo 2020). *PKLR* (1q22; NCBI 2026b), a paralogue of *PKM*, encodes the R- and L-type isozymes from two tissue-specific promoters whose mRNAs differ in their first exon (Bianchi and Fermo 2020), as in rat (Noguchi *et al.* 1987); the first 33 of the 574 R-type residues are replaced by two in the 543-residue L-type, whose numbering is therefore 31 lower (UniProt Consortium 2026). R-type PK is a homotetramer activated allosterically by fructose 1,6-bisphosphate (FBP); each subunit has N, A, B, and C domains, and FBP binds in the C domain (Valentini *et al.* 2002).

Red cells lack mitochondria and depend on glycolysis for ATP (Fattizzo *et al.* 2022, Luke *et al.* 2023); in pyruvate kinase deficiency (PKD), the most common glycolytic enzyme defect associated with congenital non-spherocytic haemolytic anaemia, insufficient red-cell PK activity lowers ATP and leads to premature removal of red cells, mainly in the spleen (Bianchi and Fermo 2020, Fattizzo *et al.* 2022). Inheritance is autosomal recessive; 142 of 177 unrelated patients of the PKD Natural History Study were compound heterozygous (Bianchi *et al.* 2020). Diagnosed PKD has a prevalence of 3.2–8.5 per million in Western populations (Secrest *et al.* 2020); among 254 patients, 48% had iron overload and 45% gallstones (Grace *et al.* 2018). Besides transfusion, splenectomy, and iron chelation (Grace and Barcellini 2020), the PK activator mitapivat, approved for adults (Luke *et al.* 2023), produced a haemoglobin response in 16 of 40 patients and in none of 40 on placebo in a phase 3 trial (Al-Samkari *et al.* 2022).

*PKLR* c.1456C>T, p.(Arg486Trp) (R486W), described in 1993 (Baronciani and Beutler 1993), was the most frequent *PKLR* variant in Italian and Spanish series (Zanella *et al.* 1997, Zarza *et al.* 1998). Its crystal structure is essentially unchanged although its catalytic efficiency is reduced (Valentini *et al.* 2002), and the conservation of position 486 across species had not been reported in the sources we read. To answer the course question of how sequence and structure changes explain disease and what comparison shows about function (IB Teaching Staff 2026), we aimed to (i) identify the cause of the disease at gene, protein, structure, and function level; (ii) describe how the mutant differs from the wild type; (iii) explain the consequences for the patient; and (iv) compare human PKLR with its rat, mouse, and dog orthologues, focusing on position 486.

## Results

### The supplied sequences differ by c.1456C>T, p.(Arg486Trp)

### Position 486 is Arg in most orthologues and Gln in rat

### Trp486 causes no change beyond model-to-model differences

### Clinical, population, and functional evidence

## Discussion

### From the variant to the patient

### Comparison with Valentini *et al.* (2002)

### What the species comparison shows about the affected function

### Population frequency and penetrance

### Limitations

## Conclusions

## Methods

### Data sources and species

Records (versions in the tables and reference list) were accessed on 22 October 2026 from the HUGO Gene Nomenclature Committee (HGNC:9020), National Center for Biotechnology Information (NCBI) Gene (NCBI 2026b), RefSeq (O'Leary *et al.* 2016) with MANE Select status, Ensembl 116 and Compara (Dyer *et al.* 2025), UniProtKB/Swiss-Prot (UniProt Consortium 2026), the PDB (RCSB PDB 2026), ClinVar (NCBI 2026a), gnomAD v4.1 and v2.1.1 (gnomAD 2026), and Online Mendelian Inheritance in Man (OMIM 2026). As the rat RefSeq Select protein (NP_036756.3) is the L-type, the rat R-type CDS was taken from XM_006232593.5, whose translation equals UniProt P12928. Rat is required by our assignment, mouse is its closest relative, and dog, the fourth reviewed entry, is the outgroup (divergence about 89 million years ago, against 84 between primates and rodents; TimeTree 5, Kumar *et al.* 2022).

### Analysis

Sequences were compared with Biopython 1.88 and named after the Human Genome Variation Society recommendations (2016 update). Pairwise global alignments used the EMBOSS Needle defaults (BLOSUM62, gap opening 10, extension 0.5) and were reproduced with EMBOSS needle 6.6.0 and BLAST+ 2.17.0; multiple alignments used Clustal Omega 1.2.4, run locally with default parameters, and ancestral states were inferred by inspection under maximum parsimony. Contacts (heavy atoms ≤ 4.0 Å; polar N/O ≤ 3.5 Å) were computed with Biopython, relative solvent accessibility with FreeSASA 2.2.1, and secondary structure with DSSP in UCSF ChimeraX 1.10.1 and with PyMOL 3.1.8 (Schrödinger), which rendered Figures 2 and 3. Cα atoms were superposed by least squares without outlier rejection, and rat residues were mapped as L-type number + 31. Expected homozygotes were computed as (AN/2) × AF² and carriers as 2pq. Parameters, scripts, and step-by-step protocols are in the README files of the `analysis/` folders.

## Abbreviations

## Declarations

### Availability of data and materials

### Authors' contributions

### Use of generative AI

## References
