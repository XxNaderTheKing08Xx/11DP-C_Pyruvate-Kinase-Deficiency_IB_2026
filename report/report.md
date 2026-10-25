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

The supplied coding sequences (CDSs) are 1725 nucleotides long, have no internal stop codon, and translate into the supplied 574-residue proteins. The wild-type CDS is identical to that of NM_000298.6 (Matched Annotation from NCBI and EMBL-EBI [MANE] Select), and the wild-type protein to NP_000289.1 and UniProt P30613-1. The mutant differs only at c.1456C>T, the first position of codon 486 (CGG→TGG), a CpG site, giving p.(Arg486Trp). c.1456 lies in exon 10 of NM_000298.6 (exon 11 in the 12-exon gene numbering; Bianchi *et al.* 2020); on the L-type transcript (NM_181871.4) the change is c.1363C>T, p.(Arg455Trp).

**Table 1.** Evidence for *PKLR* c.1456C>T and for the other substitution recorded at codon 486 (gnomAD 2026, NCBI 2026a, OMIM 2026, UniProt Consortium 2026); identifiers, coordinates, and in silico scores are in `analysis/evidence/table1_variant_evidence.md`. gnomAD, Genome Aggregation Database; AC, allele count; AN, allele number; AF, allele frequency; k_cat, turnover number; S₀.₅, substrate concentration giving half-maximal rate; WT, wild type.

| Field | c.1456C>T, p.(Arg486Trp) | c.1457G>A, p.(Arg486Gln) |
|---|---|---|
| ClinVar (criteria provided, multiple submitters, no conflicts) | VCV000001513.102: Pathogenic/Likely pathogenic; 33 submitters | VCV000873616.7: Uncertain significance; two submitters |
| gnomAD v4.1, exomes + genomes | AC 4766 / AN 1 613 440; AF 0.295%; eight homozygotes | AC 47 / AN 1 613 264; AF 0.0029%; no homozygotes |
| Highest AF by group | Finnish 0.78% (v4.1); southern European 0.52% (v2.1.1 exomes, 60/11 494) | South Asian 0.012% (v4.1) |
| Share of PKD alleles | 31 of 354 unrelated alleles (9%), second most frequent (Bianchi *et al.* 2020); seven of 30 in Italy (Zanella *et al.* 1997); seven of 22 in Spain (Zarza *et al.* 1998) | not reported in the sources read |
| Recombinant enzyme, R486W vs WT (Valentini *et al.* 2002) | k_cat 195 vs 355 s⁻¹ (ratio 0.55); S₀.₅(PEP) 1.69 vs 1.10 mM; k_cat/S₀.₅(PEP) 116 vs 323 s⁻¹ mM⁻¹ (0.36) without and 492 vs 1972 (0.25) with 1 mM FBP; k_cat/S₀.₅(ADP) ratio 0.44; thermostability at 53 °C not reduced | no data |

### Position 486 is Arg in most orthologues and Gln in rat

The only reviewed UniProt PKLR entries (Table 2) are 574-residue R-type sequences and Ensembl Compara one-to-one orthologues, each more identical to human PKLR than to human PKM (rat 91.1% vs 64.5%, mouse 91.6% vs 64.1%, dog 92.5% vs 65.3%). They align without gaps; rat differs from human at 51 residues and at 152 further, synonymous codons (mouse 48 and 155; dog 43 and 121).

In the window 476–496, 19 of 21 columns are identical in the four species; the variable columns are 476 (Thr/Lys) and 486 (Arg; Gln in rat) (Figure 1). In 21 vertebrate sequences, column 486 holds Arg in 18 and Gln in three (rat, macaque, and marmoset), and Ser485 and Tyr487–Arg490 are invariant. Among the 179 jawed-vertebrate sequences of all 190 Ensembl Compara homologues, position 486 is Arg in 156, Gln in 13 (rat, six Old World and four New World monkeys, hedgehog, and sperm whale), Tyr in one (a misaligned gene model), and unresolved in nine; Trp occurs in none of the 190. By maximum parsimony, Arg is ancestral, and Gln arose repeatedly, including on the rat lineage and twice in primates.

![Figure 1](../figures/Figure1_PKLR_MSA_486.png)

**Figure 1.** Conservation of human PKLR position 486. Clustal Omega alignment (476–496) of the Swiss-Prot R-type sequences (top) and 17 Ensembl 116 orthologues (bottom), each in its own numbering. Framed column, human Arg486; colours and Clustal symbols as in the key; bars, fraction of sequences with the majority residue.

**Table 2.** Human PKLR (R-type) compared with its orthologues. Protein: global (Needleman–Wunsch) identity and similarity as fractions of alignment columns. CDS: position-by-position identity of the gap-free 1725-nucleotide CDSs. Mouse NP_038659.2 differs from P53657 at residues 279 and 335.

| Species | Protein (UniProtKB / RefSeq) | CDS (RefSeq) | Protein identity; similarity | CDS identity | Residue at human 486 (codon) |
|---|---|---|---|---|---|
| Human | P30613 / NP_000289.1 | NM_000298.6 | — | — | Arg486 (CGG) |
| Rat | P12928 / XP_006232655.1 | XM_006232593.5 | 523/574 (91.1%); 95.1% | 1491/1725 (86.4%) | Gln486 (CAA) |
| Mouse | P53657 / NP_038659.2 | NM_013631.2 | 526/574 (91.6%); 95.5% | 1490/1725 (86.4%) | Arg486 (CGC) |
| Dog | Q29536 / NP_001243191.1 | NM_001256262.2 | 531/574 (92.5%); 95.6% | 1544/1725 (89.5%) | Arg486 (CGG) |
| 17 vertebrates (Figure 1) | Ensembl 116 orthologues | — | 62.6–99.7% identity | — | Arg in 15; Gln in macaque, marmoset |

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
