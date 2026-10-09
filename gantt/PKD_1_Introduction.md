# Pyruvate Kinase Deficiency: from a *PKLR* variant to haemolytic anaemia

**Project plan: Gantt chart and risk analysis**

Aarón Godino Martínez, Bruno Carrasquilla, Ettore Colombo, Leo Nader Al Hourani, Pol Monllau  
Group 11DP-C · BSc in Bioinformatics (UPC-UB-UAB-UPF), Barcelona, Spain  
*Introduction to Bioinformatics* (2026/27-01, FIB-2703103) · 9 October 2026

---

**Contents:** [1. Project overview and team](PKD_2_Overview_and_Team.md) · [2. Gantt chart](PKD_3_Gantt_Chart.md) · [2. Gantt chart: task details](PKD_4_Task_Details.md) · [3. Risk analysis](PKD_5_Risk_Analysis.md) · [References](PKD_6_References.md)

---

## Introduction

### Who we are

We are **Group 11DP-C**: Aarón Godino Martínez, Bruno Carrasquilla, Ettore Colombo, Leo Nader Al Hourani, and Pol Monllau, five first-year students of the BSc in Bioinformatics (UPC-UB-UAB-UPF). In the Disease Project of *Introduction to Bioinformatics*, each group is assigned a disease and investigates its molecular basis by analysing the sequences of the gene responsible and of the protein it encodes, and the three-dimensional structure of that protein. Our assigned disease is **Pyruvate Kinase Deficiency (PKD)**.

This document is our project plan: it defines our objectives, the work needed to reach them, who is responsible for each task, the schedule, and the risks we foresee.

### Pyruvate Kinase Deficiency (PKD) in brief

**The enzyme.** Pyruvate kinase (PK) catalyses the last step of glycolysis: the irreversible transfer of a phosphate group from phosphoenolpyruvate to adenosine diphosphate (ADP), which yields pyruvate and adenosine triphosphate (ATP). Its activity is regulated allosterically (Bianchi and Fermo 2020). The enzyme is a tetramer of four identical subunits, and both of its human isoforms are activated by fructose 1,6-bisphosphate, an intermediate of an earlier step of glycolysis (Valentini *et al.* 2002).

**Why red blood cells are affected.** Mature red blood cells lack a nucleus and organelles, including mitochondria, and therefore depend entirely on glycolysis for the ATP they need to maintain their integrity and function (Bianchi and Fermo 2020, Fattizzo *et al.* 2022, Luke *et al.* 2023). When PK activity is insufficient, ATP levels fall and the upstream metabolite 2,3-diphosphoglycerate (2,3-DPG) accumulates; the red cells become rigid and less flexible and are removed prematurely from the circulation, mainly in the spleen (Bianchi and Fermo 2020, Fattizzo *et al.* 2022). Because the anaemia is chronic and 2,3-DPG favours the release of oxygen to the tissues, patients may tolerate it better than other forms of anaemia (Fattizzo *et al.* 2022). The result is a chronic non-spherocytic haemolytic anaemia, whose severity ranges from mild or fully compensated haemolysis to life-threatening neonatal anaemia that requires exchange transfusion (Bianchi and Fermo 2020). In the 254 patients of the PKD Natural History Study, the most frequent complications were iron overload (48%) and gallstones (45%); pulmonary hypertension, extramedullary haematopoiesis, and leg ulcers were also observed (Grace *et al.* 2018).

**The gene.** PKD is caused by pathogenic variants in *PKLR*, located at 1q22 (HGNC 2026). The gene comprises 12 exons. Two tissue-specific promoters produce two mRNAs that differ only in their first exon: exon 1 is specific to the mRNA of the **R-type** isoform, present in red blood cells, and exon 2 to that of the **L-type** isoform, expressed in the liver; the other 10 exons are shared (Kanno *et al.* 1992, Bianchi and Fermo 2020). The R-type isoform has 574 amino acids and the L-type 543. The two differ only in their N-terminal region, where the first 33 residues of the R-type are replaced by two residues in the L-type; the same residue therefore has a number 31 positions lower in the L-type than in the R-type (UniProt Consortium 2026). Inheritance is autosomal recessive (Bianchi and Fermo 2020): affected individuals carry a pathogenic variant on both alleles of *PKLR*, either the same variant (homozygous) or two different ones (compound heterozygous). In the PKD Natural History Study, 142 of 177 unrelated patients (80%) were compound heterozygous, and 127 different pathogenic variants were identified, 84 of them missense (Bianchi *et al.* 2020); more than 300 different variants have been reported overall (Bianchi and Fermo 2020).

**Frequency and treatment.** PKD is the most common glycolytic enzyme defect associated with congenital non-spherocytic haemolytic anaemia (Bianchi and Fermo 2020). A systematic review estimated the prevalence of clinically diagnosed PKD in Western populations at 3.2–8.5 cases per million, and that of diagnosed and undiagnosed cases together at up to 51 per million, a figure extrapolated from allele frequencies (Secrest *et al.* 2020). Treatment has long been supportive: transfusions, splenectomy, and iron chelation (Grace and Barcellini 2020). Two disease-modifying treatments have since emerged: mitapivat, an oral allosteric activator of PK now approved for adults with PKD, and lentiviral gene therapy, which is being evaluated in clinical trials (Luke *et al.* 2023). In a phase 3 trial in adults not receiving regular transfusions, mitapivat increased haemoglobin levels compared with placebo (Al-Samkari *et al.* 2022).

### The data we were given

Our instructors provided four files in FASTA format:

| File                       | Content                     | Length           |
| -------------------------- | --------------------------- | ---------------- |
| `PKD_CDS_WT.fasta`         | coding sequence, wild type  | 1725 nucleotides |
| `PKD_CDS_mutant.fasta`     | coding sequence, mutant     | 1725 nucleotides |
| `PKD_protein_WT.fasta`     | protein sequence, wild type | 574 amino acids  |
| `PKD_protein_mutant.fasta` | protein sequence, mutant    | 574 amino acids  |

The headers (for example, `>PKD_CDS_WT`) state only the disease, the type of sequence, and whether the sequence is wild type or mutant. They contain no identifier for the gene, the isoform, or the variant; establishing these is our first task.

### Key concepts

**FASTA format.** A plain-text format for biological sequences, named after the FASTA sequence-comparison programs (Pearson and Lipman 1988): a header line beginning with `>`, followed by the sequence written with one letter per nucleotide or amino acid. It is the standard input format of most bioinformatics tools.

**Coding sequence and protein sequence.** The coding sequence (CDS) is the part of a gene's transcript, from the start codon to the stop codon, that is translated into protein; it is usually written as a DNA sequence. Translating it codon by codon (one codon is three consecutive nucleotides) gives the amino acid sequence of the protein. Having both sequences allows us to describe a variant twice: as a change in the coding DNA and as its consequence for the amino acid sequence. The nomenclature of the Human Genome Variation Society (HGVS) marks the first with the prefix *c.* and the second with *p.* (den Dunnen *et al.* 2016).

**Wild type and mutant.** The wild-type sequence is the reference sequence; the mutant sequence carries a disease-associated variant. Comparing the two establishes the position and the nature of the variant.

**Protein structure.** Experimentally determined three-dimensional structures of proteins are deposited in the Protein Data Bank (PDB) (Berman *et al.* 2000). For human R-type PK, the structures of the wild-type enzyme and of disease-associated mutants have been determined, which makes it possible to relate each amino acid substitution to its position relative to the active site, the allosteric site, and the interfaces between domains and between subunits (Valentini *et al.* 2002).

**Orthologs and conservation.** Orthologs are genes in different species that descend from a single gene of their last common ancestor and diverged through speciation (Koonin 2005). Aligning the sequences of orthologous proteins shows which positions have been conserved during evolution. A residue conserved over long evolutionary distances is likely to be under structural or functional constraint, and sequence conservation is one of the most widely used predictors of functionally important residues (Capra and Singh 2007).

### Our question and objectives

The Disease Project asks every group to answer the same question: *how can changes in gene/protein sequence and structure explain the molecular basis of disease, and what can comparison tell us about the affected biological function?*

Our general objective is to **explain** how the variant present in our mutant sequences changes the *PKLR* coding sequence, the sequence, structure, and function of the PK protein, and how these changes cause PKD; and to establish what comparing human PK with its orthologs in other species reveals about the affected function. We have set four specific objectives, one for each topic of the Disease Project Overview (IB Teaching Staff 2026a):

1. **Identify** the cause of the disease from the gene, protein, structure, and function perspectives.
2. **Describe** the differences between the wild-type and the mutant gene, protein, structure, and function.
3. **Explain** how these differences affect the patient.
4. **Compare** the gene, protein, structure, and function of human PK with those of its orthologs in at least two other species, including the rat, and infer what the comparison reveals about the affected function.

### Our approach

The work is organized in six steps, each with a defined output; most of them apply methods taught in a specific class session.

| Step                         | Work                                                         | Main sources and tools                       | Class session                                 | Objective |
| ---------------------------- | ------------------------------------------------------------ | -------------------------------------------- | --------------------------------------------- | --------- |
| 1. Sequence characterization | Verify that each CDS translates into the corresponding protein sequence; identify and name every difference between wild type and mutant, in the coding sequence and in the protein sequence | Sequence translation and pairwise comparison | Data sources                                  | 1, 2      |
| 2. Background review         | Review the normal function of the gene and the enzyme, and the mechanism, inheritance, and clinical features of PKD | PubMed, NCBI Gene, OMIM, UniProt             | Data sources; paper structure and citation    | 1, 3      |
| 3. Comparative analysis      | Retrieve the orthologous sequences, align them with the human sequence, and assess the conservation of the substituted position | UniProt, Ensembl, Clustal Omega              | Sequence alignment (8 and 15 Oct)             | 4         |
| 4. Structural analysis       | Locate the substituted residue in experimental structures, relative to the active and allosteric sites; superpose the wild-type and mutant structures | PDB, molecular visualization software        | PDB (1 Oct); structure superposition (22 Oct) | 1, 2      |
| 5. Pathogenicity assessment  | Combine clinical classification, population frequency, conservation, structural, and biochemical data, distinguishing measured results from hypotheses and open questions | ClinVar, gnomAD, primary literature          | Integrated bioinformatics reasoning (19 Oct)  | 2, 3      |
| 6. Synthesis and report      | Write the answer to the project question as a manuscript in the *BMC Bioinformatics* format | GitHub, Markdown                             | Project management and GitHub                 | all       |

These steps are grouped into four work packages: **WP1**, project management and GitHub, for the whole duration of the project; **WP2**, data sources, literature, and sequences (steps 1 and 2); **WP3**, sequence and structure analysis (steps 3 and 4); and **WP4**, report, delivery, and defence (steps 5 and 6, and the preparation of the oral defence). The tasks are distributed among the five members so that each has a similar workload (10 to 12 working days of planned effort; Section 1), each task has a single owner, and each member commits their own work, so that the history of the repository records individual contributions.

### Species compared

The comparison must include the rat (*Rattus norvegicus*), as specified in our disease assignment. The rat is also the species in which the production of both PK isoforms from a single gene — the L-type, expressed in the liver, and the R-type, present in red blood cells — was demonstrated (Noguchi *et al.* 1987), and experimental structures of its L-type isoform are available (Gassaway *et al.* 2019).

As the second species we propose the mouse (*Mus musculus*). Like those of human and rat, its PK sequence is a manually reviewed UniProt entry (UniProt Consortium 2025), and mice of the CBA/N strain carrying a spontaneous mutation at the *Pk-1* locus (the former symbol of *Pklr*) develop a non-spherocytic haemolytic anaemia, as in the human disease (Morimoto *et al.* 1995). Because rat and mouse are closely related rodents, comparing them shows whether a difference between the rat and human sequences is specific to the rat lineage. To infer which residue is ancestral, we will include the dog (*Canis lupus familiaris*) in the alignment as an outgroup: the dog lineage separated from the common ancestor of primates and rodents about 89 million years ago, before primates and rodents separated from each other, about 84 million years ago (Kumar *et al.* 2022).

### Deliverables and schedule

The project is developed in a **GitHub repository**, to which each member commits their own contributions as the work progresses. At the end, the repository must contain the final **report** (a Markdown manuscript in the *BMC Bioinformatics* format, of at most 7 pages or 5 000 words, whichever is shorter, with at most five figures or tables), a live **Project Status document**, this **Gantt chart** (the version frozen on 9 October and a live version updated until the end of the project), and all the **references** used (IB Teaching Staff 2026b, 2026c).

The project runs from 17 September 2026, when it was presented in class and our repository was created, to the oral defence on 5 November, Monday to Friday, excluding the holidays of 24 September (local) and 12 October (national). The key dates are:

- **9 October:** submission of this plan (18:00 in the Disease Project Overview, 22:00 on the course's Moodle page; we work to the earlier time);
- **22 October:** discussion of our Project Status document in class;
- **28 October, 18:00:** delivery of the final report on GitHub (the Overview gives "Monday October 28th", but 28 October 2026 is a Wednesday, so we have asked the instructors to confirm the date);
- **5 November:** oral defence.

### Risks

We identified our risks from the requirements stated in the course instructions. The five scored as high are committing work in bulk or too late, an unequal contribution among members, selecting an incorrect database record or misinterpreting a result, developing the parts of the project in isolation, and failing to meet the report format or the delivery deadline. Each risk has an owner and a contingency plan.

### Structure of this document

Section 1 details the objectives, the milestones, and the team; Section 2 presents the Gantt chart and the task specifications; Section 3 presents the risk analysis. References follow the Mini Oxford style required by the *BMC Bioinformatics* template.

---

[Index](README.md) · [1. Project overview and team →](PKD_2_Overview_and_Team.md)
