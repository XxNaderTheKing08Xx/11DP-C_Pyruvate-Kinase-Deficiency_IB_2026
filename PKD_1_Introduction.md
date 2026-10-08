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

We are **Group 11DP-C**: Aarón Godino Martínez, Bruno Carrasquilla, Ettore Colombo, Leo Nader Al Hourani, and Pol Monllau, first-year students of the BSc in Bioinformatics (UPC-UB-UAB-UPF). In the Disease Project of *Introduction to Bioinformatics*, each group investigates the molecular basis of an inherited disease by analysing the sequences of the gene responsible and of the protein it encodes, and the three-dimensional structure of that protein. Our assigned disease is **Pyruvate Kinase Deficiency (PKD)**. This document is our project plan: objectives, work, responsibilities, schedule, and risks.

### Pyruvate Kinase Deficiency (PKD) in brief

Pyruvate kinase (PK) catalyses the last step of glycolysis: the transfer of a phosphate group from phosphoenolpyruvate to ADP, which yields pyruvate and ATP. Mature red blood cells have no mitochondria and depend entirely on glycolysis for ATP; when PK activity is insufficient, they are removed prematurely from the circulation, mainly in the spleen (van Wijk and van Solinge 2005, Bianchi and Fermo 2020). The result is a chronic non-spherocytic haemolytic anaemia, ranging from fully compensated haemolysis to life-threatening anaemia in the newborn; the most frequent complications are iron overload and gallstones (Grace *et al.* 2018).

PKD is caused by pathogenic variants in *PKLR*, which encodes two isoforms from tissue-specific promoters: the **R-type** of red blood cells (574 amino acids) and the **L-type** of the liver (543 amino acids) (Kanno *et al.* 1992, UniProt Consortium 2026). Inheritance is autosomal recessive. More than 300 different variants have been reported, and most patients carry two different ones (compound heterozygous) (Bianchi *et al.* 2020). The prevalence of diagnosed PKD is 3.2–8.5 cases per million in Western populations (Secrest *et al.* 2020). Treatment has long been supportive (transfusions, splenectomy, and iron chelation); the oral PK activator mitapivat is now approved for adults, and gene therapy is in clinical trials (Luke *et al.* 2023).

### The data we were given

Our instructors provided four FASTA files:

| File | Content | Length |
|---|---|---|
| `PKD_CDS_WT.fasta` | coding sequence, wild type | 1 725 nucleotides |
| `PKD_CDS_mutant.fasta` | coding sequence, mutant | 1 725 nucleotides |
| `PKD_protein_WT.fasta` | protein sequence, wild type | 574 amino acids |
| `PKD_protein_mutant.fasta` | protein sequence, mutant | 574 amino acids |

The headers identify neither the gene, the isoform, nor the variant. Comparing wild-type and mutant sequences will reveal the variant, which we will name at the DNA (*c.*) and protein (*p.*) levels in HGVS nomenclature (den Dunnen *et al.* 2016).

### Our question and objectives

The project asks: *how can changes in gene/protein sequence and structure explain the molecular basis of disease, and what can comparison tell us about the affected biological function?* Our general objective is to **explain** how the variant in our mutant sequences alters PK and causes PKD, and what comparison with other species reveals about the affected function. Our four specific objectives follow the topics of the Disease Project Overview:

1. **Identify** the cause of the disease from the gene, protein, structure, and function perspectives.
2. **Describe** the differences between the wild-type and the mutant gene, protein, structure, and function.
3. **Explain** how these differences affect the patient.
4. **Compare** human PK with its orthologs in at least two other species, including the rat.

### Our approach

| Step | Work | Class session | Objective |
|---|---|---|---|
| 1. Sequence characterization | Translate both CDSs; identify and name the variant | Data sources | 1, 2 |
| 2. Background review | Gene and enzyme function; mechanism and clinical features of PKD | Data sources; citation | 1, 3 |
| 3. Comparative analysis | Align orthologous sequences; assess conservation of the substituted position | Sequence alignment (8, 15 Oct) | 4 |
| 4. Structural analysis | Locate the substituted residue; superpose wild-type and mutant structures | PDB (1 Oct); superposition (22 Oct) | 1, 2 |
| 5. Pathogenicity assessment | Combine clinical, population, conservation, and structural evidence | — | 2, 3 |
| 6. Synthesis and report | Write the manuscript in *BMC Bioinformatics* format | Project management and GitHub | all |

The steps are grouped into four work packages: **WP1** project management and GitHub; **WP2** data sources, literature, and sequences (steps 1–2); **WP3** sequence and structure analysis (steps 3–4); **WP4** report, delivery, and defence (steps 5–6). Tasks are distributed equally among the five members; each task has a single owner, and each member commits their own work.

### Species compared

The rat (*Rattus norvegicus*) is required by our disease assignment; it is also the species in which the production of both isoforms from a single gene was demonstrated (Noguchi *et al.* 1987). As the second species we propose the mouse (*Mus musculus*), a close relative of the rat in which spontaneous *Pklr* mutations cause a haemolytic anaemia similar to the human disease (Morimoto *et al.* 1995); the comparison shows whether a difference between the rat and human sequences is specific to the rat lineage. The dog (*Canis lupus familiaris*), whose lineage separated from that of primates and rodents before these two diverged (Kumar *et al.* 2022), will serve as an outgroup to establish which residue is ancestral.

### Deliverables and schedule

The project is developed in a GitHub repository, which must contain the final **report** (Markdown, *BMC Bioinformatics* format, at most 7 pages or 5 000 words and five figures or tables), a live **Project Status document**, this **Gantt chart** (frozen on 9 October, plus a live version), and all **references** (IB Teaching Staff 2026a, 2026b). The project runs from 17 September 2026, when it was presented and our repository was created, to the oral defence on 5 November, Monday to Friday, excluding the holidays of 24 September and 12 October:

- **9 October:** submission of this plan;
- **22 October:** discussion of the Project Status document;
- **28 October, 18:00:** delivery of the final report;
- **5 November:** oral defence.

### Risks and structure of this document

Our most serious risks, derived from the course requirements, are bulk or late commits, unequal distribution of work, incorrect database records, non-compliance with the report format, and missed deadlines. Section 1 details the objectives, milestones, and team; Section 2 the Gantt chart and task specifications; Section 3 the risk analysis with contingency plans.

---
[Index](README.md) · [1. Project overview and team →](PKD_2_Overview_and_Team.md)
