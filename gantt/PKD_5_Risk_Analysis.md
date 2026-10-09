## 3. Risk analysis

We identified our risks from the requirements of the course documents (Disease Project Overview, Formal Requirements, *BMC Bioinformatics* template, Use of GenAI in IB) and from the descriptors of the marking rubric, and we checked them against the current state of our repository on 8 October. Following the method taught in class, each risk is scored as **Risk = Likelihood × Impact**, both from 1 (low) to 3 (high): **1–2 Low**, **3–5 Medium**, **6–9 High**. Each risk has an owner, who monitors it and activates the contingency plan when needed.

### Risk matrix

| Likelihood ↓ / Impact → | 1 (low) | 2 (moderate) | 3 (high) |
|---|---|---|---|
| **3 (likely)** | — | — | R1, R2 (9) |
| **2 (possible)** | — | R4, R6, R9 (4) | R5, R8, R10 (6) |
| **1 (unlikely)** | — | R7 (2) | R3 (3) |

### Risk register

| ID | Risk | L | I | Score | Level | Owner | Status on 8 Oct | Source |
|---|---|---|---|---|---|---|---|---|
| R1 | Commits made in bulk or near the deadlines | 3 | 3 | **9** | High | Ettore | Partly occurred: 11 commits, none since 28 Sep | Formal Requirements I; rubric, Project Documentation |
| R2 | Unequal contribution among members | 3 | 3 | **9** | High | Pol | Partly occurred: three of five members have committed | Formal Requirements I; rubric, Project Documentation and Project Execution; Use of GenAI in IB |
| R3 | Repository structure not compliant with the requirements | 1 | 3 | **3** | Medium | Bruno | Occurred: required folders missing; correction planned in task 1.1 | Formal Requirements I–II; rubric, Project Documentation |
| R4 | Project Status document or live Gantt chart not kept up to date | 2 | 2 | **4** | Medium | Ettore | Open | Formal Requirements III–IV; rubric, Project Scheduling |
| R5 | Incorrect database record selected, or database result misinterpreted | 2 | 3 | **6** | High | Bruno | Open | Rubric, Bioinformatic database literacy |
| R6 | Alignment or structural results misinterpreted | 2 | 2 | **4** | Medium | Pol | Open | Rubric, Bioinformatic analysis; Structural and evolutionary reasoning |
| R7 | Species comparison insufficiently justified, or evolutionary changes inferred without an outgroup | 1 | 2 | **2** | Low | Leo | Open | Overview, Scope Expectations; rubric, Bioinformatic analysis |
| R8 | Parts of the project developed in isolation; members aware only of their own tasks | 2 | 3 | **6** | High | Leo | Open | Rubric, Knowledge integration; Oral scientific communication |
| R9 | Statements supported by weak or unverified sources, including unverified GenAI output | 2 | 2 | **4** | Medium | Aarón | Open | Rubric, Evidence gathering; Use of GenAI in IB |
| R10 | Report format, length limits, or delivery deadline not met | 2 | 3 | **6** | High | Aarón | Open | Overview, Timelines; Formal Requirements I; *BMC Bioinformatics* template |

### Analysis and contingency plans

**High risks**

- **R1 (9).** Our identified risk is *committing our work in bulk or near the deadlines*, we have realized because the Formal Requirements state that a repository committed once or very few times before the deadline will be heavily penalized, and our history shows 11 commits made on two days (17 and 28 September) and none since; our contingency plan consists of committing the plan section by section on 9 October (tasks 1.1–1.5), committing at the end of every work session thereafter, and reviewing the commit history at each weekly meeting (task 1.7).
- **R2 (9).** Our identified risk is *an unequal contribution among members*, we have realized because the Formal Requirements state that contributions concentrated in one member severely penalize the whole team, only human commits are assessed (Use of GenAI in IB), and two of us have not yet committed; our contingency plan consists of an equal distribution of the planned work (seven tasks per member, Section 1), each member committing their own work from 9 October, recording each member's contributions in the Project Status document, and an Authors' contributions section in the report.
- **R5 (6).** Our identified risk is *selecting an incorrect database record or misinterpreting a database result*, we have realized because several records exist for the same object: *PKLR* encodes two isoforms, UniProt contains reviewed and unreviewed entries, the gene symbol (*PKLR*) differs from the UniProt entry name (KPYR_HUMAN), the PDB codes printed in Valentini *et al.* (2002) have been superseded, and only four *PKLR* entries in UniProt are manually reviewed (human, rat, mouse, and dog); our contingency plan consists of using reviewed entries only, which is why we compare the rat and the mouse and use the dog as outgroup, recording every accession number with its access date and the reason for its choice (tasks 2.2, 2.5, 2.6), and having a second member verify each record and the variant table (task 2.8).
- **R8 (6).** Our identified risk is *developing the parts of the project in isolation, so that each member understands only their own tasks*, we have realized because the rubric marks as unattained a project generated in a disconnected way and a defence showing "niche awareness"; our contingency plan consists of a synthesis meeting before writing (task 3.8), a full read-through of the report by every member (task 4.7), each member preparing and presenting their own part of the defence (tasks 4.10 to 4.14), and a rehearsal in which each member answers questions on the parts of others (task 4.16).
- **R10 (6).** Our identified risk is *failing to meet the report format, its length limits, or the delivery deadline*, we have realized because the report must follow the *BMC Bioinformatics* structure with at most 7 pages or 5 000 words and five figures or tables, and contributions committed after 18:00 on 28 October are discarded; our contingency plan consists of a word budget per section, a format and reference check on 27 October (task 4.8), assembling the report two days before the deadline (task 4.6), and committing the final version by 12:00 on 28 October, six hours before the deadline (task 4.9).

**Medium risks**

- **R3 (3).** Our identified risk is *a repository structure that does not comply with the requirements*, we have realized because the Formal Requirements require the folders `report/`, `status/`, `gantt/`, and `references/`, which our repository does not yet contain; our contingency plan consists of creating them on 9 October (task 1.1) and checking the structure at each weekly meeting (task 1.7).
- **R4 (4).** Our identified risk is *not keeping the Project Status document and the live Gantt chart up to date*, we have realized because both must be maintained until the end of the project and the rubric assesses the comparison between planned and real dates; our contingency plan consists of creating the status document on 13 October (task 1.6), updating both at each weekly meeting (task 1.7), and producing the final chart with real dates (task 1.9).
- **R6 (4).** Our identified risk is *misinterpreting alignment or structural results*, we have realized because we learn these methods during the project (labs on 8, 15, and 22 October) and the rubric marks incorrect interpretation as unattained; our contingency plan consists of scheduling each analysis after the lab that teaches it, asking the instructors during the labs, having a second member repeat the multiple alignment (task 3.9), and repeating the superposition in the 22 October lab (task 3.4).
- **R9 (4).** Our identified risk is *supporting statements with weak or unverified sources, including GenAI output*, we have realized because the rubric requires trusted, relevant sources deposited in the repository, and the GenAI guidelines make us liable for every error we include; our contingency plan consists of using primary articles for facts and reviews for context, verifying every statement suggested by GenAI against its original source, and archiving every cited source in `references/` (task 2.7).

**Low risk**

- **R7 (2).** Our identified risk is *an insufficiently justified species comparison, or evolutionary changes inferred without an outgroup*, we have realized because the rubric requires the chosen species to be justified and useful, and three sequences alone cannot show which residue is ancestral; our contingency plan consists of stating the reason for each species in the report, aligning the dog sequence as an outgroup (tasks 2.5, 3.2), and making claims about the lineage of a change only when the outgroup supports them.

### Monitoring

The register is copied into the Project Status document (task 1.6). At each weekly meeting (task 1.7), the owner of each risk reports its status, any change in its score, the actions taken, and whether the contingency plan was activated.

---
[← 2. Gantt chart: task details](PKD_4_Task_Details.md) · [Index](README.md) · [References →](PKD_6_References.md)
