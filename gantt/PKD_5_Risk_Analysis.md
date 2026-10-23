## 3. Risk analysis

Risk management means identifying uncertain events that may affect the project, assessing their importance, preparing a response and monitoring them over time. A risk is not a problem: a problem has already happened, whereas a risk may happen in the future. For this reason we carried out this analysis while planning the project, following four steps: (1) identify the risk, (2) assess its likelihood and impact, (3) define a contingency plan, and (4) monitor its status.

This section only lists **external risks**, that is, events outside our control that could affect the project (the availability of services and data, the course calendar, or personal circumstances). Risks that depend on our own work habits, such as the distribution of tasks or the frequency of commits, are managed through the project plan and the Project Status document. We identified the risks from the Disease Project Overview, from the course schedule, and from the tools and databases our analysis relies on. Each risk is scored as **Risk = Likelihood × Impact**, both from 1 (low) to 3 (high): **1–2 Low**, **3–5 Medium**, **6–9 High**. Each risk has an owner, who monitors it and activates the contingency plan when needed.

### Risk matrix

| Likelihood ↓ / Impact → | 1 (low) | 2 (moderate) | 3 (high) |
| --- | --- | --- | --- |
| **3 (likely)** | — (3) | R4 (6) | — (9) |
| **2 (possible)** | R9 (2) | R5, R6 (4) | R1, R2, R3 (6) |
| **1 (unlikely)** | — (1) | R8 (2) | R7 (3) |

### Risk register

| ID | External risk | L | I | Score | Level | Owner | Status on 18 Oct | Source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | GitHub unavailable or access problems near a deadline | 2 | 3 | **6** | High | Ettore | Open | Overview, Timelines |
| R2 | Public databases (UniProt, PDB, NCBI) unavailable, or their records updated or superseded during the project | 2 | 3 | **6** | High | Bruno | Open | Overview, Content Integration |
| R3 | Illness or unforeseen personal circumstances make a member unavailable | 2 | 3 | **6** | High | Pol | Open | Overview, Goal (team project) |
| R4 | Exams and deadlines from other courses overlap with the project period | 3 | 2 | **6** | High | Leo | Open | Overview, Timelines |
| R5 | Online analysis tools or servers (alignment, structure viewers) unavailable, overloaded or changed | 2 | 2 | **4** | Medium | Aarón | Open | Overview, Content Integration |
| R6 | Little published evidence, or no experimental structures, for the variants or species studied; articles behind paywalls | 2 | 2 | **4** | Medium | Pol | Open | Overview, Scope Expectations |
| R7 | Course requirements, deadlines or instructions are changed or clarified during the project | 1 | 3 | **3** | Medium | Ettore | Open | Overview, Timelines |
| R8 | A practical session (8, 15 or 22 October) is cancelled or rescheduled | 1 | 2 | **2** | Low | Leo | Open | Overview, Content Integration |
| R9 | Failure of a member's computer or internet connection | 2 | 1 | **2** | Low | Aarón | Open | Overview, Timelines |

### Analysis and contingency plans

**High risks**

- **R1 (6).** Our identified risk is *GitHub being unavailable or inaccessible near a deadline*, we have realized because the whole project is developed and assessed in GitHub, content committed after the deadline is discarded, the plan was delivered on 9 October, and the final report must be submitted by 26 October; our contingency plan consists of finishing each delivery at least one day before its deadline, keeping a local clone of the repository on every member's computer, committing at the end of every work session so that little work is ever pending, and contacting the instructors immediately with the evidence of the outage if it affects a deadline.
- **R2 (6).** Our identified risk is *the public databases being unavailable, or their records being updated or superseded while we work*, we have realized because our analysis depends on UniProt, the PDB and genomic databases, whose entries and PDB codes are revised over time, and some codes printed in the literature have already been superseded; our contingency plan consists of recording every accession number with its access date, downloading and archiving in the repository the sequences and structure files we use, checking the current record before citing it, and using an alternative database or mirror if one is unavailable.
- **R3 (6).** Our identified risk is *a member becoming unavailable because of illness or unforeseen personal circumstances*, we have realized because the project lasts several weeks, is developed by five members and requires every member to defend their own part; our contingency plan consists of assigning to every task a second member who knows its status, documenting each task in the Project Status document so that others can continue it, redistributing the pending tasks in the next weekly meeting, and updating the Gantt chart with the new dates.
- **R4 (6).** Our identified risk is *exams and deadlines from other courses overlapping with the project*, we have realized because the project runs through October with fixed dates (the plan on 9 October, the discussion of the Project Status on 22 October, and the final report on 26 October) and all members have other academic obligations; our contingency plan consists of collecting each member's known exam and delivery dates at the first meeting, scheduling the heaviest tasks away from those dates, leaving a margin of at least two days before each delivery, and rebalancing tasks between members when a conflict appears.

**Medium risks**

- **R5 (4).** Our identified risk is *online analysis tools or servers being unavailable, overloaded or changed*, we have realized because the multiple alignment and the structural analysis rely on web services and software we do not control, and the labs we need are concentrated in October; our contingency plan consists of running each analysis soon after the lab that teaches it, knowing an alternative tool for each analysis (for example, a second alignment server), saving the input files and the parameters used, and exporting every result to the repository as soon as it is obtained.
- **R6 (4).** Our identified risk is *finding little published evidence or experimental structures for the variants or species we study, or being unable to read articles behind paywalls*, we have realized because the Overview asks us to compare humans with at least two other species and explain the impact of the differences, and the literature on non-human proteins is usually smaller; our contingency plan consists of searching for sources early, using the university library access, preferring open-access articles and reviews for context, using predicted structures when no experimental structure exists and stating this clearly in the report, and justifying any topic that must be left out of the scope, as the Overview allows.
- **R7 (3).** Our identified risk is *course requirements, deadlines or instructions changing or being clarified during the project*, we have realized because the instructors may update the course documents (Formal Requirements, template) and the dates in different documents may not coincide; our contingency plan consists of checking the course announcements and documents at every weekly meeting, asking the instructors in class about any ambiguity, recording the confirmed dates in the Project Status document, and adapting the plan and the Gantt chart as soon as a change is confirmed.

**Low risks**

- **R8 (2).** Our identified risk is *a practical session being cancelled or rescheduled*, we have realized because the labs on 8, 15 and 22 October teach the methods we apply; our contingency plan consists of using the course materials and the instructors' tutoring hours to learn the method on our own, and moving the dependent task to the first available date.
- **R9 (2).** Our identified risk is *the failure of a member's computer or internet connection*, we have realized because all work is done in GitHub and online tools; our contingency plan consists of keeping all work in the remote repository, using the university computers or a teammate's device as an alternative, and informing the team so that the affected task can be reassigned if needed.

### Monitoring

Risk monitoring is documented in this register, which is copied into the Project Status document and updated at each weekly meeting. At each meeting, the owner of every risk reports its current status, the actions taken, any change in its score, and whether the contingency plan was activated. The register is a shared, evidence-based plan that we review and improve as the project develops.

| ID | Owner | Current status | Actions taken | Change in score | Contingency plan activated |
| --- | --- | --- | --- | --- | --- |
| R1 | Ettore | Open | — | — | No |
| R2 | Bruno | Open | — | — | No |
| R3 | Pol | Open | — | — | No |
| R4 | Leo | Open | — | — | No |
| R5 | Aarón | Open | — | — | No |
| R6 | Pol | Open | — | — | No |
| R7 | Ettore | Open | — | — | No |
| R8 | Leo | Open | — | — | No |
| R9 | Aarón | Open | — | — | No |

---

[← 2. Gantt chart: task details](PKD_4_Task_Details.md) · [Index](README.md) · [References →](../references/PKD_6_References.md)