# Structure analysis: p.Arg486Trp in human R-type pyruvate kinase

These notes document how we located residue 486 in the crystal structures of human R-type pyruvate kinase (PK-R), compared the wild-type and R486W structures, and compared the human site with the equivalent residue of rat L-type pyruvate kinase. All database records were accessed on 22 October 2026. All numbers were computed from the deposited coordinates with the scripts in `scripts/`. Paths are relative to this folder (`analysis/structure/`).

Numbering: human R-type isoform (UniProt P30613-1, 574 residues) unless stated otherwise. L-type numbering = R-type numbering − 31, confirmed by pairwise alignment for human P30613-2 and rat P12928-2 (`data/rat_mapping.txt`). The author residue numbering of 2VGB and 2VGI equals the R-type numbering (`_struct_ref_seq`: database 47–574 = author 47–574), so no conversion is needed for these entries.

## 1. Aim

1. Select the experimental structures of the wild-type and R486W enzyme and of the rat orthologue.
2. Describe the environment of Arg486 in the wild-type structure: domain, secondary structure, solvent accessibility, interactions, subunit interfaces, and distances to the functional sites.
3. Compare the wild-type and R486W structures.
4. Compare the human site with rat Gln455 (= rat R-type Gln486).
5. State what the structures can and cannot show, and how the results compare with Valentini *et al.* (2002).

## 2. Data sources

| Source | Records (version) | Files |
|---|---|---|
| RCSB PDB, human R-type wild type | 2VGB (revision 2.0, 1 May 2024) | `raw/2VGB.cif`, `raw/2vgb_validation.xml.gz` |
| RCSB PDB, human R486W | 2VGI (revision 2.0, 1 May 2024) | `raw/2VGI.cif`, `raw/2vgi_validation.xml.gz` |
| RCSB PDB, rat L-type | 6ECH, 6ECK (revision 1.4, 9 and 23 October 2024) | `raw/6ECH.cif`, `raw/6ECK.cif`, `raw/6ech_validation.xml.gz` |
| RCSB PDB, independent check | 8TBS (human R-type construct with activator AG-946, 2.35 Å) | not distributed (Section 7.1); results in `data/site_8TBS_486.json` |
| RCSB PDB, survey | all entries mapped to UniProt P30613 and P12928 | `data/rcsb_search_*.json`, `data/rcsb_survey*.tsv`, `data/rcsb_survey.json` |
| UniProtKB | P30613 (entry version 260), P30613-2; P12928 (entry version 209, sequence version 2), P12928-2 | `raw/uniprot_*.json`, `raw/uniprot_*.fasta` |
| InterPro | domain annotations of P30613 and P12928 | `raw/interpro_*.json` |
| Literature | Valentini *et al.* 2002, J Biol Chem 277:23807–14, doi:10.1074/jbc.M202107200 (CC BY) | not distributed |

The wwPDB validation reports (PDF) are not distributed here. They are linked from the entry pages (https://www.rcsb.org/structure/2VGB, https://www.rcsb.org/structure/2VGI, https://www.rcsb.org/structure/6ECH) and can be downloaded directly (gzip-compressed):

- https://files.rcsb.org/pub/pdb/validation_reports/vg/2vgb/2vgb_full_validation.pdf.gz
- https://files.rcsb.org/pub/pdb/validation_reports/vg/2vgi/2vgi_full_validation.pdf.gz
- https://files.rcsb.org/pub/pdb/validation_reports/ec/6ech/6ech_full_validation.pdf.gz

Software: Biopython 1.88 (all geometry unless stated); FreeSASA 2.2.1 (accessibility); PyMOL 3.1.8 and UCSF ChimeraX 1.10.1 (cross-checks and figures).

## 3. Steps and parameters

1. **Metadata** (`scripts/01_metadata.py`): method, resolution, R-work/R-free, entities, ligands, chains, observed residues, revisions.
2. **PDB survey** (`scripts/02_rcsb_survey.py`, `scripts/02b_isoform_classify.py`): all entries for P30613 and P12928 through the RCSB GraphQL API; classification by the isoform-specific N-terminus present in the construct; `pdbx_mutation` field of every entry.
3. **Site analysis** (`scripts/03_site_analysis.py`, helpers in `scripts/pkcommon.py`): domain (boundaries of Valentini *et al.* 2002, p. 23809: N-terminal 57–84, A 85–159 and 263–431, B 160–262, C 432–574); polar contacts ≤ 3.5 Å (heavy atoms); salt bridges (Arg N to Asp/Glu O ≤ 4.0 Å); all contacts ≤ 4.0 Å; relative SASA (FreeSASA, Lee–Richards algorithm, ProtOr radii, probe 1.4 Å, protein atoms only, relative to the FreeSASA Ala-X-Ala reference); area buried by other subunits; minimum distances to ligands and to other chains.
4. **Validation** (`scripts/04_validation_local.py`): per-residue wwPDB metrics (RSCC, RSRZ, rotamer, Ramachandran, clashes) for residue 486 and 13 neighbouring residues.
5. **Superposition of wild type and R486W** (`scripts/05_superpose_wt_mut.py`): Biopython `Superimposer` (least-squares, Kabsch), atoms paired by chain identifier, author residue number, and atom name, no outlier rejection; global, per-domain, and local (residues with any atom within 10 Å of residue 486) fits. Cross-check with PyMOL `align`/`super` (`scripts/07_pymol_align.py`).
6. **Baseline** (`scripts/06_local_baseline.py`): the same local 10 Å measure centred on every N-, A-, and C-domain residue, to compare the differences around residue 486 with those elsewhere in the two models.
7. **Secondary structure**: PyMOL `dss` (`scripts/08_secondary_structure.py`) and DSSP in ChimeraX (`scripts/09_chimerax_dssp.cxc`); deposited HELX_P records.
8. **Rat comparison**: mapping of human 486 onto rat L-type and the 6ECH/6ECK SEQRES (Biopython `PairwiseAligner`, global, BLOSUM62, gap open −10, gap extension −0.5; `scripts/10_rat_mapping.py`); superposition with Cα pairs matched by rat L-type number + 31 = human R-type number (`scripts/11_rat_vs_human.py`).
9. **Cation–π geometry** (`scripts/12_cation_pi_geometry.py`): Arg486 CZ to the Tyr487 ring centroid, and the angle to the ring normal.
10. **Figures** (`scripts/fig2_pymol.py`, `scripts/fig3_pymol.py`, `scripts/compose_figures.py`).

## 4. Results

A summary table of all measurements is in `TABLE_structural_measurements.md`.

### 4.1 Structures

| | 2VGB (wild type) | 2VGI (R486W) |
|---|---|---|
| Title | Human erythrocyte pyruvate kinase | Human erythrocyte pyruvate kinase: R486W mutant |
| Method; resolution (mmCIF) | X-ray diffraction; 2.73 Å | X-ray diffraction; 2.87 Å |
| R-work / R-free (mmCIF) | 0.227 / 0.271 | 0.257 / 0.311 |
| R / R-free in the paper (Table I, p. 23808) | 23.0% / 27.9%, 2.7 Å | 24.3% / 30.9%, 2.9 Å |
| Space group; asymmetric unit | P2₁; one tetramer (chains A–D) = biological assembly 1 | same |
| Ligands per subunit | FBP (fructose 1,6-bisphosphate), PGA (2-phosphoglycolate, a PEP analogue), K⁺, Mn²⁺ | same |
| Waters | 65 | 0 |
| Construct (SEQRES) | P30613 residues 47–574 (528 residues), expressed in *Escherichia coli* | same, with R486W in all four chains ("engineered mutation") |
| Residues observed | 57–573 in every chain; gaps only in the B domain (chain B 167–174, 187–196, 229–236; chain D 167–171) | 57–573; gaps in chain B (eight segments within 167–260), chain C (258–259), chain D (168–171) |
| Deposition; current revision | 12 November 2007; 2.0 (1 May 2024) | 13 November 2007; 2.0 (1 May 2024) |

- The 2002 paper deposited 1LIU (wild type), 1LIW (T384M), 1LIX (R486W), and 1LIY (R479H) (p. 23807). RCSB lists these as removed on 20 November 2007 and replaced by 2VGB, 2VGF, 2VGI, and 2VGG (`data/obsolete_1LIx.txt`). The refinement statistics of the current entries differ slightly from Table I of the paper, so the current entries are not identical to the 2002 models.
- The paper states that the crystallized protein lacks the first 49 residues (Met + residues 50–574; pp. 23808–9), whereas the SEQRES of 2VGB/2VGI lists residues 47–574. Residues 47–56 are not modelled, so the difference does not affect any analysis here.
- The modelled residues of every 2VGB chain are identical to the supplied WT protein at all observed positions (517, 491, 517, and 512 residues compared); 2VGI differs only by R486W in each chain.
- In 2VGI all residual B-factors are 2.0 Å² (refinement with TLS groups; `_refine.B_iso_mean` 2.00), so per-atom B-factors of the mutant carry no local information. In 2VGB the B-factors are identical in the four chains, consistent with the tight non-crystallographic symmetry restraints stated in Table I of the paper; the four chains of each entry are therefore not independent observations.

**Other entries** (RCSB search for P30613 and P12928, 81 entries; `data/rcsb_survey.tsv`, `data/rcsb_survey_isoform.tsv`): 2VGF (T384M, 2.75 Å) and 2VGG (R479H, 2.74 Å) from the same study; R-type constructs (residues ≥ 50) with activators, 8TBS (AG-946, 2.35 Å), 8TBT, 8TBU; cryo-EM structures of full-length R-type, 9W4Q, 9W7Y, 9W81, 9WDO, 9WDQ, 9WDR, 9WDT (mitapivat), 9WEZ (FBP), 2.57–2.93 Å, and 9VXP (2.48 Å; construct lacking the isoform-specific N-terminus); about 65 X-ray entries of the human L-type isoform (e.g. 4IMA, 4IP7, 6NN4–6NN8, 7QDN, 8XFD with mitapivat), in which the residue equivalent to R-type Arg486 is Arg455. No entry other than 2VGI carries a substitution at residue 486. AlphaFold models were not needed, because experimental structures cover residues 57–573 of the R-type and the whole rat L-type region analysed.

**Rat L-type structures** (Gassaway *et al.* 2019, Cell Rep 29:3394–404.e9; only the PDB records were used):

| | 6ECH | 6ECK |
|---|---|---|
| Title | L-type PK with phosphorylated Ser12 (pS12) in complex with FBP | L-type PK with phosphorylated Ser113 (pS113) in complex with FBP |
| Resolution; R-work / R-free | 2.19 Å; 0.181 / 0.238 | 2.36 Å; 0.204 / 0.275 |
| Asymmetric unit; assembly | tetramer (A–D) = assembly 1 | dimer (A, B); tetramer generated by a crystallographic operator |
| Ligands | FBP (4), Mg²⁺ (13), ATP (3), oxalate (4, a PEP/pyruvate analogue), ethylene glycol and PEG fragments, Cl⁻ | FBP (2), citrate, glycerol, PEG fragments |
| Phosphorylation | Ser12 phosphorylated in SEQRES, but residues 1–17 are not modelled | pSer113 modelled in both chains |
| Construct | rat L-type 1–543 (P12928-2) with His-tag | same |

We used 6ECH for the comparison (higher resolution, tetramer in the asymmetric unit, active-site ligands oxalate, Mg²⁺, and ATP in addition to FBP) and 6ECK as a check of the Gln455 environment.

### 4.2 Arg486 in the wild-type structure (2VGB)

| Property | Result |
|---|---|
| Domain | C domain (432–574). CATH-Gene3D 3.40.1380.20 (429–566) and Pfam PF02887 PK_C (454–571) in InterPro agree |
| Secondary structure | C-terminal end of helix 479–487: DSSP (ChimeraX) α-helix 479–484 followed by a 3₁₀ segment 485–487; PyMOL `dss` helix 479–487; deposited HELX_P 478–486. The paper describes this helix as 480–486 (p. 23810) and 479–486 (p. 23812), with its N-terminus pointing to the 1-phosphate of FBP |
| Relative SASA | 0.11–0.12 (side chain 0.09–0.11): largely buried within its own subunit |
| Polar partners ≤ 3.5 Å | only the main-chain carbonyl O of Leu362 (A domain, same chain): NH2–O 2.74–2.90 Å in all four chains; NH1–O 3.36–3.49 Å (chains A, B) |
| Weaker contacts 3.5–4.0 Å | Asn361 O, Ala363 O, Gly364 N |
| Salt bridge | none (nearest carboxylate Glu453, 4.3–4.4 Å) |
| Tyr487 ring | centroid 5.1–5.2 Å from Arg486 CZ at 28–31° from the ring normal: at the long end of a cation–π geometry; not counted as an established interaction |
| Waters | none modelled within 4 Å |
| Interdomain contact | all side-chain contacts ≤ 4 Å are with A-domain residues 361–364 or C-domain residues 441, 445, 482–488: Arg486 lies at the interface between the A and C domains of one subunit (as stated by Valentini *et al.* 2002, p. 23811) |
| Subunit interfaces | no contact with another subunit: 0 Å² buried by other subunits; nearest atoms of other chains 6.1–6.5 Å (Lys465 NZ of the C/C′ partner) and 6.6–7.0 Å (Leu70 CD2 of the A/A′ partner). Tyr487 and Leu362 are themselves interface residues (C/C′ and A/A′) |
| Minimum distances, same subunit | FBP 12.4–12.8 Å; phosphoglycolate (active site) 24.3–24.6 Å; Mn²⁺ 26.6–26.9 Å; K⁺ 27.5–27.9 Å; Lys313 (transition-state stabilizer, UniProt) 19.5–19.6 Å; Arg116 and Thr371 22.5–22.9 Å |
| Model quality | wwPDB validation: rotamer outlier in all four chains; Ramachandran outlier in chain A; RSCC 0.83–0.89 (neighbouring residues 0.76–0.98, median 0.92) |

Tetramer interfaces: A–B and C–D are C/C′ interfaces (26–29 C-domain residues); A–C and B–D are A/A′ interfaces (mostly A- and N-domain residues).

**Independent check, 8TBS** (2.35 Å, eight chains): Arg486 lies at the C-terminal end of helix 478–486, buried (relative SASA 0.11–0.19), 5.6–6.7 Å from the nearest other subunit, 12.4–12.6 Å from FBP, and 26.4–27.7 Å from Mn²⁺. Its guanidinium hydrogen-bonds Leu362 O in seven of eight chains (2.47–3.03 Å) and Glu453 OE2 in chains A and H (3.09 and 2.55 Å). The Arg486–Leu362 carbonyl hydrogen bond is therefore reproducible; a salt bridge with Glu453 occurs only in some chains of 8TBS.

### 4.3 Wild type versus R486W

| Measure | Result |
|---|---|
| Global Cα RMSD (2VGI fitted on 2VGB) | chain A 0.587 Å (517 Cα), B 0.594 Å (463), C 0.523 Å (515), D 0.401 Å (512); tetramer 0.765 Å (2007 Cα) |
| Per-domain fits | N, A, and C domains 0.21–0.24 Å; B domain 0.42–0.66 Å |
| PyMOL cross-check | `align`/`super` with `cycles=0` reproduce the Biopython values; with the default outlier rejection (5 cycles, 2 Å) 0.35–0.40 Å per chain and 0.52 Å for the tetramer (1741–1758 Cα retained) |
| Local RMSD (residues within 10 Å of residue 486; 58–60 residues, three chains) | Cα 0.30–0.36 Å after local fit (0.41–0.55 Å after the tetramer fit); all common heavy atoms 0.43–0.48 Å; backbone and Cβ of residue 486 move ≤ 0.56 Å |
| Largest local shifts | side chains of Glu453 (OE2 1.56–1.75 Å, towards Trp486), Leu441, Arg443, Val506, Arg510 (1.0–1.5 Å) |
| Baseline | the same local measure centred on every N/A/C-domain residue: median Cα RMSD 0.24–0.25 Å, 95th percentile 0.34–0.41 Å per chain; the sphere around residue 486 ranks at the 74th–95th percentile, and its largest atom shift at the 59th–87th percentile |
| Trp486 interactions | indole NE1 2.68–3.40 Å from Asn361 O and 2.82–3.37 Å from Leu362 O; CZ3 2.99–3.18 Å from Glu453 OE1 (C···O); NE1 2.89–3.28 Å from Gly364 N; in chain D, CZ2 2.41 Å from Leu362 O (clash flagged by wwPDB, overlap 0.73 Å) |
| Trp486 burial | relative SASA 0.14–0.18; 0.5–1.2 Å² buried by neighbouring subunits (Lys465 NZ of the C/C′ partner at 5.2–5.7 Å): no new subunit contact |
| Torsions | χ1 of residue 486 −143° to −152° (Arg) → −169° to −174° (Trp); χ1 of Asn361 shifts by 10–18° and of Glu453 by 10–13°; φ of residue 486 −34° to −41° (WT) → −52° to −65° (mutant), all helical |
| B-factors and density | not comparable (2VGI residual B fixed at 2.0 Å²; 2VGI RSCC ≈ 0.99 everywhere). No residue near 486 is missing in either model |

The global values are dominated by the B domain, which the paper describes as flexible (orientations differing by about 4° between subunits, p. 23809). The differences around residue 486 are within the range of differences found throughout the two models (independent refinements at 2.7 and 2.9 Å resolution) and cannot be attributed specifically to the substitution. Trp486 is accommodated at the A/C interdomain contact with conserved backbone; the hydrogen bond to the 361–362 carbonyls is retained with a different donor, and neighbouring side chains (Glu453) adjust.

### 4.4 Comparison with Valentini *et al.* (2002)

Kinetic parameters of the recombinant enzymes (Table II, p. 23810; 37 °C; means ± S.E.; PEP parameters at 1.5 mM ADP from Hill plots; ADP parameters at 5 mM PEP). The paper reports S₀.₅, not K_m; the ratios were computed by us.

| Parameter | Wild type | R486W | R486W / WT |
|---|---|---|---|
| PEP, no FBP: k_cat (s⁻¹) | 355 ± 12 | 195 ± 4 | 0.55 |
| S₀.₅ PEP (mM) | 1.10 ± 0.04 | 1.69 ± 0.06 | 1.5 |
| k_cat/S₀.₅ (s⁻¹ mM⁻¹) | 323 | 116 | 0.36 |
| n_H | 1.60 ± 0.16 | 2.07 ± 0.11 | – |
| PEP, 1 mM FBP: k_cat (s⁻¹) | 355 ± 11 | 203 ± 10 | 0.57 |
| S₀.₅ PEP (mM) | 0.18 ± 0.020 | 0.40 ± 0.07 | 2.2 |
| k_cat/S₀.₅ (s⁻¹ mM⁻¹) | 1972 | 492 | 0.25 |
| n_H | 1.05 ± 0.07 | 1.32 ± 0.03 | – |
| ADP: k_cat (s⁻¹) | 355 ± 13 | 218 ± 13 | 0.61 |
| S₀.₅ ADP (mM) | 0.17 ± 0.01 | 0.24 ± 0.02 | 1.4 |
| k_cat/S₀.₅ (s⁻¹ mM⁻¹) | 2080 | 908 | 0.44 |

Statements of the paper (paraphrased): Arg486 lies in a region of close association between the A and C domains and is hydrogen-bonded to the carbonyl oxygen of Leu362 (p. 23811). In the mutant the Trp side chain is accommodated without structural perturbation, no atomic movement larger than 0.25 Å is detected, and the indole nitrogen hydrogen-bonds the carbonyl oxygens of Leu362 and "Asp361" (p. 23812); residue 361 is Asn in P30613, in the supplied sequence, and in both PDB models. The thermostability of R486W at 53 °C (Fig. 3B, p. 23811) is "even higher" than that of the wild type (p. 23812). The allosteric properties are essentially unmodified; the only significant change is the catalytic efficiency, which drops to about 30% of the wild-type value (p. 23812); FBP still lowers S₀.₅ for PEP (1.69 to 0.40 mM) and n_H (2.07 to 1.32) in R486W. Arg486 is ≥ 20 Å from the catalytic centre, which is unperturbed; the authors suggest that the kinetic effect "might reflect altered dynamic properties" (p. 23812; a hypothesis). Patients homozygous for R486W generally have mild anaemia (p. 23812).

- **Agreement**: C domain; A/C interdomain location; Arg486 NH2–Leu362 O hydrogen bond (2.74–2.90 Å); Trp486 NE1 contacts with the 361 and 362 carbonyls; distance ≥ 19.5 Å from the active site; unchanged active site.
- **Difference**: in the current deposited models (2VGB, 2VGI, deposited in 2007) local displacements exceed 0.25 Å (local Cα RMSD 0.30–0.36 Å; side-chain shifts up to 1.75 Å). They are of the same size as differences elsewhere in the two models, so they do not contradict the paper's conclusion that there is no substitution-specific structural perturbation. The paper does not state how its 0.25 Å figure was computed, and its models (1LIU, 1LIX) are no longer distributed.
- **Additional observations**: Arg486 is not at a subunit interface (no contact ≤ 4 Å with other chains; 6.1–7.0 Å from the nearest atoms); its distance to bound FBP is 12.4–12.8 Å; the Glu453 side chain is displaced towards Trp486; the deposited Arg486 conformation is a rotamer outlier.

### 4.5 Rat L-type Gln455 (= rat R-type Gln486)

- **Equivalent residue.** Global alignment of the human R-type sequence with rat L-type P12928-2, rat R-type P12928-1, and the SEQRES of 6ECH/6ECK aligns human Arg486 with rat L-type Gln455 (rat R-type Gln486); author residue 455 of 6ECH/6ECK is Gln in all chains (author number = SEQRES index − 1). Window: human R 476–496 TTGRSAQLLS**R**YRPRAAVIAV; rat L 445–465 KTGRSAQLLS**Q**YRPRAAVIAV. Differences in human 456–516 versus rat: A463S, T476K, R486Q, C509S. Mouse and dog carry Arg at this position (`../sequence/`).
- **Fold and position.** Gln455 is in the C domain at the C-terminal end of helix 447–455 (equivalent to human 478–486), relative SASA 0.16–0.20, 0 Å² buried by other subunits, 6.0–6.4 Å from the C/C′ partner and 7.9–8.4 Å from the A/A′ partner, 11.9–12.0 Å from FBP, and 23.9–24.2 Å from oxalate.
- **Superposition** (Cα pairs by rat L number + 31 = human R number): subunit A of 2VGB against the rat subunits, all 517 Cα 0.98–1.95 Å; N+A+C core 0.35–0.52 Å; B domain alone 0.69–0.74 Å; best tetramer mapping 1.73 Å (2037 Cα). Local 10 Å sphere 0.42 Å (59 Cα); the Cα of rat 455 is 0.69 Å from the Cα of human 486. The higher whole-subunit RMSD reflects different B-domain orientations relative to the core.
- **Local environment.** Gln455 NE2 hydrogen-bonds Glu422 OE2 (2.86–3.10 Å; = human Glu453) and the main-chain O of Leu452 (2.88–3.25 Å; = human Leu483, i−3 in the same helix); OE1 contacts a water (2.30–2.54 Å). The A-domain carbonyls of Asn330/Leu331 (= human Asn361/Leu362) are 5.5–6.7 Å from the Gln side chain. The same pattern is present in both chains of 6ECK (NE2–Glu422 OE2 2.73–3.09 Å; NE2–Leu452 O 3.20–3.37 Å). The amide O/N assignment of Gln is not determined by the electron density at 2.2 Å; we inferred it from the hydrogen-bond pattern.
- **Interpretation (inference).** In rat L-type PK the Gln side chain does not form the A/C interdomain hydrogen bond made by human Arg486; it is satisfied within the C domain (Glu422, Leu452 O) and by solvent. A non-Arg residue at this position is therefore compatible with a folded, FBP-bound, catalytically competent rat L-type enzyme (the structure contains ATP, Mg²⁺, and oxalate in the active site). This does not show that Trp is tolerated in the human R-type enzyme: rat L-type differs from human R-type at other positions and in its N-terminus, and the kinetic consequences of Gln at this position have not been measured in the sources we read.

## 5. What the structures can and cannot show

**Can show (measured):** the location of residue 486 in the C domain at the A/C interdomain contact of one subunit; its burial; its hydrogen-bond partners; its distances from bound FBP (≈ 12 Å), the active site (≈ 24 Å from phosphoglycolate, ≈ 20 Å from Lys313), and other subunits (≥ 6 Å); that the mutant backbone and active site are unchanged within the precision of 2.7–2.9 Å structures; the side-chain contacts of Trp486 and of rat Gln455.

**Cannot show:** thermodynamic or kinetic stability (Valentini *et al.* measured R486W as more thermostable than wild type at 53 °C); conformational dynamics or the populations of inactive (T-state) conformations, because both structures are of FBP- and substrate-analogue-bound enzyme; effects in erythrocytes (protein turnover, expression); coordinate shifts smaller than the coordinate error at 2.7–2.9 Å; B-factor changes (2VGI has no informative B-factors). The four chains of each entry were refined with tight NCS restraints and are not independent replicates.

## 6. Notes

1. The 0.25 Å statement of Valentini *et al.* cannot be reproduced exactly from the current entries; the original 1LIU and 1LIX coordinates (PDB obsolete archive) would be needed.
2. The deposited Arg486 conformation of 2VGB is a rotamer outlier; 8TBS (2.35 Å) supports the Arg486–Leu362 O hydrogen bond, with a Glu453 salt bridge in some chains.
3. The cryo-EM structures of full-length R-type (9WEZ with FBP, 9WDO with mitapivat) were not analysed.
4. The kinetic effect of Gln at this position (rat) was not examined in the literature.

## 7. How to reproduce

### 7.1 Command line

Requirements: Python 3 with the packages in `requirements.txt`; PyMOL 3.1.8 (command-line binary `pymol`; on macOS inside `PyMOL.app/Contents/MacOS/`) for scripts 07, 08, and the figures; UCSF ChimeraX 1.10.1 for script 09. Run everything from this folder.

```sh
cd analysis/structure
python -m pip install -r requirements.txt
# coordinate files that are not distributed (needed for 8TBS checks and metadata of 8XFD)
curl -L -o raw/8TBS.cif https://files.rcsb.org/download/8TBS.cif
curl -L -o raw/8XFD.cif https://files.rcsb.org/download/8XFD.cif
# the distributed files can be re-downloaded in the same way, e.g.
#   curl -L -o raw/2VGB.cif https://files.rcsb.org/download/2VGB.cif
#   curl -L -o raw/2vgb_validation.xml.gz https://files.rcsb.org/pub/pdb/validation_reports/vg/2vgb/2vgb_validation.xml.gz

python scripts/01_metadata.py 2VGB 2VGI 6ECH 6ECK > data/metadata_raw.txt
python scripts/01_metadata.py 8XFD 8TBS > data/metadata_raw_extra.txt
python scripts/02_rcsb_survey.py            # RCSB GraphQL API; reads data/rcsb_search_*.json
python scripts/02b_isoform_classify.py
python scripts/03_site_analysis.py 2VGB ARG
python scripts/03_site_analysis.py 2VGI TRP
python scripts/03_site_analysis.py 8TBS ARG
python scripts/03_site_analysis.py 6ECH GLN 455 31
python scripts/04_validation_local.py 2vgb 2vgi
python scripts/05_superpose_wt_mut.py
python scripts/06_local_baseline.py
python scripts/10_rat_mapping.py             # reads ../sequence/supplied/PKD_protein_WT.fasta
python scripts/11_rat_vs_human.py
python scripts/12_cation_pi_geometry.py > data/cation_pi_geometry.txt
pymol -cq scripts/07_pymol_align.py
pymol -cq scripts/08_secondary_structure.py
ChimeraX --nogui --exit --script scripts/09_chimerax_dssp.cxc

# figures: PyMOL panels (written to figures/panels/), then the composed figures (written to ../../figures/)
pymol -cq scripts/fig2_pymol.py
pymol -cq scripts/fig3_pymol.py
python scripts/compose_figures.py
```

`data/rcsb_search_P30613.json` and `data/rcsb_search_P12928.json` are the RCSB Search API result sets of the entries mapped to each UniProt accession, which `scripts/02_rcsb_survey.py` reads.

### 7.2 With PyMOL or ChimeraX (course tools)

PyMOL (commands typed in the PyMOL command line):

```
fetch 2VGB, type=cif
fetch 2VGI, type=cif
remove solvent
hide everything
show cartoon
color orange, resi 57-84
color blue, resi 85-159 or resi 263-431
color pink, resi 160-262
color green, resi 432-574
show sticks, resi 486 or resn FBP+PGA
show spheres, resn MN+K
# hydrogen bond Arg486 - Leu362 carbonyl (expect 2.7 A in chain A)
distance hb1, 2VGB and chain A and resi 486 and name NH2, 2VGB and chain A and resi 362 and name O
# minimum distances used in Figure 2 (expect 12.4 and 24.4 A)
distance dFBP, 2VGB and chain A and resi 486 and name CB, 2VGB and chain A and resn FBP and name O6P
distance dPGA, 2VGB and chain A and resi 486 and name CA, 2VGB and chain A and resn PGA and name C2
# residues of other chains near residue 486 (expect none within 5 A)
select near486, (2VGB and not chain A) within 5 of (2VGB and chain A and resi 486)
# secondary structure
dss 2VGB
iterate 2VGB and chain A and resi 478-488 and name CA, print(resi, resn, ss)
# superposition without outlier rejection (expect 0.587 A, 517 atoms; tetramer 0.765 A)
align 2VGI and chain A and name CA, 2VGB and chain A and name CA, cycles=0
align 2VGI and name CA, 2VGB and name CA, cycles=0
# mutant contacts (expect 3.0 and 3.4 A)
distance hb2, 2VGI and chain A and resi 486 and name NE1, 2VGI and chain A and resi 361 and name O
distance hb3, 2VGI and chain A and resi 486 and name NE1, 2VGI and chain A and resi 362 and name O
# rat L-type: Gln455 corresponds to human 486
fetch 6ECH, type=cif
align 6ECH and chain D and name CA, 2VGB and chain A and name CA, cycles=0
show sticks, 6ECH and chain D and resi 455+422+452
distance hb4, 6ECH and chain D and resi 455 and name NE2, 6ECH and chain D and resi 422 and name OE2
```

ChimeraX (equivalent):

```
open 2VGB
open 2VGI
dssp #1
info residues #1/A:478-488 attribute ss_type
hbonds #1/A:486 reveal true
contacts #1/A:486 restrict cross reveal true
measure sasa #1/A & protein sum #1/A:486
matchmaker #2/A to #1/A cutoffDistance none
distance #1/A:486@NH2 #1/A:362@O
distance #2/A:486@NE1 #2/A:361@O
```

With the default criterion of `contacts`, no partners from other chains are expected. Small differences from the values reported here are expected for SASA (different radii and reference values) and for superpositions that use sequence alignment or outlier rejection (PyMOL `cycles` > 0; ChimeraX default `cutoffDistance 2.0`).

## 8. Figures

The final figures are in the repository folder `figures/`.

- **Figure 2** (`../../figures/Figure2_PKR_overview_R486.png`, 300 dpi). Legend: Structure of human R-type pyruvate kinase (PDB 2VGB). (a) Tetramer, domains coloured (N-terminal orange, A blue, B reddish purple, C green); black spheres, residue 486 side chains (those of subunits A and C are labelled; B and D are on the far face). (b) Subunit A with FBP, phosphoglycolate (PGA), Mn²⁺ (purple), and K⁺ (grey); dotted lines, minimum interatomic distances from Arg486 to FBP (12.4 Å) and PGA (24.4 Å).
- **Figure 3** (`../../figures/Figure3_R486W_closeup.png`, 300 dpi). Legend: Environment of residue 486. (a) Wild-type Arg486 hydrogen-bonded to the Leu362 carbonyl (2VGB). (b) Trp486 of R486W (2VGI) after local superposition on the wild type (Cα RMSD 0.29 Å over 62 Cα of residues 355–370, 436–460, and 475–495, subunit A); thin dark sticks, wild-type Arg486. (c) Rat L-type Gln455 (6ECH, subunit D) superposed on the human wild type (0.41 Å over the 62 equivalent Cα, rat residues 324–339, 405–429, and 444–464); rat numbering = human R-type − 31. Distances in Å; all panels in the same orientation.

## 9. File index

| Path | Content |
|---|---|
| `README.md` | these notes |
| `requirements.txt` | Python packages |
| `TABLE_structural_measurements.md` | summary table of the measurements for residue 486 (2VGB, 2VGI) and rat Gln455 (6ECH) |
| `scripts/pkcommon.py` | shared helpers (domains, contacts, salt bridges, cation–π) |
| `scripts/01_metadata.py` | mmCIF metadata → `data/metadata_raw.txt`, `data/metadata_raw_extra.txt` |
| `scripts/02_rcsb_survey.py`, `scripts/02b_isoform_classify.py` | PDB survey → `data/rcsb_survey.json`, `data/rcsb_survey.tsv`, `data/rcsb_survey_isoform.tsv` |
| `scripts/03_site_analysis.py` | environment of residue 486 (455 in rat) → `data/site_2VGB_486.json`, `data/site_2VGI_486.json`, `data/site_8TBS_486.json`, `data/site_6ECH_455.json` |
| `scripts/04_validation_local.py` | wwPDB validation around 486 → `data/validation_local.tsv` |
| `scripts/05_superpose_wt_mut.py` | superposition 2VGB/2VGI → `data/superposition_2VGB_2VGI.json` |
| `scripts/06_local_baseline.py` | local baseline → `data/local_baseline.tsv`, `data/local_baseline_summary.json` |
| `scripts/07_pymol_align.py` | PyMOL `align`/`super` cross-check → `data/pymol_align.tsv` |
| `scripts/08_secondary_structure.py` | PyMOL `dss` → `data/ss_pymol.tsv` |
| `scripts/09_chimerax_dssp.cxc` | DSSP in ChimeraX (printed report) |
| `scripts/10_rat_mapping.py` | human 486 → rat 455 → `data/rat_mapping.txt` |
| `scripts/11_rat_vs_human.py` | rat vs human superposition and environment → `data/rat_vs_human.json` |
| `scripts/12_cation_pi_geometry.py` | Arg486–Tyr487 geometry → `data/cation_pi_geometry.txt` |
| `scripts/fig2_pymol.py`, `scripts/fig3_pymol.py`, `scripts/pymol_project.py`, `scripts/compose_figures.py` | Figures 2 and 3 |
| `data/obsolete_1LIx.txt` | replacement of 1LIU–1LIY by 2VGB, 2VGF, 2VGI, 2VGG |
| `data/rcsb_search_P30613.json`, `data/rcsb_search_P12928.json` | RCSB search result sets |
| `raw/` | mmCIF of 2VGB, 2VGI, 6ECH, 6ECK; wwPDB validation XML of 2VGB, 2VGI, 6ECH; UniProt JSON/FASTA of P30613, P30613-2, P12928, P12928-2; InterPro JSON of P30613, P12928 |
| `sessions/fig2b_subunit.pse` | PyMOL session of Figure 2b (subunit A of 2VGB with ligands); a copy of `figures/panels/fig2b_subunit.pse`, which `scripts/fig2_pymol.py` writes |
