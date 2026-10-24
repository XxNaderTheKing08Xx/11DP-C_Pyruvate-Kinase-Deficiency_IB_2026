"""Tasks 2.5/3.1: assemble the core and broad PKLR ortholog panels and confirm orthology.

Inputs (already downloaded into raw/):
  - UniProt Swiss-Prot entries P30613 (human), P12928 (rat), P53657 (mouse), Q29536 (dog)
  - human PKM (P14618 = isoform M2; P14618-2 = isoform M1) and chicken P00548 (KPYM_CHICK)
  - RefSeq proteins of rat/mouse/dog PKLR and the two Pkm records NP_445749.1, NP_035229.2
  - Ensembl Compara (release 116) orthologues of ENSG00000143627: raw/ensembl/compara_PKLR_orthologues.json
Steps:
  1. Writes results/core_panel.fasta (four Swiss-Prot R-type sequences) and
     results/broad_panel.fasta (core four + Ensembl Compara one-to-one orthologues).
  2. Looks up the Ensembl gene symbol of each broad-panel gene (REST POST /lookup/id).
  3. Orthology/paralogy check: global identity (Needleman-Wunsch, BLOSUM62, gap open 10,
     extend 0.5) of each sequence to human PKLR (P30613) and to human PKM (M1 and M2);
     an ortholog must be more similar to human PKLR than to human PKM.
  4. Compares UniProt and RefSeq sequences of rat, mouse and dog.
Output: results/03_orthologs.txt, results/panel_metadata.tsv, results/*.fasta,
        raw/ensembl/lookup_panel_genes.json
"""
import json
import time
from pathlib import Path
import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

BASE = Path(__file__).resolve().parent.parent
RAW, RES = BASE / "raw", BASE / "results"
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


def read1(p):
    return str(next(SeqIO.parse(p, "fasta")).seq).upper()


aligner = PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.end_open_gap_score = 0  # EMBOSS Needle default: end gaps not penalized
aligner.end_extend_gap_score = 0


def pid(a, b):
    """Percentage identity = identical positions / alignment length (EMBOSS Needle definition)."""
    aln = aligner.align(a, b)[0]
    s1, s2 = aln[0], aln[1]
    L = len(s1)
    ident = sum(1 for x, y in zip(s1, s2) if x == y and x != "-")
    return 100 * ident / L


human = read1(RAW / "uniprot/P30613.fasta")
pkm2 = read1(RAW / "uniprot/P14618.fasta")
pkm1 = read1(RAW / "uniprot/P14618-2.fasta")

core = [  # (label, species, accession, file)
    ("Human", "Homo sapiens", "P30613", RAW / "uniprot/P30613.fasta"),
    ("Rat", "Rattus norvegicus", "P12928", RAW / "uniprot/P12928.fasta"),
    ("Mouse", "Mus musculus", "P53657", RAW / "uniprot/P53657.fasta"),
    ("Dog", "Canis lupus familiaris", "Q29536", RAW / "uniprot/Q29536.fasta"),
]
broad_species = [  # (label, Ensembl species key, clade)
    ("Chimpanzee", "pan_troglodytes", "Primates"),
    ("Macaque", "macaca_mulatta", "Primates"),
    ("Marmoset", "callithrix_jacchus", "Primates"),
    ("Guinea pig", "cavia_porcellus", "Rodentia"),
    ("Rabbit", "oryctolagus_cuniculus", "Lagomorpha"),
    ("Cat", "felis_catus", "Carnivora"),
    ("Cow", "bos_taurus", "Artiodactyla"),
    ("Pig", "sus_scrofa", "Artiodactyla"),
    ("Horse", "equus_caballus", "Perissodactyla"),
    ("Opossum", "monodelphis_domestica", "Marsupialia"),
    ("Platypus", "ornithorhynchus_anatinus", "Monotremata"),
    ("Painted turtle", "chrysemys_picta_bellii", "Testudines"),
    ("Wall lizard", "podarcis_muralis", "Squamata"),
    ("Xenopus", "xenopus_tropicalis", "Amphibia"),
    ("Coelacanth", "latimeria_chalumnae", "Actinistia"),
    ("Zebrafish", "danio_rerio", "Actinopterygii"),
    ("Medaka", "oryzias_latipes", "Actinopterygii"),
]
core_clade = {"Human": "Primates", "Rat": "Rodentia", "Mouse": "Rodentia", "Dog": "Carnivora"}

comp = json.loads((RAW / "ensembl/compara_PKLR_orthologues.json").read_text())["data"][0]["homologies"]
by_sp = {}
for h in comp:
    by_sp.setdefault(h["target"]["species"], []).append(h)

records = []  # dicts
for lab, sp, acc, f in core:
    records.append(dict(label=lab, species=sp, clade=core_clade[lab], source="UniProtKB/Swiss-Prot",
                        accession=acc, gene_id="", seq=read1(f), compara_type="(core; see UniProt/RefSeq)"))
for lab, key, clade in broad_species:
    hs = by_sp[key]
    assert len(hs) == 1 and hs[0]["type"] == "ortholog_one2one", (key, [h["type"] for h in hs])
    t = hs[0]["target"]
    records.append(dict(label=lab, species=key.replace("_", " ").capitalize(), clade=clade,
                        source="Ensembl 116 Compara", accession=t["protein_id"], gene_id=t["id"],
                        seq=t["seq"].replace("-", ""), compara_type=hs[0]["type"],
                        compara_perc_id_target=round(t["perc_id"], 1)))
# Compara annotation of the core species (for the record)
for lab, key in [("Rat", "rattus_norvegicus"), ("Mouse", "mus_musculus"), ("Dog", "canis_lupus_familiaris")]:
    h = by_sp[key][0]
    rec = next(r for r in records if r["label"] == lab)
    rec["compara_type"] = h["type"]
    rec["gene_id"] = h["target"]["id"]
    log(f"Compara {lab}: {h['type']} gene {h['target']['id']} protein {h['target']['protein_id']} "
        f"length {len(h['target']['seq'].replace('-', ''))}")

# Ensembl gene symbols of panel genes
ids = [r["gene_id"] for r in records if r["gene_id"]]
for attempt in range(5):
    resp = requests.post("https://rest.ensembl.org/lookup/id", json={"ids": ids},
                         headers={"Content-Type": "application/json", "Accept": "application/json"})
    if resp.ok:
        break
    time.sleep(4)
look = resp.json()
(RAW / "ensembl/lookup_panel_genes.json").write_text(json.dumps(look, indent=1))
for r in records:
    if r["gene_id"]:
        L = look.get(r["gene_id"]) or {}
        r["ens_symbol"] = L.get("display_name", "NA")
        r["ens_description"] = (L.get("description") or "NA").split(" [")[0]
    else:
        r["ens_symbol"] = "PKLR (human ENSG00000143627)"
        r["ens_description"] = "pyruvate kinase L/R"

log("\n== Orthology check: identity to human PKLR vs human PKM (Needle-equivalent global alignment)")
log("label\tlength\t%id_vs_PKLR(P30613)\t%id_vs_PKM2(P14618)\t%id_vs_PKM1(P14618-2)\tcloser_to")
for r in records:
    a, b, c = pid(r["seq"], human), pid(r["seq"], pkm2), pid(r["seq"], pkm1)
    r.update(id_pklr=round(a, 1), id_pkm2=round(b, 1), id_pkm1=round(c, 1), length=len(r["seq"]),
             closer="PKLR" if a > max(b, c) else "PKM")
    log(f"{r['label']}\t{len(r['seq'])}\t{a:.1f}\t{b:.1f}\t{c:.1f}\t{r['closer']}")

log("\n== Control: records that are NOT PKLR orthologs")
for lab, f in [("rat Pkm NP_445749.1", RAW / "ncbi/NP_445749.1.fasta"),
               ("mouse Pkm (M2) NP_035229.2", RAW / "ncbi/NP_035229.2.fasta"),
               ("chicken PKM P00548 (NCBI Gene 396456 symbol 'PKLR')", RAW / "uniprot/P00548.fasta")]:
    s = read1(f)
    a, b, c = pid(s, human), pid(s, pkm2), pid(s, pkm1)
    log(f"{lab}\tlength {len(s)}\tvs PKLR {a:.1f}\tvs PKM2 {b:.1f}\tvs PKM1 {c:.1f}\t-> closer to {'PKLR' if a > max(b, c) else 'PKM'}")

log("\n== UniProt vs RefSeq/Ensembl sequences of the core orthologs")
pairs = [
    ("Rat P12928 (R-type, canonical)", RAW / "uniprot/P12928.fasta", "XP_006232655.1 (RefSeq model, isoform X1)", RAW / "ncbi/XP_006232655.1.fasta"),
    ("Rat P12928 (R-type) residues 32-574", None, "NP_036756.3 (L-type) residues 1-543", RAW / "ncbi/NP_036756.3.fasta"),
    ("Mouse P53657", RAW / "uniprot/P53657.fasta", "NP_038659.2 (isoform 1)", RAW / "ncbi/NP_038659.2.fasta"),
    ("Dog Q29536", RAW / "uniprot/Q29536.fasta", "NP_001243191.1 (isoform 1)", RAW / "ncbi/NP_001243191.1.fasta"),
]
for n1, f1, n2, f2 in pairs:
    if f1 is None:
        s1 = read1(RAW / "uniprot/P12928.fasta")
        s2 = read1(f2)
        log(f"{n1} vs {n2}: rat R-type length {len(s1)}, L-type length {len(s2)}; "
            f"R[34:] == L[3:]: {s1[33:] == s2[2:]}; L-type N-terminus {s2[:2]}")
        continue
    s1, s2 = read1(f1), read1(f2)
    d = [(i + 1, x, y) for i, (x, y) in enumerate(zip(s1, s2)) if x != y] if len(s1) == len(s2) else "length differs"
    log(f"{n1} ({len(s1)} aa) vs {n2} ({len(s2)} aa): identical {s1 == s2}; differences {d}")
for lab, key in [("Rat", "rattus_norvegicus"), ("Mouse", "mus_musculus"), ("Dog", "canis_lupus_familiaris")]:
    t = by_sp[key][0]["target"]
    s = t["seq"].replace("-", "")
    up = next(r["seq"] for r in records if r["label"] == lab)
    log(f"{lab}: Ensembl {t['protein_id']} ({len(s)} aa) identical to UniProt: {s == up}; "
        f"is C-terminal suffix of UniProt: {up.endswith(s)}")

# write FASTA and metadata
with open(RES / "core_panel.fasta", "w") as fh:
    for r in records[:4]:
        fh.write(f">{r['label']}_{r['accession']}\n{r['seq']}\n")
with open(RES / "broad_panel.fasta", "w") as fh:
    for r in records:
        fh.write(f">{r['label'].replace(' ', '_')}_{r['accession']}\n{r['seq']}\n")
cols = ["label", "species", "clade", "source", "accession", "gene_id", "ens_symbol", "ens_description",
        "compara_type", "length", "id_pklr", "id_pkm2", "id_pkm1", "closer"]
with open(RES / "panel_metadata.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in records:
        fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
log("\n== Panel metadata")
log((RES / "panel_metadata.tsv").read_text())
(RES / "03_orthologs.txt").write_text("\n".join(out) + "\n")
