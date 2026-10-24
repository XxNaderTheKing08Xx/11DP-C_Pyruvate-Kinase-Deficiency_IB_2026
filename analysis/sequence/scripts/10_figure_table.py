"""Task 3.2: Figure 1 (MSA window around human PKLR Arg486) and Table 2 (pairwise identities and
residue at 486 per species).

Input: the broad-panel Clustal Omega alignment (EBI output if present in raw/clustalo/broad_panel,
otherwise results/msa_local/broad_panel.clustalo.aln), results/panel_metadata.tsv,
results/pairwise_summary.tsv, results/09_verify.txt (codons), raw CDS files.
Output: ../../figures/Figure1_PKLR_MSA_486.png (300 dpi) and .svg (repository-level figures folder),
        results/Table2_pairwise_and_486.md, results/figure1_window_data.tsv
Run: python scripts/10_figure_table.py
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import AlignIO, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

B = Path(__file__).resolve().parent.parent
RAW, RES, FIG = B / "raw", B / "results", B.parent.parent / "figures"  # repository-level figures/
FIG.mkdir(exist_ok=True)
HPOS, WIN = 486, 10

ebi = sorted((RAW / "clustalo/broad_panel").glob("*.aln-clustal_num.txt")) if (RAW / "clustalo/broad_panel").exists() else []
alnfile = ebi[-1] if ebi else RES / "msa_local/broad_panel.clustalo.aln"
aln = AlignIO.read(alnfile, "clustal")
meta = {}
for line in (RES / "panel_metadata.tsv").read_text().splitlines()[1:]:
    f = line.split("\t")
    meta[f"{f[0].replace(' ', '_')}_{f[4]}"] = dict(label=f[0], species=f[1], clade=f[2], acc=f[4],
                                                    length=int(f[9]), id_pklr=float(f[10]))
order = ["Human", "Rat", "Mouse", "Dog", "Chimpanzee", "Macaque", "Marmoset", "Guinea pig", "Rabbit",
         "Cat", "Cow", "Pig", "Horse", "Opossum", "Platypus", "Painted turtle", "Wall lizard", "Xenopus",
         "Coelacanth", "Zebrafish", "Medaka"]
recs = {meta[r.id]["label"]: (r.id, str(r.seq)) for r in aln}
hseq = recs["Human"][1]
col_of, n = {}, 0
for c, ch in enumerate(hseq):
    if ch != "-":
        n += 1
        col_of[n] = c
c0, c1 = col_of[HPOS - WIN], col_of[HPOS + WIN]
cols = list(range(c0, c1 + 1))
assert len(cols) == 2 * WIN + 1, "insertion columns inside the window; adapt the figure"
c486 = col_of[HPOS]


def ownpos(seq, col):
    return sum(1 for ch in seq[:col + 1] if ch != "-")


# Clustal conservation symbols (strong / weak groups as in Clustal Omega/ClustalW)
strong = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
weak = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY"]


def symbol(colres):
    s = set(colres)
    if "-" in s:
        return " "
    if len(s) == 1:
        return "*"
    if any(s <= set(g) for g in strong):
        return ":"
    if any(s <= set(g) for g in weak):
        return "."
    return " "


rows = [(lab, recs[lab][1]) for lab in order]
colres = {c: [s[c] for _, s in rows] for c in cols}
majority = {c: max(set(colres[c]), key=colres[c].count) for c in cols}

# ---- figure
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
ink, muted, grid = "#1f1f1e", "#5f5e5a", "#e6e5e0"
R_FILL, Q_FILL = "#9ec5f4", "#f6b89c"   # tints of categorical blue #2a78d6 / orange #eb6834
cw, rh = 0.36, 0.30            # cell width/height in inches
left_w, num_w, right_w = 2.75, 0.42, 0.42
nrows = len(rows) + 1          # + gap row between core and broad blocks
fig_w = left_w + num_w + len(cols) * cw + right_w + 0.15
fig_h = 0.75 + (nrows + 1 + 2.2) * rh + 0.45
fig = plt.figure(figsize=(fig_w, fig_h))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, fig_w)
ax.set_ylim(fig_h, 0)
ax.axis("off")
x0 = left_w + num_w
y = 0.62

# ruler (human R-type numbering)
ax.text(x0 - 0.08, y - 0.02, "Human R-type position", ha="right", va="center", fontsize=7.5, color=muted)
for k, c in enumerate(cols):
    p = HPOS - WIN + k
    if p % 5 == 1 or p == HPOS:
        ax.text(x0 + (k + 0.5) * cw, y - 0.02, str(p), ha="center", va="center", fontsize=7.5,
                color=ink if p == HPOS else muted, fontweight="bold" if p == HPOS else "normal")
k486 = cols.index(c486)
ax.text(x0 + (k486 + 0.5) * cw, y - 0.30, "human p.Arg486Trp (c.1456C>T)", fontsize=7.5, color=ink,
        ha="center", va="center")
y += 0.22
top = y + rh
ax.text(0.12, y + rh * 0.55, "Core panel (UniProtKB/Swiss-Prot, R-type)", fontsize=7.5, color=muted,
        va="center", style="italic")
y += rh
data_rows = []
for i, (lab, s) in enumerate(rows):
    if i == 4:
        ax.text(0.12, y + rh * 0.55, "Additional Ensembl Compara one-to-one orthologues", fontsize=7.5,
                color=muted, va="center", style="italic")
        y += rh
    m = meta[recs[lab][0]]
    first = ownpos(s, c0) if s[c0] != "-" else ownpos(s, c0) + 1
    last = ownpos(s, c1)
    ax.text(0.12, y + rh / 2, f"{lab} ", fontsize=8, color=ink, va="center", ha="left")
    ax.text(1.12, y + rh / 2, m["species"], fontsize=7.5, color=ink, va="center", ha="left", style="italic")
    ax.text(x0 - 0.06, y + rh / 2, str(first), fontsize=7, color=muted, va="center", ha="right")
    for k, c in enumerate(cols):
        aa = s[c]
        xx = x0 + k * cw
        if c == c486:
            fill = R_FILL if aa == "R" else Q_FILL if aa == "Q" else "#ffffff"
        else:
            fill = grid if aa == majority[c] else "#ffffff"
        ax.add_patch(Rectangle((xx + 0.01, y + 0.01), cw - 0.02, rh - 0.02, facecolor=fill, edgecolor="none"))
        ax.text(xx + cw / 2, y + rh / 2, aa, ha="center", va="center", fontsize=8.5, family="DejaVu Sans Mono",
                color=ink, fontweight="bold" if c == c486 else "normal")
    ax.text(x0 + len(cols) * cw + 0.06, y + rh / 2, str(last), fontsize=7, color=muted, va="center", ha="left")
    own486 = ownpos(s, c486)
    data_rows.append((lab, m["species"], m["acc"], m["source"] if "source" in m else "", first, last,
                      s[c486], own486, "".join(s[c] for c in cols)))
    y += rh
# conservation row
y += 0.06
ax.text(x0 - 0.06, y + rh / 2, "Conservation", fontsize=7.5, color=muted, va="center", ha="right")
for k, c in enumerate(cols):
    ax.text(x0 + (k + 0.5) * cw, y + rh / 2, symbol(colres[c]), ha="center", va="center", fontsize=9,
            family="DejaVu Sans Mono", color=ink)
y += rh
ax.text(x0 - 0.06, y + rh / 2, "Majority-residue fraction", fontsize=7.5, color=muted, va="center", ha="right")
for k, c in enumerate(cols):
    frac = colres[c].count(majority[c]) / len(rows)
    ax.add_patch(Rectangle((x0 + k * cw + 0.06, y + rh - 0.04 - frac * (rh - 0.06)), cw - 0.12,
                           frac * (rh - 0.06), facecolor="#5f5e5a", edgecolor="none"))
y += rh
# column 486 frame
ax.add_patch(Rectangle((x0 + k486 * cw, top), cw, y - top - 2 * rh - 0.06, fill=False, edgecolor=ink, lw=1.1))
# legend
y += 0.12
lx = x0
for fill, txt in [(R_FILL, "Arg (R) at human 486 column"), (Q_FILL, "Gln (Q) at human 486 column"),
                  (grid, "residue = column majority")]:
    ax.add_patch(Rectangle((lx, y), 0.2, 0.16, facecolor=fill, edgecolor="none"))
    ax.text(lx + 0.27, y + 0.08, txt, fontsize=7, color=ink, va="center")
    lx += 2.25 if fill != Q_FILL else 2.25
ax.text(0.12, y + 0.08, "* identical  : strongly similar  . weakly similar", fontsize=7, color=muted, va="center")
fig.savefig(FIG / "Figure1_PKLR_MSA_486.png", dpi=300)
fig.savefig(FIG / "Figure1_PKLR_MSA_486.svg")
print("saved figure from", alnfile)

with open(RES / "figure1_window_data.tsv", "w") as fh:
    fh.write("label\tspecies\taccession\tsource\twindow_first_own\twindow_last_own\tres486\town_pos486\twindow\n")
    for r in data_rows:
        fh.write("\t".join(map(str, r)) + "\n")
syms = "".join(symbol(colres[c]) for c in cols)
print("conservation row:", repr(syms))

# ---- Table 2
pw = {}
for line in (RES / "pairwise_summary.tsv").read_text().splitlines()[1:]:
    f = line.split("\t")
    pw[f[0].split(" vs ")[1].split()[0]] = f
codon = {}
for line in (RES / "09_verify.txt").read_text().splitlines():
    if "Swiss-Prot residue 486" in line:
        sp = line.split()[0]
        codon[sp] = line.split("= ")[-1].strip()
cds_acc = {"Human": "NM_000298.6", "Rat": "XM_006232593.5", "Mouse": "NM_013631.2", "Dog": "NM_001256262.2"}
prot_acc = {"Human": "NP_000289.1", "Rat": "XP_006232655.1 (L-type NP_036756.3)", "Mouse": "NP_038659.2",
            "Dog": "NP_001243191.1"}
manual_cds = {}
for line in (RES / "09_verify.txt").read_text().splitlines():
    if line.startswith("Human vs"):
        sp = line.split()[2].rstrip(":")
        manual_cds[sp] = line.split("CDS ")[1].split(" identical")[0] + " (" + line.split("identical (")[2].split(")")[0] + ")"
L = []
L.append("**Table 2.** Pairwise comparison of human PKLR (R-type) with orthologues, and the residue at the position "
         "equivalent to human Arg486. Protein: global alignment (Needleman–Wunsch; BLOSUM62, gap open 10, gap "
         "extension 0.5, end gaps not penalized), identity and similarity as fraction of alignment columns. CDS: "
         "position-by-position identity of the 1725-nt R-type coding sequences (no indels at protein level). "
         "Residue numbering is that of each species' own sequence. Lower block: broad panel (Ensembl release 116 "
         "Compara one-to-one orthologues), used only for the conservation of position 486.\n")
L.append("| Species | Protein (UniProtKB / RefSeq) | CDS (RefSeq) | Protein identity to human | Protein similarity | CDS identity | Residue at human 486 (codon) |")
L.append("|---|---|---|---|---|---|---|")
L.append(f"| Human (*Homo sapiens*) | P30613 / NP_000289.1 | NM_000298.6 | — | — | — | Arg486 (CGG) |")
for sp, lat in [("Rat", "Rattus norvegicus"), ("Mouse", "Mus musculus"), ("Dog", "Canis lupus familiaris")]:
    f = pw[sp]
    aa = {"R": "Arg", "Q": "Gln"}[next(r[6] for r in data_rows if r[0] == sp)]
    L.append(f"| {sp} (*{lat}*) | {meta[next(k for k in meta if meta[k]['label']==sp)]['acc']} / {prot_acc[sp]} | {cds_acc[sp]} | "
             f"{f[2]}/{f[1]} ({f[3]}%) | {f[4]}/{f[1]} ({f[5]}%) | {manual_cds[sp]} | {aa}486 ({codon[sp]}) |")
L.append("")
L.append("| Species (broad panel) | Ensembl protein | Length (aa) | Protein identity to human (%) | Residue at human 486 (own number) |")
L.append("|---|---|---|---|---|")
three = {"R": "Arg", "Q": "Gln", "K": "Lys"}
for r in data_rows[4:]:
    m = meta[recs[r[0]][0]]
    L.append(f"| {r[0]} (*{r[1]}*) | {r[2]} | {m['length']} | {m['id_pklr']:.1f} | {three.get(r[6], r[6])}{r[7]} |")
(RES / "Table2_pairwise_and_486.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
