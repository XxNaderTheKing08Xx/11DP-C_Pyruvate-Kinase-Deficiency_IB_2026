"""Compose Figure 2 and Figure 3 (PNG, 300 dpi) from the PyMOL panels and the projected label positions.
Run from the analysis/structure folder after scripts/fig2_pymol.py and scripts/fig3_pymol.py:  python scripts/compose_figures.py
Input: figures/panels/ (written by the two PyMOL scripts). Output: the final figures in the repository folder figures/ (../../figures/)."""
import os
import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from PIL import Image
P = "figures/panels/"
OUT = "../../figures/"  # repository-level figures folder
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})

def crop(path, ann, margin=40):
    im = np.asarray(Image.open(path).convert("RGB"))
    mask = (im < 245).any(axis=2); ys, xs = np.where(mask)
    x0, x1 = max(xs.min() - margin, 0), min(xs.max() + margin, im.shape[1]); y0, y1 = max(ys.min() - margin, 0), min(ys.max() + margin, im.shape[0])
    return im[y0:y1, x0:x1], {k: (v[0] - x0, v[1] - y0) for k, v in ann.items()}

def lab(ax, xy, text, dx=0, dy=0, **kw):
    ax.annotate(text, xy=xy, xytext=(xy[0] + dx, xy[1] + dy), textcoords="data", ha="center", va="center", fontsize=kw.pop("fs", 8),
                fontweight=kw.pop("fw", "normal"), bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8),
                arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3") if (dx or dy) else None, **kw)

# ---------------- Figure 2 ----------------
A2 = json.load(open(P + "fig2_annotations.json"))
ima, aa = crop(P + "fig2a_tetramer.png", A2["a"]); imb, ab = crop(P + "fig2b_subunit.png", A2["b"])
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [ima.shape[1] / ima.shape[0], imb.shape[1] / imb.shape[0]]})
for ax, im in zip(axs, (ima, imb)): ax.imshow(im); ax.set_axis_off()
for ch in "ABCD":
    lab(axs[0], aa[f"subunit_{ch}"], f"Subunit {ch}", fs=8, fw="bold")
lab(axs[0], aa["R486_A"], "R486 (A)", dx=230, dy=-260, fs=7); lab(axs[0], aa["R486_C"], "R486 (C)", dx=230, dy=260, fs=7)
lab(axs[1], ab["Arg486"], "Arg486", dx=-120, dy=-170, fw="bold")
lab(axs[1], ab["FBP"], "FBP", dx=-150, dy=110)
lab(axs[1], ab["PGA"], "PGA + Mn$^{2+}$ + K$^{+}$\n(active site)", dx=80, dy=330)
lab(axs[1], ab["mid_fbp"], "12.4 Å", dx=-170, dy=0, fs=8)
lab(axs[1], ab["mid_pga"], "24.4 Å", dx=0, dy=-130, fs=8)
axs[0].text(0.0, 1.0, "a", transform=axs[0].transAxes, fontsize=12, fontweight="bold", va="top")
axs[1].text(0.0, 1.0, "b", transform=axs[1].transAxes, fontsize=12, fontweight="bold", va="top")
handles = [Patch(color=(0.90, 0.62, 0.00), label="N-terminal domain (57–84)"), Patch(color=(0.00, 0.45, 0.70), label="A domain (85–159, 263–431)"),
           Patch(color=(0.80, 0.47, 0.65), label="B domain (160–262)"), Patch(color=(0.00, 0.62, 0.45), label="C domain (432–574)"),
           Line2D([], [], marker="o", ls="", color="black", label="residue 486"), Line2D([], [], ls=":", color="black", label="minimum interatomic distance")]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.01))
plt.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig(OUT + "Figure2_PKR_overview_R486.png", dpi=300, facecolor="white")
plt.close(fig)

# ---------------- Figure 3 ----------------
A3 = json.load(open(P + "fig3_annotations.json"))
panels = [("fig3a_wt.png", "fig3a", "a  Human wild type, Arg486 (2VGB)"), ("fig3b_mut.png", "fig3b", "b  Human R486W, Trp486 (2VGI)"),
          ("fig3c_rat.png", "fig3c", "c  Rat L-type, Gln455 (6ECH)")]
dist_txt = {"fig3a": ["2.7 Å", "3.5 Å"], "fig3b": ["3.0 Å", "3.4 Å", "3.0 Å*"], "fig3c": ["2.9 Å", "2.9 Å", "2.5 Å"]}
offsets = {"Arg486": (-60, -330), "Trp486": (0, -420), "Gln455": (-80, -420), "Asn361": (110, 0), "Asn330": (110, 0),
           "Leu362": (60, 120), "Leu331": (60, 120), "Gly364": (-200, 200), "Gly333": (-60, 150), "Ala363": (70, 100), "Ala332": (70, 100),
           "Glu453": (-80, 120), "Glu422": (-80, 120), "Leu483": (-150, 0), "Leu452": (-150, 0), "Tyr487": (180, -60), "Tyr456": (180, -60),
           "water": (0, 90)}
dist_off = {("fig3a", "dist0"): (200, -120), ("fig3a", "dist1"): (180, 120), ("fig3b", "dist2"): (-60, -150),
            ("fig3b", "dist0"): (150, -60), ("fig3b", "dist1"): (170, 40)}
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.75))
for ax, (png, key, title) in zip(axs, panels):
    im = np.asarray(Image.open(P + png).convert("RGB")); y0, y1, x0, x1 = 150, 1650, 150, 1650
    ax.imshow(im[y0:y1, x0:x1]); ax.set_axis_off(); ax.set_title(title, fontsize=7, loc="left")
    ann = {k: (v[0] - x0, v[1] - y0) for k, v in A3[key].items()}
    for k, v in ann.items():
        if k.startswith("dist"):
            dx, dy = dist_off.get((key, k), (0, 0))
            lab(ax, v, dist_txt[key][int(k[4:])], dx=dx, dy=dy, fs=6.5)
        else:
            dx, dy = offsets.get(k, (0, 0))
            bold = k in ("Arg486", "Trp486", "Gln455")
            lab(ax, v, k, dx=dx, dy=dy, fs=6.5, fw="bold" if bold else "normal")
fig.text(0.01, 0.015, "Dashed lines: heavy-atom polar contacts ≤ 3.5 Å. *Carbon–oxygen contact (Trp486 CZ3–Glu453 OE1). Thin dark sticks in b and c: wild-type Arg486 after local superposition.",
         fontsize=5.8)
plt.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(OUT + "Figure3_R486W_closeup.png", dpi=300, facecolor="white")
print("done")
