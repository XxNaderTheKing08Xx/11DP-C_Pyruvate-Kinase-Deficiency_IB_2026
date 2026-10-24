# Figure 2 panels (PyMOL 3.x, headless), run from the analysis/structure folder:  pymol -cq scripts/fig2_pymol.py
# Human R-type PK, PDB 2VGB (wild type). Domains after Valentini et al. 2002: N 57-84, A 85-159 + 263-431, B 160-262, C 432-574.
from pymol import cmd
import os, sys, json
sys.path.insert(0, os.path.join(os.getcwd(), "scripts"))  # run from the analysis/structure folder
from pymol_project import project
W, H = 2400, 2000
ANN = {}
BASE = os.getcwd()
os.makedirs("figures/panels", exist_ok=True)
# Okabe-Ito colour-blind-safe palette
def common():
    cmd.set_color("dN", [0.90, 0.62, 0.00])   # orange
    cmd.set_color("dA", [0.00, 0.45, 0.70])   # blue
    cmd.set_color("dB", [0.80, 0.47, 0.65])   # reddish purple
    cmd.set_color("dC", [0.00, 0.62, 0.45])   # bluish green
    cmd.set_color("mut", [0.0, 0.0, 0.0])  # black for residue 486
    cmd.bg_color("white"); cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_fancy_helices", 1); cmd.set("cartoon_highlight_color", "grey70")
    cmd.set("ray_shadows", 0); cmd.set("ambient", 0.45); cmd.set("specular", 0.15)
    cmd.set("depth_cue", 0); cmd.set("antialias", 2); cmd.set("orthoscopic", 1)
    cmd.set("label_font_id", 7); cmd.set("label_color", "black"); cmd.set("float_labels", 1)
    cmd.set("dash_color", "black"); cmd.set("dash_gap", 0.35); cmd.set("dash_radius", 0.12)
    cmd.set("label_distance_digits", 1)
def colour_domains(obj):
    cmd.color("dN", f"{obj} and polymer.protein and resi 57-84")
    cmd.color("dA", f"{obj} and polymer.protein and (resi 85-159 or resi 263-431)")
    cmd.color("dB", f"{obj} and polymer.protein and resi 160-262")
    cmd.color("dC", f"{obj} and polymer.protein and resi 432-574")

# ---------- panel a: tetramer ----------
cmd.reinitialize(); common()
cmd.load("raw/2VGB.cif", "pk"); cmd.remove("solvent or hydro")
cmd.hide("everything"); cmd.show("cartoon", "pk and polymer.protein"); colour_domains("pk")
cmd.show("spheres", "pk and resi 486 and not name N+C+O"); cmd.color("mut", "pk and resi 486"); cmd.set("sphere_scale", 1.2, "pk and resi 486")
cmd.show("sticks", "pk and resn FBP+PGA"); cmd.set("stick_radius", 0.3)
cmd.color("grey30", "pk and resn FBP+PGA and elem C")
cmd.show("spheres", "pk and resn MN+K"); cmd.set("sphere_scale", 0.6, "pk and resn MN+K")
cmd.color("purple", "pk and resn MN"); cmd.color("grey50", "pk and resn K")
cmd.orient("pk and polymer.protein")
cmd.turn("x", 90)   # look down one molecular 2-fold axis so that the four subunits are side by side
cmd.zoom("pk and polymer.protein", 2)
ANN["a"] = {f"subunit_{ch}": project(f"pk and chain {ch} and polymer.protein and resi 160-262 and name CA", W, H) for ch in "ABCD"}
ANN["a"].update({f"R486_{ch}": project(f"pk and chain {ch} and resi 486 and name CZ", W, H) for ch in "ABCD"})
cmd.png("figures/panels/fig2a_tetramer.png", width=W, height=H, dpi=300, ray=1)
cmd.save("figures/panels/fig2a_tetramer.pse")

# ---------- panel b: one subunit with ligands and distances ----------
cmd.reinitialize(); common()
cmd.load("raw/2VGB.cif", "pk"); cmd.remove("solvent or hydro"); cmd.remove("pk and not chain A")
cmd.hide("everything"); cmd.show("cartoon", "pk and polymer.protein"); colour_domains("pk")
cmd.set("cartoon_transparency", 0.15)
cmd.show("sticks", "pk and resi 486 and not name N+C+O"); cmd.color("mut", "pk and resi 486 and elem C")
cmd.show("spheres", "pk and resi 486 and name CA"); cmd.set("sphere_scale", 0.5, "pk and resi 486")
cmd.show("sticks", "pk and resn FBP+PGA"); cmd.color("grey30", "pk and resn FBP+PGA and elem C")
cmd.show("spheres", "pk and resn MN+K"); cmd.set("sphere_scale", 0.6, "pk and resn MN+K")
cmd.color("purple", "pk and resn MN"); cmd.color("grey50", "pk and resn K")
cmd.set("stick_radius", 0.28)
cmd.util.cnc("pk and (resn FBP+PGA or resi 486)")
# minimum atom-atom distances (values from scripts/03_site_analysis.py)
cmd.distance("d_fbp", "pk and resi 486 and name CB", "pk and resn FBP and name O6P")
cmd.distance("d_pga", "pk and resi 486 and name CA", "pk and resn PGA and name C2")
cmd.hide("labels", "d_*")
for obj in ("d_fbp", "d_pga"):
    cmd.color("black", obj)
cmd.orient("pk and polymer.protein")
cmd.zoom("pk and polymer.protein", 1)
ANN["b"] = {"Arg486": project("pk and resi 486 and name CZ", W, H), "FBP": project("pk and resn FBP", W, H),
            "PGA": project("pk and resn PGA", W, H), "Mn": project("pk and resn MN", W, H), "K": project("pk and resn K", W, H),
            "mid_fbp": project("(pk and resi 486 and name CB) or (pk and resn FBP and name O6P)", W, H),
            "mid_pga": project("(pk and resi 486 and name CA) or (pk and resn PGA and name C2)", W, H)}
cmd.png("figures/panels/fig2b_subunit.png", width=W, height=H, dpi=300, ray=1)
cmd.save("figures/panels/fig2b_subunit.pse")

json.dump(ANN, open("figures/panels/fig2_annotations.json", "w"), indent=1)
