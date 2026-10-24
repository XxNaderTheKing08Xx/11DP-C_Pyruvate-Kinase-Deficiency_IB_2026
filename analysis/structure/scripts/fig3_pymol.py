# Figure 3 panels (PyMOL 3.x, headless), run from the analysis/structure folder: pymol -cq scripts/fig3_pymol.py
# (a) WT Arg486 (2VGB chain A); (b) Trp486 (2VGI chain A) superposed on WT; (c) rat L-type Gln455 (6ECH chain D) superposed on WT.
# Superposition: least-squares fit (cmd.align, cycles=0, i.e. no outlier rejection) on CA atoms of the local window
# R-type residues 355-370, 436-460, 475-495 (rat L-type numbering = R-type - 31).
from pymol import cmd
import os, sys, json
BASE = os.getcwd()  # run from the analysis/structure folder (paths below are relative to it)
sys.path.insert(0, os.path.join(BASE, "scripts"))
from pymol_project import project
os.makedirs("figures/panels", exist_ok=True)
W, H = 1800, 1800
WIN_R = "resi 355-370+436-460+475-495"
WIN_L = "resi 324-339+405-429+444-464"
cmd.load("raw/2VGB.cif", "wt"); cmd.load("raw/2VGI.cif", "mut"); cmd.load("raw/6ECH.cif", "rat")
cmd.remove("hydro")
fits = {}
r = cmd.align(f"mut and chain A and name CA and {WIN_R}", f"wt and chain A and name CA and {WIN_R}", cycles=0)
fits["2VGI_A_on_2VGB_A"] = {"rmsd": round(r[0], 3), "n_CA": r[1]}
r = cmd.align(f"rat and chain D and name CA and {WIN_L}", f"wt and chain A and name CA and {WIN_R}", cycles=0)
fits["6ECH_D_on_2VGB_A"] = {"rmsd": round(r[0], 3), "n_CA": r[1]}
cmd.set_color("cWT", [0.00, 0.45, 0.70]); cmd.set_color("cMUT", [0.84, 0.37, 0.00]); cmd.set_color("cRAT", [0.00, 0.62, 0.45])
cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("orthoscopic", 1); cmd.set("ray_shadows", 0)
cmd.set("depth_cue", 0); cmd.set("antialias", 2); cmd.set("ambient", 0.5); cmd.set("specular", 0.1)
cmd.set("cartoon_transparency", 0.75); cmd.set("stick_radius", 0.22)
cmd.set("dash_color", "black"); cmd.set("dash_radius", 0.07); cmd.set("dash_gap", 0.3); cmd.set("dash_length", 0.3)
NB_R = "resi 361+362+363+364+453+483+487"
NB_L = "resi 330+331+332+333+422+452+456"
VIEW_SEL = f"wt and chain A and (resi 486 or {NB_R}) and not name N+C+O"
ANN = {"fits": fits}

def base(obj, chain, nb, res486, colour):
    cmd.hide("everything")
    cmd.show("cartoon", f"{obj} and chain {chain} and polymer.protein and {WIN_R if obj != 'rat' else WIN_L}")
    cmd.color("grey80", f"{obj} and chain {chain} and polymer.protein")
    cmd.show("sticks", f"{obj} and chain {chain} and ({nb})")
    cmd.color("grey55", f"{obj} and chain {chain} and ({nb}) and elem C")
    cmd.show("sticks", f"{obj} and chain {chain} and resi {res486} and not name N+C+O")
    cmd.color(colour, f"{obj} and chain {chain} and resi {res486} and elem C")
    cmd.color("red", f"{obj} and chain {chain} and elem O"); cmd.color("blue", f"{obj} and chain {chain} and elem N")
    cmd.set("stick_radius", 0.30, f"{obj} and chain {chain} and resi {res486}")

def ghost_wt():
    cmd.show("sticks", "wt and chain A and resi 486 and not name N+C+O"); cmd.color("grey20", "wt and chain A and resi 486")
    cmd.set("stick_radius", 0.12, "wt and chain A and resi 486")

def render(name, dists, labels):
    cmd.orient(VIEW_SEL)
    cmd.set_view(VIEW)
    for i, (a1, a2) in enumerate(dists):
        cmd.distance(f"{name}_d{i}", a1, a2); cmd.hide("labels", f"{name}_d{i}")
    ann = {k: project(v, W, H) for k, v in labels.items()}
    for i, (a1, a2) in enumerate(dists):
        ann[f"dist{i}"] = project(f"({a1}) or ({a2})", W, H)
    cmd.png(f"figures/panels/{name}.png", width=W, height=H, dpi=300, ray=1)
    cmd.delete(f"{name}_d*")
    return ann

# common view, defined once on the WT environment
cmd.orient(VIEW_SEL); cmd.zoom(VIEW_SEL, 3.0)
cmd.turn("y", -25)
VIEW = list(cmd.get_view()); VIEW[15] -= 40; VIEW[16] += 40   # widen the clipping slab

base("wt", "A", NB_R, 486, "cWT")
ANN["fig3a"] = render("fig3a_wt", [("wt and chain A and resi 486 and name NH2", "wt and chain A and resi 362 and name O"),
                                   ("wt and chain A and resi 486 and name NH1", "wt and chain A and resi 362 and name O")],
                      {"Arg486": "wt and chain A and resi 486 and name CZ", "Asn361": "wt and chain A and resi 361 and name CG",
                       "Leu362": "wt and chain A and resi 362 and name CG", "Gly364": "wt and chain A and resi 364 and name CA",
                       "Ala363": "wt and chain A and resi 363 and name CB",
                       "Glu453": "wt and chain A and resi 453 and name CD", "Leu483": "wt and chain A and resi 483 and name CG",
                       "Tyr487": "wt and chain A and resi 487 and name CZ"})
base("mut", "A", NB_R, 486, "cMUT"); ghost_wt()
ANN["fig3b"] = render("fig3b_mut", [("mut and chain A and resi 486 and name NE1", "mut and chain A and resi 361 and name O"),
                                    ("mut and chain A and resi 486 and name NE1", "mut and chain A and resi 362 and name O"),
                                    ("mut and chain A and resi 486 and name CZ3", "mut and chain A and resi 453 and name OE1")],
                      {"Trp486": "mut and chain A and resi 486 and name CZ2", "Asn361": "mut and chain A and resi 361 and name CG",
                       "Leu362": "mut and chain A and resi 362 and name CG", "Gly364": "mut and chain A and resi 364 and name CA",
                       "Glu453": "mut and chain A and resi 453 and name CD", "Leu483": "mut and chain A and resi 483 and name CG",
                       "Tyr487": "mut and chain A and resi 487 and name CZ", "Ala363": "mut and chain A and resi 363 and name CB"})
base("rat", "D", NB_L, 455, "cRAT"); ghost_wt()
cmd.show("nb_spheres", "rat and chain D and resn HOH and resi 717"); cmd.set("nb_spheres_size", 0.35); cmd.color("red", "rat and resn HOH")
ANN["fig3c"] = render("fig3c_rat", [("rat and chain D and resi 455 and name NE2", "rat and chain D and resi 422 and name OE2"),
                                    ("rat and chain D and resi 455 and name NE2", "rat and chain D and resi 452 and name O"),
                                    ("rat and chain D and resi 455 and name OE1", "rat and chain D and resn HOH and resi 717")],
                      {"Gln455": "rat and chain D and resi 455 and name CD", "Asn330": "rat and chain D and resi 330 and name CG",
                       "Leu331": "rat and chain D and resi 331 and name CG", "Gly333": "rat and chain D and resi 333 and name CA",
                       "Glu422": "rat and chain D and resi 422 and name CD", "Leu452": "rat and chain D and resi 452 and name CG",
                       "Tyr456": "rat and chain D and resi 456 and name CZ", "Ala332": "rat and chain D and resi 332 and name CB",
                       "water": "rat and chain D and resn HOH and resi 717"})
# measured distances (A) printed for the record
for obj, pairs in [("wt", [("486/NH2", "362/O"), ("486/NH1", "362/O")]), ("mut", [("486/NE1", "361/O"), ("486/NE1", "362/O"), ("486/CZ3", "453/OE1")])]:
    for a, b in pairs:
        ANN.setdefault("dist_values", {})[f"{obj} A {a}-{b}"] = round(cmd.get_distance(f"{obj} and chain A and resi {a.split('/')[0]} and name {a.split('/')[1]}", f"{obj} and chain A and resi {b.split('/')[0]} and name {b.split('/')[1]}"), 2)
for a, b in [("455/NE2", "422/OE2"), ("455/NE2", "452/O")]:
    ANN["dist_values"][f"rat D {a}-{b}"] = round(cmd.get_distance(f"rat and chain D and resi {a.split('/')[0]} and name {a.split('/')[1]}", f"rat and chain D and resi {b.split('/')[0]} and name {b.split('/')[1]}"), 2)
ANN["dist_values"]["rat D 455/OE1-HOH717"] = round(cmd.get_distance("rat and chain D and resi 455 and name OE1", "rat and chain D and resn HOH and resi 717"), 2)
json.dump(ANN, open("figures/panels/fig3_annotations.json", "w"), indent=1)
print(json.dumps(ANN["fits"]), json.dumps(ANN["dist_values"]))
