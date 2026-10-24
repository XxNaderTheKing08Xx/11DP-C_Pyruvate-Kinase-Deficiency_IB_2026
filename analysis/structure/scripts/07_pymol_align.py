# Run from the analysis/structure folder: pymol -cq scripts/07_pymol_align.py
# PyMOL align/super comparison for 2VGB (WT) vs 2VGI (R486W). Results to data/pymol_align.tsv
from pymol import cmd
import os
BASE = os.getcwd()  # run from the analysis/structure folder (paths below are relative to it)
cmd.load("raw/2VGB.cif", "wt"); cmd.load("raw/2VGI.cif", "mut")
cmd.remove("hydro")
rows = []
for scope, sel in [("chain A", "chain A"), ("chain B", "chain B"), ("chain C", "chain C"), ("chain D", "chain D"), ("tetramer", "polymer.protein")]:
    for method in ("align", "super"):
        for cycles in (0, 5):
            f = getattr(cmd, method)
            r = f(f"mut and {sel} and name CA", f"wt and {sel} and name CA", cycles=cycles)
            rows.append((scope, method, cycles, round(r[0], 3), r[1], round(r[3], 3), r[4]))
# local: residues within 10 A of residue 486 (chain A), CA atoms, fitted locally
cmd.select("locA", "byres (wt and polymer.protein within 10 of (wt and chain A and resi 486))")
cmd.select("locA_mut", "mut and polymer.protein and (chain+resi) in locA" if False else "none")
for method in ("align", "super"):
    r = getattr(cmd, method)("mut and name CA and (byres (mut and polymer.protein within 10 of (mut and chain A and resi 486)))",
                             "wt and name CA and locA", cycles=0)
    rows.append(("local 10 A sphere around A486 (CA)", method, 0, round(r[0], 3), r[1], round(r[3], 3), r[4]))
with open("data/pymol_align.tsv", "w") as fh:
    fh.write("scope\tmethod\tcycles\trmsd_after_refinement_A\tn_atoms_after\trmsd_before_refinement_A\tn_atoms_before\n")
    for x in rows:
        fh.write("\t".join(map(str, x)) + "\n")
for x in rows: print(x)
