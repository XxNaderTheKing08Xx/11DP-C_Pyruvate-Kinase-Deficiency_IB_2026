# Run from the analysis/structure folder: pymol -cq scripts/08_secondary_structure.py
# Requires raw/8TBS.cif (not distributed; download command in README.md).
# Secondary structure assigned by PyMOL `dss` (PyMOL's own algorithm, not DSSP) for residues 470-495 (human R-type numbering)
# and the equivalent rat L-type window; output data/ss_pymol.tsv
from pymol import cmd, stored
import os
BASE = os.getcwd()  # run from the analysis/structure folder (paths below are relative to it)
out = []
for pid, rng in [("2VGB", "470-495"), ("2VGI", "470-495"), ("8TBS", "470-495"), ("6ECH", "439-464")]:
    cmd.load(f"raw/{pid}.cif", pid)
    cmd.dss(pid)
    for ch in ("A", "B", "C", "D"):
        stored.l = []
        cmd.iterate(f"{pid} and chain {ch} and resi {rng} and name CA", "stored.l.append((resi, resn, ss))")
        s = "".join(x[2] if x[2] else "-" for x in stored.l)
        out.append((pid, ch, rng, s, " ".join(f"{x[1]}{x[0]}" for x in stored.l)))
    cmd.delete(pid)
with open("data/ss_pymol.tsv", "w") as f:
    f.write("pdb\tchain\trange\tss(H=helix,S=strand,-/L=loop)\tresidues\n")
    for x in out: f.write("\t".join(x) + "\n")
for x in out: print(x[0], x[1], x[2], x[3])
