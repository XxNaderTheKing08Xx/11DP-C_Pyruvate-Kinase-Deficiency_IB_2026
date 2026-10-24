"""Shared helpers for the PKLR structure analysis (Biopython)."""
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

# Domain boundaries of human R-type PK (R-type numbering) as reported by Valentini et al. 2002, J Biol Chem 277:23807, p. 23809
DOMAINS_R = [("N", 57, 84), ("A", 85, 159), ("B", 160, 262), ("A", 263, 431), ("C", 432, 574)]

def domain_of(resnum, offset=0):
    """Domain of a residue; offset converts the structure's numbering to R-type numbering (R = n + offset)."""
    r = resnum + offset
    for name, a, b in DOMAINS_R:
        if a <= r <= b:
            return name
    return "outside"

def load(pid, path=None):
    return MMCIFParser(QUIET=True).get_structure(pid, path or f"raw/{pid}.cif")[0]

def protein_residues(chain):
    return [r for r in chain if r.id[0] == " " or r.resname in ("SEP", "TPO")]

def hetero(model, resname):
    return [r for ch in model for r in ch if r.resname == resname]

def min_dist(res_a, res_b, atoms_a=None):
    A = [a for a in res_a if (atoms_a is None or a.get_id() in atoms_a)]
    best = (1e9, None, None)
    for a in A:
        for b in res_b:
            d = a - b
            if d < best[0]:
                best = (d, a, b)
    return best

def fmt_atom(a):
    r = a.get_parent(); c = r.get_parent()
    return f"{c.id}:{r.resname}{r.id[1]}{r.id[2].strip()}:{a.get_id()}"

ARG_DONORS = {"NE", "NH1", "NH2"}
TRP_DONORS = {"NE1"}
ACIDIC = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
AROMATIC_RINGS = {"PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                  "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                  "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
                  "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}

def polar_contacts(res, model, donors, cutoff=3.5):
    """Side-chain donor atoms of `res` to any N/O atom (protein, ligand, water) within cutoff (heavy-atom distance)."""
    ns = NeighborSearch([a for a in model.get_atoms()])
    out = []
    for d in res:
        if d.get_id() not in donors:
            continue
        for b in ns.search(d.coord, cutoff):
            if b.get_parent() is res:
                continue
            if b.element not in ("N", "O"):
                continue
            out.append((d.get_id(), fmt_atom(b), round(float(d - b), 2)))
    return sorted(out, key=lambda x: x[2])

def salt_bridges(res, model, cutoff=4.0):
    ns = NeighborSearch([a for a in model.get_atoms()])
    out = []
    for d in res:
        if d.get_id() not in ARG_DONORS:
            continue
        for b in ns.search(d.coord, cutoff):
            p = b.get_parent()
            if (p.resname, b.get_id()) in ACIDIC or (p.id[0].startswith("H_") and b.element == "O" and p.resname in ("FBP", "PGA", "ATP", "OXL", "FLC")):
                out.append((d.get_id(), fmt_atom(b), round(float(d - b), 2)))
    return sorted(out, key=lambda x: x[2])

def cation_pi(res, model, cutoff=6.0):
    """Arg CZ to aromatic ring centroid distance (cation-pi candidate if <= 6.0 A)."""
    if "CZ" not in res:
        return []
    cz = res["CZ"].coord
    out = []
    for ch in model:
        for r in ch:
            if r.resname in AROMATIC_RINGS and r is not res:
                try:
                    cen = np.mean([r[n].coord for n in AROMATIC_RINGS[r.resname]], axis=0)
                except KeyError:
                    continue
                d = float(np.linalg.norm(cz - cen))
                if d <= cutoff:
                    out.append((f"{ch.id}:{r.resname}{r.id[1]}", round(d, 2)))
    return out

def all_contacts(res, model, cutoff=4.0, sidechain_only=False):
    ns = NeighborSearch([a for a in model.get_atoms()])
    out = {}
    bb = {"N", "CA", "C", "O"}
    for a in res:
        if sidechain_only and a.get_id() in bb:
            continue
        for b in ns.search(a.coord, cutoff):
            p = b.get_parent()
            if p is res:
                continue
            # skip covalent neighbours (peptide bond partners)
            key = f"{p.get_parent().id}:{p.resname}{p.id[1]}"
            d = float(a - b)
            if key not in out or d < out[key][0]:
                out[key] = (round(d, 2), a.get_id(), b.get_id())
    return dict(sorted(out.items(), key=lambda kv: kv[1][0]))
