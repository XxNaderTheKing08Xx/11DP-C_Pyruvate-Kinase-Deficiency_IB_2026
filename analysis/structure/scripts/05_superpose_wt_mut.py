"""Superposition of wild-type (2VGB) and R486W (2VGI) human R-type PK with Biopython Superimposer (least-squares, Kabsch).
Atoms are paired by chain ID, author residue number and atom name. Outputs data/superposition_2VGB_2VGI.json."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from Bio.PDB import Superimposer, NeighborSearch
from pkcommon import *

wt = load("2VGB"); mu = load("2VGI")
CH = ["A", "B", "C", "D"]
out = {}

def pairs(chains, sel):
    P, Q = [], []
    for c in chains:
        for r in protein_residues(wt[c]):
            n = r.id[1]
            if n not in mu[c]: continue
            r2 = mu[c][n]
            for a in r:
                nm = a.get_id()
                if sel(r, a) and nm in r2 and (r.resname == r2.resname or nm in ("N", "CA", "C", "O", "CB")):
                    P.append(a); Q.append(r2[nm])
    return P, Q

def rmsd(P, Q):
    return float(np.sqrt(np.mean(np.sum((np.array([a.coord for a in P]) - np.array([b.coord for b in Q])) ** 2, axis=1))))

def fit(P, Q):
    s = Superimposer(); s.set_atoms(P, Q)   # fixed = WT, moving = mutant
    return s

CA = lambda r, a: a.get_id() == "CA"
HEAVY = lambda r, a: a.element != "H"
# per chain
out["per_chain_CA"] = {}
for c in CH:
    P, Q = pairs([c], CA); s = fit(P, Q)
    out["per_chain_CA"][c] = {"n_CA": len(P), "rmsd_A": round(s.rms, 3)}
    # per domain (own fit)
    for dom in ("N", "A", "B", "C"):
        P2, Q2 = pairs([c], lambda r, a, dom=dom: a.get_id() == "CA" and domain_of(r.id[1]) == dom)
        s2 = fit(P2, Q2)
        out["per_chain_CA"][c][f"domain_{dom}"] = {"n_CA": len(P2), "rmsd_A": round(s2.rms, 3)}
    Ph, Qh = pairs([c], HEAVY); sh = fit(Ph, Qh)
    out["per_chain_CA"][c]["all_heavy_common_atoms"] = {"n": len(Ph), "rmsd_A": round(sh.rms, 3)}
# tetramer
P, Q = pairs(CH, CA); s = fit(P, Q)
out["tetramer_CA"] = {"n_CA": len(P), "rmsd_A": round(s.rms, 3)}
rot, tran = s.rotran
# cross-check of chain correspondence: RMSD of WT chain A vs mutant chains after tetramer fit
def apply(atoms):
    return np.array([a.coord for a in atoms]) @ rot + tran

# local environment of 486 in each chain
out["local_486"] = {}
nsw = NeighborSearch([a for a in wt.get_atoms() if a.get_parent().id[0] == " "])
for c in CH:
    r486 = wt[c][486]
    near = set()
    for a in r486:
        for b in nsw.search(a.coord, 10.0):
            near.add((b.get_parent().get_parent().id, b.get_parent().id[1]))
    near_sorted = sorted(near)
    sel = lambda chain_ids, fn: None
    def sel_local(r, a, which="CA"):
        key = (r.get_parent().id, r.id[1])
        if key not in near: return False
        if r.id[1] == 486 and r.get_parent().id == c and a.get_id() not in ("N", "CA", "C", "O", "CB"): return False
        return a.get_id() == "CA" if which == "CA" else a.element != "H"
    chains_involved = sorted({k[0] for k in near})
    Pl, Ql = pairs(chains_involved, lambda r, a: sel_local(r, a, "CA"))
    Plh, Qlh = pairs(chains_involved, lambda r, a: sel_local(r, a, "heavy"))
    # (a) after global tetramer fit
    g_ca = float(np.sqrt(np.mean(np.sum((np.array([p.coord for p in Pl]) - apply(Ql)) ** 2, axis=1))))
    g_h = float(np.sqrt(np.mean(np.sum((np.array([p.coord for p in Plh]) - apply(Qlh)) ** 2, axis=1))))
    # (b) local fit on the local CA atoms
    sl = fit(Pl, Ql); lr, lt = sl.rotran
    Xh = np.array([q.coord for q in Qlh]) @ lr + lt
    d = np.linalg.norm(np.array([p.coord for p in Plh]) - Xh, axis=1)
    per_res = {}
    for p, dd in zip(Plh, d):
        r = p.get_parent(); k = f"{r.get_parent().id}:{r.resname}{r.id[1]}"
        per_res.setdefault(k, []).append((p.get_id(), float(dd)))
    shifts = {k: {"max_shift_A": round(max(x[1] for x in v), 2), "atom": max(v, key=lambda x: x[1])[0]} for k, v in per_res.items()}
    big = {k: v for k, v in sorted(shifts.items(), key=lambda kv: -kv[1]["max_shift_A"]) if v["max_shift_A"] > 0.5}
    # 486 backbone+CB shift
    out["local_486"][c] = {"n_residues_within_10A": len(near), "chains_in_sphere": chains_involved,
                           "n_CA": len(Pl), "n_heavy": len(Plh),
                           "CA_rmsd_after_tetramer_fit": round(g_ca, 3), "heavy_rmsd_after_tetramer_fit": round(g_h, 3),
                           "CA_rmsd_local_fit": round(sl.rms, 3), "heavy_rmsd_after_local_CA_fit": round(float(np.sqrt(np.mean(d ** 2))), 3),
                           "max_atom_shift_local_fit": round(float(d.max()), 2),
                           "residue_486_backbone_CB_shift": shifts.get(f"{c}:ARG486"),
                           "residues_with_atom_shift_gt_0.5A": big,
                           "n_atoms_shift_gt_0.25A": int((d > 0.25).sum()), "n_atoms_shift_gt_0.5A": int((d > 0.5).sum())}
# side-chain chi angles of 486 and neighbours (WT vs mutant)
from Bio.PDB.vectors import calc_dihedral
def chi1(r):
    g = {"ARG": "CG", "TRP": "CG", "GLU": "CG", "ASN": "CG", "LEU": "CG", "TYR": "CG", "GLN": "CG"}.get(r.resname)
    if not g or g not in r: return None
    return round(float(np.degrees(calc_dihedral(r["N"].get_vector(), r["CA"].get_vector(), r["CB"].get_vector(), r[g].get_vector()))), 1)
out["chi1"] = {c: {n: {"WT": (wt[c][n].resname, chi1(wt[c][n])), "MUT": (mu[c][n].resname, chi1(mu[c][n]))} for n in (361, 362, 441, 453, 482, 483, 486, 487, 488)} for c in CH}
# phi/psi of 486
from Bio.PDB import PPBuilder
def phipsi(model, c, n):
    for pp in PPBuilder().build_peptides(model[c]):
        for r, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            if r.id[1] == n:
                return (round(np.degrees(phi), 1) if phi else None, round(np.degrees(psi), 1) if psi else None)
out["phi_psi_486"] = {c: {"WT": phipsi(wt, c, 486), "MUT": phipsi(mu, c, 486)} for c in CH}
out["phi_psi_364"] = {c: {"WT": phipsi(wt, c, 364), "MUT": phipsi(mu, c, 364)} for c in CH}
json.dump(out, open("data/superposition_2VGB_2VGI.json", "w"), indent=1)
print(json.dumps(out, indent=1))
