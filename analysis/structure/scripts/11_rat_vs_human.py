"""Superpose rat L-type PK (6ECH, author numbering = rat L-type numbering) onto human R-type PK (2VGB, R-type numbering).
Residue correspondence: rat L n <-> human R n+31 (verified gapless by scripts/10_rat_mapping.py: same length in the aligned region).
Outputs data/rat_vs_human.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from Bio.PDB import Superimposer, NeighborSearch
from pkcommon import *
H = load("2VGB"); R = load("6ECH"); RK = load("6ECK")
OFF = 31
out = {}
def ca_pairs(hc, rc, model_r, keep=lambda n: True):
    P, Q, same = [], [], 0
    for r in protein_residues(H[hc]):
        n = r.id[1]; m = n - OFF
        if m in model_r[rc] and keep(n) and "CA" in r and "CA" in model_r[rc][m]:
            P.append(r["CA"]); Q.append(model_r[rc][m]["CA"]); same += r.resname == model_r[rc][m].resname
    return P, Q, same
for rc in "ABCD":
    P, Q, same = ca_pairs("A", rc, R)
    s = Superimposer(); s.set_atoms(P, Q)
    d = {"n_CA": len(P), "identical_resnames": same, "rmsd_A": round(s.rms, 3)}
    P3, Q3, _ = ca_pairs("A", rc, R, lambda n: domain_of(n) in "NAC")
    s3 = Superimposer(); s3.set_atoms(P3, Q3); d["core_N+A+C"] = {"n": len(P3), "rmsd_A": round(s3.rms, 3)}
    for dom in "NABC":
        P2, Q2, _ = ca_pairs("A", rc, R, lambda n, dom=dom: domain_of(n) == dom)
        s2 = Superimposer(); s2.set_atoms(P2, Q2); d[f"domain_{dom}"] = {"n": len(P2), "rmsd_A": round(s2.rms, 3)}
    out[f"2VGB_A_vs_6ECH_{rc}_CA"] = d
# tetramer (chain correspondence chosen by matching interface types: human A-B C/C', A-C A/A'; rat A-C C/C', A-B A/A')
def tet_pairs(mapping):
    P, Q = [], []
    for hc, rc in mapping.items():
        p, q, _ = ca_pairs(hc, rc, R); P += p; Q += q
    return P, Q
best = None
import itertools
for perm in itertools.permutations("ABCD"):
    mp = dict(zip("ABCD", perm))
    P, Q = tet_pairs(mp); s = Superimposer(); s.set_atoms(P, Q)
    if best is None or s.rms < best[0]: best = (s.rms, mp, len(P))
out["tetramer_best_chain_mapping"] = {"rmsd_A": round(best[0], 3), "mapping_human_to_rat": best[1], "n_CA": best[2]}

# local environment: 10 A sphere around human A:486 (any atom), CA local fit, then key distances in both
ns = NeighborSearch([a for a in H.get_atoms() if a.get_parent().id[0] == " "])
near = sorted({(b.get_parent().get_parent().id, b.get_parent().id[1]) for a in H["A"][486] for b in ns.search(a.coord, 10.0)})
mp = best[1]
P, Q, diffs = [], [], []
for hc, n in near:
    rc = mp[hc]; m = n - OFF
    if m in R[rc] and "CA" in R[rc][m]:
        P.append(H[hc][n]["CA"]); Q.append(R[rc][m]["CA"])
        if H[hc][n].resname != R[rc][m].resname: diffs.append(f"{hc}:{H[hc][n].resname}{n} -> rat {rc}:{R[rc][m].resname}{m}")
s = Superimposer(); s.set_atoms(P, Q)
out["local_sphere_486"] = {"n_CA": len(P), "CA_rmsd_local_fit_A": round(s.rms, 3), "residue_differences_in_sphere": diffs}
rot, tran = s.rotran
def xyz(atom, rat=False):
    return atom.coord @ rot + tran if rat else atom.coord
# key distances
def dist(a, b): return round(float(np.linalg.norm(a - b)), 2)
key = {}
RC = mp["A"]  # rat chain occupying the position of human chain A in the superposed tetramer
for label, model, ch, off, rat in [("human 2VGB A (Arg486)", H, "A", 0, False)] + [(f"rat 6ECH {x} (Gln455)", R, x, OFF, True) for x in "ABCD"]:
    c = model[ch]; r = c[486 - off]
    tips = [a for a in r if a.get_id() in ("NE", "NH1", "NH2", "OE1", "NE2")]
    k = {}
    for nm, (n, atoms) in {"L362 O / rat L331 O": (362, ["O"]), "N361 O / rat N330 O": (361, ["O"]), "G364 N / rat G333 N": (364, ["N"]),
                           "E453 carboxylate / rat E422": (453, ["OE1", "OE2"]), "L483 O / rat L452 O": (483, ["O"]),
                           "Y487 OH / rat Y456 OH": (487, ["OH"])}.items():
        rr = c[n - off]
        best_ = min(((t, rr[x]) for t in tips for x in atoms if x in rr), key=lambda p: p[0] - p[1])
        k[nm] = (best_[0].get_id(), best_[1].get_id(), round(float(best_[0] - best_[1]), 2))
    k["CA(486/455) to FBP (same chain, min any atom)"] = round(float(min(r["CA"] - b for x in hetero(model, "FBP") if x.get_parent().id == ch for b in x)), 2)
    key[label] = k
# displacement of the side-chain tip and CA after local fit
h486 = H["A"][486]; q455 = R[RC][455]
key["rat chain superposed on human chain A"] = RC
key["CA displacement human486 vs rat455 after local fit"] = dist(h486["CA"].coord, xyz(q455["CA"], True))
key["CB displacement"] = dist(h486["CB"].coord, xyz(q455["CB"], True))
key["CD(Arg) vs CD(Gln) displacement"] = dist(h486["CD"].coord, xyz(q455["CD"], True))
# per-residue backbone shift in sphere (CA)
sh = {}
for (hc, n), p, q in zip([x for x in near if (x[1] - OFF) in R[mp[x[0]]]], P, Q):
    sh[f"{hc}:{n}"] = dist(p.coord, xyz(q, True))
key["max CA shift in sphere"] = max(sh.items(), key=lambda kv: kv[1])
out["key_distances"] = key
# 6ECK (pS113, FBP, no active-site ligands): same residue, chain A and B
k6 = {}
for ch in "AB":
    r = RK[ch][455]
    k6[ch] = {"resname": r.resname, "polar_contacts_le_3.5A": polar_contacts(r, RK, {"NE2", "OE1"}, 3.5)}
out["6ECK_Gln455"] = k6
json.dump(out, open("data/rat_vs_human.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
