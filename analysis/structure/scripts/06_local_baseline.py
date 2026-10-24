"""Baseline for the local WT-vs-R486W comparison: the same 10 A local-sphere superposition repeated for every residue
of the N, A and C domains of each chain, to see whether the differences around residue 486 exceed the differences
seen elsewhere between the two independently refined models (2VGB, 2VGI). Output: data/local_baseline.tsv + summary."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from Bio.PDB import NeighborSearch, Superimposer
from pkcommon import *
wt = load("2VGB"); mu = load("2VGI")
atoms_wt = [a for a in wt.get_atoms() if a.get_parent().id[0] == " "]
ns = NeighborSearch(atoms_wt)
rows = []
for c in "ABCD":
    for r in protein_residues(wt[c]):
        n = r.id[1]
        if domain_of(n) == "B" or n not in mu[c]: continue
        near = set()
        for a in r:
            for b in ns.search(a.coord, 10.0):
                near.add((b.get_parent().get_parent().id, b.get_parent().id[1]))
        P, Q, Ph, Qh = [], [], [], []
        for (cc, nn) in near:
            if nn not in mu[cc]: continue
            r1, r2 = wt[cc][nn], mu[cc][nn]
            for a in r1:
                nm = a.get_id()
                if nm not in r2: continue
                if r1.resname != r2.resname and nm not in ("N", "CA", "C", "O", "CB"): continue
                if (cc, nn) == (c, n) and nm not in ("N", "CA", "C", "O", "CB"): continue
                Ph.append(a); Qh.append(r2[nm])
                if nm == "CA": P.append(a); Q.append(r2[nm])
        if len(P) < 10: continue
        s = Superimposer(); s.set_atoms(P, Q); R, t = s.rotran
        d = np.linalg.norm(np.array([a.coord for a in Ph]) - (np.array([b.coord for b in Qh]) @ R + t), axis=1)
        rows.append((c, n, r.resname, len(P), s.rms, float(np.sqrt(np.mean(d ** 2))), float(d.max()), int((d > 0.5).sum()) / len(d)))
import csv
with open("data/local_baseline.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["chain", "resnum", "resname_WT", "n_CA", "CA_rmsd_local", "heavy_rmsd_local", "max_shift", "frac_atoms_gt_0.5A"])
    for x in rows: w.writerow([x[0], x[1], x[2], x[3], round(x[4], 3), round(x[5], 3), round(x[6], 2), round(x[7], 3)])
summ = {}
for c in "ABCD":
    R = [x for x in rows if x[0] == c]
    ca = np.array([x[4] for x in R]); hv = np.array([x[5] for x in R]); mx = np.array([x[6] for x in R])
    me = [x for x in R if x[1] == 486][0]
    summ[c] = {"n_centres": len(R), "CA_rmsd_median": round(float(np.median(ca)), 3), "CA_rmsd_p95": round(float(np.percentile(ca, 95)), 3),
               "heavy_rmsd_median": round(float(np.median(hv)), 3), "heavy_rmsd_p95": round(float(np.percentile(hv, 95)), 3),
               "max_shift_median": round(float(np.median(mx)), 2),
               "res486_CA_rmsd": round(me[4], 3), "res486_CA_percentile": round(float((ca < me[4]).mean() * 100), 1),
               "res486_heavy_rmsd": round(me[5], 3), "res486_heavy_percentile": round(float((hv < me[5]).mean() * 100), 1),
               "res486_max_shift": round(me[6], 2), "res486_max_shift_percentile": round(float((mx < me[6]).mean() * 100), 1)}
json.dump(summ, open("data/local_baseline_summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
