"""Geometry of the Arg486 guanidinium relative to the Tyr487 ring (2VGB): CZ-centroid distance and angle between the
ring normal and the centroid->CZ vector (0 deg = cation on the ring axis, typical of cation-pi; > ~60 deg = in-plane).
Run from the analysis/structure folder: python scripts/12_cation_pi_geometry.py  (requires raw/8TBS.cif; see README.md)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from pkcommon import *
for pid in ("2VGB", "8TBS"):
    m = load(pid)
    for ch in m:
        if 486 not in ch or ch[486].resname != "ARG": continue
        r, y = ch[486], ch[487]
        ring = np.array([y[n].coord for n in AROMATIC_RINGS["TYR"]]); cen = ring.mean(0)
        u, s, vt = np.linalg.svd(ring - cen); nrm = vt[2]
        v = r["CZ"].coord - cen; d = np.linalg.norm(v)
        ang = np.degrees(np.arccos(abs(np.dot(v, nrm)) / d))
        print(f"{pid} {ch.id}: CZ-centroid {d:.2f} A, angle to ring normal {ang:.0f} deg")
