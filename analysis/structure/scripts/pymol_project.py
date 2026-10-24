"""Helper: project atom coordinates to pixel positions of a PyMOL image rendered with cmd.png(width, height)
in orthoscopic mode, so that labels can be drawn afterwards with matplotlib (scripts/compose_figures.py)."""
import math
import numpy as np
from pymol import cmd

def project(sel, W, H):
    xyz = np.array(cmd.get_coords(sel)).mean(axis=0)
    v = cmd.get_view()
    R = np.array(v[0:9]).reshape(3, 3)
    cam = np.dot(xyz - np.array(v[12:15]), R) + np.array(v[9:12])
    half_h = abs(v[11]) * math.tan(math.radians(cmd.get_setting_float("field_of_view") / 2.0))
    half_w = half_h * W / H
    return float((cam[0] / half_w + 1) / 2 * W), float((1 - (cam[1] / half_h + 1) / 2) * H)
