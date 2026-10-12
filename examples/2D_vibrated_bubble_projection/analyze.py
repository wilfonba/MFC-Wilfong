#!/usr/bin/env python3
"""
Depth below the free surface and area of each bubble over time, against the critical depth h_c = 2 gamma p_amb / (rho g
Gamma^2), from the restart files of a run of case.py (code units: depth in H, time in forcing periods). Gas connected to
the headspace is excluded, so a bubble that reaches the surface leaves the record.

usage: analyze.py RUN_DIR [--out FILE]
"""

import argparse
import glob
import json
import math
import os
import re

import matplotlib
import numpy as np
from scipy import ndimage

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", help="case directory holding simulation.inp and restart_data/")
parser.add_argument("--out", default="result.png", help="output image (default: %(default)s)")
args = parser.parse_args()

sim = open(os.path.join(args.dir, "simulation.inp")).read()


def value(key):
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", sim, re.M).group(1))


nx, ny, t_save = int(value("m")) + 1, int(value("n")) + 1, value("t_save")
w = value("w_y")
# The free surface's height and h_c, which case.py records next to itself
geo = json.load(open(os.path.join(args.dir, "vibrated_bubble.json")))
H, h_c = geo["H"], geo["h_c"]
xb = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_x_cb.dat"))
yb = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_y_cb.dat"))
x, y = 0.5 * (xb[1:] + xb[:-1]), 0.5 * (yb[1:] + yb[:-1])
dA = np.outer(np.diff(yb), np.diff(xb))
steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")))

t, depth, area = [], [], []
for s in steps:
    a2 = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{s}.dat")).reshape(-1, ny, nx)[6]  # alpha of the gas
    # Drop the headspace: gas connected to the top row, with its smeared interface
    lab, _ = ndimage.label(a2 > 0.5)
    head = np.isin(lab, np.unique(lab[-1][lab[-1] > 0]))
    a2 = np.where(ndimage.binary_dilation(head, iterations=4), 0.0, a2)
    row_d, row_a = [], []
    for half in (x < 0, x > 0):
        a = a2[:, half] * dA[:, half]
        row_a.append(a.sum())
        row_d.append(H - (a.sum(1) @ y) / max(a.sum(), 1e-300))
    t.append(s * t_save * w / (2 * math.pi))
    depth.append(row_d)
    area.append(row_a)
t, depth, area = map(np.array, (t, depth, area))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
for i, name in enumerate(["shallow", "deep"]):
    ax1.plot(t, depth[:, i], label=f"{name} bubble")
    ax2.plot(t, area[:, i] / area[0, i], label=f"{name} bubble")
ax1.axhline(h_c, color="k", ls="--", label=r"$h_c = 2\gamma p_{amb}/(\rho g \Gamma^2)$")
ax1.invert_yaxis()
ax1.set_ylabel("depth below the free surface / H")
ax1.legend()
ax2.set_xlabel("t / forcing period")
ax2.set_ylabel("area / initial (sampled)")
fig.tight_layout()
fig.savefig(args.out, dpi=150)
for i, name in enumerate(["shallow", "deep"]):
    print(f"{name}: depth {depth[0, i]:.4f} -> {depth[-1, i]:.4f} H ({'sank' if depth[-1, i] > depth[0, i] else 'rose'}); h_c = {h_c:.4f} H")
print(f"wrote {args.out}")
