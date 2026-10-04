#!/usr/bin/env python3
"""
Snapshots of the bubble column's gas volume fraction, and the gas volume of the bubbles in the water against time, from
the restart files of a run of case.py. The bubbles' volume grows as they rise into lower hydrostatic pressure; with
adiabatic air it follows (p_bottom/p)^(1/1.4) of the depth each bubble has reached.

usage: snapshots.py RUN_DIR [--out FILE] [--panels N]
"""

import argparse
import glob
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", help="case directory holding pre_process.inp, simulation.inp and restart_data/")
parser.add_argument("--out", default="result.png", help="output image (default: %(default)s)")
parser.add_argument("--panels", type=int, default=6, help="snapshots to show (default: %(default)s)")
args = parser.parse_args()


def value(inp, key):
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", inp, re.M).group(1))


inp = open(os.path.join(args.dir, "simulation.inp")).read()
nx, ny, t_save = int(value(inp, "m")) + 1, int(value(inp, "n")) + 1, value(inp, "t_save")
x = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_x_cb.dat"))
y = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_y_cb.dat"))
H = value(open(os.path.join(args.dir, "pre_process.inp")).read(), "patch_icpp(1)%length_y")  # the water depth, below the headspace
steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")))

t, vol = [], []
frames = []
for s in steps:
    q = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{s}.dat")).reshape(-1, ny, nx)
    a_g = q[6]  # alpha_rho(1:2), mom(1:2), E, alpha(1:2), color function
    t.append(s * t_save)
    vol.append(a_g[0.5 * (y[1:] + y[:-1]) < H - 0.01 * H].sum())  # bubbles still in the water
    frames.append(a_g)
t, vol = np.array(t), np.array(vol)

idx = np.linspace(0, len(frames) - 1, min(args.panels, len(frames))).astype(int)
fig = plt.figure(figsize=(2.0 * len(idx) + 4, 6.5))
gs = fig.add_gridspec(1, len(idx) + 1, width_ratios=[1] * len(idx) + [2.4])
for n, i in enumerate(idx):
    ax = fig.add_subplot(gs[n])
    ax.imshow(frames[i], origin="lower", extent=(x[0], x[-1], y[0], y[-1]), cmap="Blues", vmin=0, vmax=1, interpolation="nearest")
    ax.axhline(H, color="0.5", lw=0.5)
    ax.set_title(f"t = {t[i]:.2f} s", fontsize=9)
    ax.set_xticks([])
    if n:
        ax.set_yticks([])
    else:
        ax.set_ylabel("y [m]")
ax = fig.add_subplot(gs[-1])
ax.plot(t, vol / vol[0], "k-")
ax.set_xlabel("t [s]")
ax.set_ylabel("gas volume in the water / initial")
ax.set_title("bubble expansion", fontsize=9)
fig.tight_layout()
fig.savefig(args.out, dpi=150)
print(f"{len(steps)} saves; gas volume in the water {vol[0]:.4g} -> {vol[-1]:.4g} cells ({vol[-1] / vol[0] - 1:+.2%}); wrote {args.out}")
