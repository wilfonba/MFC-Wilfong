#!/usr/bin/env python3
"""
Temperature snapshots and the Nusselt number of a run of case.py, from its restart files. Nu at each plate is its conductive
heat flux over the pure-conduction value k dT/H, from the wall temperature and the first cell's; the two agree once the rolls
are statistically steady.

usage: analyze.py RUN_DIR [--out FILE] [--panels N]
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
parser.add_argument("dir", help="case directory holding simulation.inp and restart_data/")
parser.add_argument("--out", default="result.png", help="output image (default: %(default)s)")
parser.add_argument("--panels", type=int, default=4, help="temperature snapshots to show (default: %(default)s)")
args = parser.parse_args()


def value(txt, key):
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", txt, re.M).group(1))


sim = open(os.path.join(args.dir, "simulation.inp")).read()
nx, ny, t_save = int(value(sim, "m")) + 1, int(value(sim, "n")) + 1, value(sim, "t_save")
Th, Tc, G, cv = value(sim, "bc_y%Twall_in"), value(sim, "bc_y%Twall_out"), value(sim, "fluid_pp(1)%gamma"), value(sim, "fluid_pp(1)%cv")
x = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_x_cb.dat"))
y = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_y_cb.dat"))
H, dy = y[-1] - y[0], np.diff(y)
steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")))


def temperature(s):
    q = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{s}.dat")).reshape(-1, ny, nx)
    rho, ke = q[0], 0.5 * (q[1] ** 2 + q[2] ** 2) / q[0]
    return (q[3] - ke) / G / (rho * cv / G)  # p = (E - ke)/Gamma, T = p/(rho R), R = cv/Gamma


t, nu_b, nu_t = [], [], []
for s in steps:
    T = temperature(s)
    t.append(s * t_save)
    nu_b.append(np.mean(Th - T[0]) / (dy[0] / 2) * H / (Th - Tc))
    nu_t.append(np.mean(T[-1] - Tc) / (dy[-1] / 2) * H / (Th - Tc))
t, nu_b, nu_t = map(np.array, (t, nu_b, nu_t))

idx = np.linspace(0, len(steps) - 1, min(args.panels, len(steps))).astype(int)
fig = plt.figure(figsize=(3.2 * len(idx), 5.2))
for c, i in enumerate(idx):
    ax = fig.add_subplot(2, len(idx), c + 1)
    im = ax.imshow(temperature(steps[i]), origin="lower", extent=(x[0] / H, x[-1] / H, 0, 1), cmap="RdBu_r", vmin=Tc, vmax=Th)
    ax.set_title(f"t = {t[i]:.1f} s", fontsize=9)
    ax.set_xlabel("x/H")
    if c == 0:
        ax.set_ylabel("y/H")
fig.colorbar(im, ax=fig.axes, shrink=0.5, label="T [K]")
ax = fig.add_subplot(2, 1, 2)
ax.plot(t, nu_b, label="hot plate")
ax.plot(t, nu_t, "--", label="cold plate")
ax.set_xlabel("t [s]")
ax.set_ylabel("Nu")
ax.legend()
fig.savefig(args.out, dpi=150, bbox_inches="tight")

late = t >= t[-1] / 2
print(f"Nu over the second half: hot plate {nu_b[late].mean():.3f}, cold plate {nu_t[late].mean():.3f}; wrote {args.out}")
