#!/usr/bin/env python3
"""
Snapshots of the flame-vortex interaction (temperature, with vorticity contours) and the flame's length, the T = 1000 K contour,
over time, from the restart files of a run of case.py. --compare overlays a second run's flame length (e.g. an --explicit one).

usage: analyze.py RUN_DIR [--compare RUN_DIR] [--out FILE] [--panels N]
"""

import argparse
import glob
import os
import re

import cantera as ct
import contourpy
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", help="case directory holding simulation.inp and restart_data/")
parser.add_argument("--compare", help="second run whose flame length to overlay")
parser.add_argument("--out", default="result.png", help="output image (default: %(default)s)")
parser.add_argument("--panels", type=int, default=4, help="snapshots to show (default: %(default)s)")
args = parser.parse_args()
gas = ct.Solution("h2o2.yaml")


def value(txt, key):
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", txt, re.M).group(1))


def load(d):
    sim = open(os.path.join(d, "simulation.inp")).read()
    nx, ny, t_save = int(value(sim, "m")) + 1, int(value(sim, "n")) + 1, value(sim, "t_save")
    x = np.fromfile(os.path.join(d, "restart_data", "lustre_x_cb.dat"))
    y = np.fromfile(os.path.join(d, "restart_data", "lustre_y_cb.dat"))
    steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(d, "restart_data", "lustre_[0-9]*.dat")))
    return nx, ny, t_save, 0.5 * (x[1:] + x[:-1]), 0.5 * (y[1:] + y[:-1]), steps


def fields(d, nx, ny, s):
    """Temperature and vorticity: rho, momentum, E, alpha, then the partial densities"""
    q = np.fromfile(os.path.join(d, "restart_data", f"lustre_{s}.dat")).reshape(-1, ny, nx)
    rY = q[5:]
    rho = rY.sum(0)
    u, v = q[1] / rho, q[2] / rho
    e = (q[3] - 0.5 * rho * (u**2 + v**2)) / rho
    states = ct.SolutionArray(gas, shape=rho.size)
    states.UVY = e.ravel(), 1 / rho.ravel(), (rY / rho).reshape(len(rY), -1).T
    return states.T.reshape(ny, nx), u, v


def flame_length(x, y, T):
    return sum(np.sum(np.hypot(*np.diff(line, axis=0).T)) for line in contourpy.contour_generator(x, y, T).lines(1000.0))


nx, ny, t_save, x, y, steps = load(args.dir)
idx = np.linspace(0, len(steps) - 1, min(args.panels, len(steps))).astype(int)
fig = plt.figure(figsize=(3.6 * len(idx), 6.4))
lengths = []
for i, s in enumerate(steps):
    T, u, v = fields(args.dir, nx, ny, s)
    lengths.append(flame_length(x, y, T))
    if i in idx:
        c = list(idx).index(i)
        ax = fig.add_subplot(2, len(idx), c + 1)
        im = ax.imshow(T, origin="lower", extent=(x[0] * 1e3, x[-1] * 1e3, y[0] * 1e3, y[-1] * 1e3), cmap="inferno", vmin=300, vmax=1900)
        ax.set_xlim(0, 10)
        w = np.gradient(v, x, axis=1) - np.gradient(u, y, axis=0)
        ax.contour(x * 1e3, y * 1e3, w, levels=[-2e3, -5e2, 5e2, 2e3], colors=["c", "c", "w", "w"], linewidths=0.6)
        ax.set_title(f"t = {s * t_save * 1e3:.2f} ms", fontsize=9)
        ax.set_xlabel("x [mm]")
        if c == 0:
            ax.set_ylabel("y [mm]")
fig.colorbar(im, ax=fig.axes, shrink=0.45, label="T [K]")
ax = fig.add_subplot(2, 1, 2)
ax.plot(np.array(steps) * t_save * 1e3, np.array(lengths) * 1e3, "o-", label=os.path.basename(os.path.normpath(args.dir)))
if args.compare:
    nx2, ny2, t2, x2, y2, s2 = load(args.compare)
    l2 = [flame_length(x2, y2, fields(args.compare, nx2, ny2, s)[0]) for s in s2]
    ax.plot(np.array(s2) * t2 * 1e3, np.array(l2) * 1e3, "s--", label=os.path.basename(os.path.normpath(args.compare)))
ax.set_xlabel("t [ms]")
ax.set_ylabel("flame length (T = 1000 K) [mm]")
ax.legend()
fig.savefig(args.out, dpi=150, bbox_inches="tight")
print(f"flame length {lengths[0] * 1e3:.2f} -> {lengths[-1] * 1e3:.2f} mm; wrote {args.out}")
