#!/usr/bin/env python3
"""
Temperature snapshots over the last flicker cycle, a temperature probe on the axis, and its spectrum, whose peak is the flicker
frequency (compared with 0.5 sqrt(g/d)), from the restart files of a run of case.py.

usage: analyze.py RUN_DIR [--out FILE] [--probe Y] [--d D] [--panels N]
"""

import argparse
import glob
import os
import re

import cantera as ct
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dir", help="case directory holding simulation.inp and restart_data/")
parser.add_argument("--out", default="result.png", help="output image (default: %(default)s)")
parser.add_argument("--probe", type=float, default=0.05, help="probe height on the axis [m] (default: %(default)s)")
parser.add_argument("--d", type=float, default=0.01, help="slot width the run used [m] (default: %(default)s)")
parser.add_argument("--panels", type=int, default=5, help="snapshots over the last flicker cycle (default: %(default)s)")
parser.add_argument("--settle", type=float, default=0.1, help="start-up excluded from the spectrum [s] (default: %(default)s)")
args = parser.parse_args()
gas = ct.Solution("h2o2.yaml")

sim = open(os.path.join(args.dir, "simulation.inp")).read()


def value(key):
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", sim, re.M).group(1))


nx, ny, t_save, g = int(value("m")) + 1, int(value("n")) + 1, value("t_save"), abs(value("g_y"))
xb = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_x_cb.dat"))
yb = np.fromfile(os.path.join(args.dir, "restart_data", "lustre_y_cb.dat"))
x, y = 0.5 * (xb[1:] + xb[:-1]), 0.5 * (yb[1:] + yb[:-1])
steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(args.dir, "restart_data", "lustre_[0-9]*.dat")))


def temperature(s):
    """rho, momentum, E, alpha, then the partial densities"""
    q = np.fromfile(os.path.join(args.dir, "restart_data", f"lustre_{s}.dat")).reshape(-1, ny, nx)
    rY = q[5:]
    rho = rY.sum(0)
    e = (q[3] - 0.5 * (q[1] ** 2 + q[2] ** 2) / rho) / rho
    states = ct.SolutionArray(gas, shape=rho.size)
    states.UVY = e.ravel(), 1 / rho.ravel(), (rY / rho).reshape(len(rY), -1).T
    return states.T.reshape(ny, nx)


i_x, i_y = np.argmin(abs(x - 0.5 * (xb[0] + xb[-1]))), np.argmin(abs(y - args.probe))
t = np.array(steps) * t_save
probe = np.array([temperature(s)[i_y, i_x] for s in steps])

late = t >= args.settle
f = np.fft.rfftfreq(late.sum(), t_save)
amp = np.abs(np.fft.rfft(probe[late] - probe[late].mean()))
f_peak = f[1:][np.argmax(amp[1:])]
f_ref = 0.5 * np.sqrt(g / args.d)

period = 1.0 / f_peak
show = [s for s, ts in zip(steps, t) if ts >= t[-1] - period]
show = [show[i] for i in np.linspace(0, len(show) - 1, min(args.panels, len(show))).astype(int)]
fig = plt.figure(figsize=(2.6 * len(show), 8.5))
for c, s in enumerate(show):
    ax = fig.add_subplot(2, len(show), c + 1)
    im = ax.imshow(temperature(s), origin="lower", extent=(xb[0] * 100, xb[-1] * 100, yb[0] * 100, yb[-1] * 100), cmap="inferno", vmin=300, vmax=2100)
    ax.plot(x[i_x] * 100, y[i_y] * 100, "c+")
    ax.set_title(f"t = {s * t_save * 1e3:.0f} ms", fontsize=9)
    ax.set_xlabel("x [cm]")
    if c == 0:
        ax.set_ylabel("y [cm]")
fig.colorbar(im, ax=fig.axes[: len(show)], shrink=0.6, label="T [K]")
ax = fig.add_subplot(4, 1, 3)
ax.plot(t * 1e3, probe)
ax.set_xlabel("t [ms]")
ax.set_ylabel(f"T at y = {args.probe * 100:.0f} cm [K]")
ax = fig.add_subplot(4, 1, 4)
ax.plot(f, amp)
ax.axvline(f_ref, color="k", ls="--", label=r"$0.5\sqrt{g/d}$")
ax.set_xlim(0, 60)
ax.set_xlabel("f [Hz]")
ax.set_ylabel("|FFT|")
ax.legend()
fig.tight_layout()
fig.savefig(args.out, dpi=150)
print(f"flicker frequency {f_peak:.1f} Hz (0.5 sqrt(g/d) = {f_ref:.1f} Hz, resolution {f[1]:.1f} Hz); wrote {args.out}")
