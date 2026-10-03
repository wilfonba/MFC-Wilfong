#!/usr/bin/env python3
"""
Permeability of the square cylinder array from runs of case.py, against Sangani & Acrivos (1982), a low solid-fraction
expansion that turns unphysical past phi ~ 0.4, and Gebart (1992), the near-contact lubrication limit. Each run is compared
with the one that holds at its phi (Sangani & Acrivos below 0.35, Gebart above).

The superficial velocity U is the flow rate through the domain's x = 0 face, which lies between cylinders, over its
height; at steady state K = mu U/(rho g). The references hold for the ordered array (no --jitter). Reads each run's
restart files and simulation.inp.

usage: analyze.py RUN_DIR... [--plot FILE]
"""

import argparse
import glob
import math
import os
import re

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("dirs", nargs="+", help="case directories holding simulation.inp and restart_data/")
parser.add_argument("--plot", help="write K/a^2 against the solid fraction to this file")
args = parser.parse_args()


def sangani_acrivos(phi):
    return (-np.log(phi) - 1.476 + 2 * phi - 1.774 * phi**2 + 4.076 * phi**3) / (8 * phi)


def gebart(phi):
    return 16 / (9 * math.pi * math.sqrt(2)) * (np.sqrt(math.pi / (4 * phi)) - 1) ** 2.5


def value(inp, key):
    return float(re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", inp, re.M).group(1))


def run(d):
    inp = open(os.path.join(d, "simulation.inp")).read()
    nx, ny = int(value(inp, "m")) + 1, int(value(inp, "n")) + 1
    a, g, mu, t_save = value(inp, "patch_ib(1)%radius"), value(inp, "g_x"), 1 / value(inp, "fluid_pp(1)%Re(1)"), value(inp, "t_save")
    xcb = np.fromfile(os.path.join(d, "restart_data", "lustre_x_cb.dat"))
    ycb = np.fromfile(os.path.join(d, "restart_data", "lustre_y_cb.dat"))
    area = (xcb[-1] - xcb[0]) * (ycb[-1] - ycb[0])
    steps = sorted(int(f.split("_")[-1][:-4]) for f in glob.glob(os.path.join(d, "restart_data", "lustre_[0-9]*.dat")))
    t, U, rho = [], [], []
    for s in steps:
        q = np.fromfile(os.path.join(d, "restart_data", f"lustre_{s}.dat")).reshape(-1, ny, nx)  # alpha_rho, mom, E, alpha
        u = q[1] / q[0]
        t.append(s * t_save)
        U.append(0.5 * (u[:, 0].mean() + u[:, -1].mean()))  # flux through the x = 0 seam, both sides of it
        rho.append(q[0][:, 0].mean())
    t, U = np.array(t), np.array(U)
    K = mu * U[-1] / (rho[-1] * g)
    tail = U[int(0.9 * len(U)) :]
    return {"phi": value(inp, "num_ibs") * math.pi * a**2 / area, "K": K / a**2, "Re": rho[-1] * U[-1] * 2 * a / mu, "drift": (tail.max() - tail.min()) / abs(U[-1]), "t": t, "U": U}


runs = {os.path.basename(os.path.normpath(d)): run(d) for d in args.dirs}
print(f"{'run':<28} {'phi':>6} {'Re':>9} {'K/a^2':>10} {'S&A':>10} {'Gebart':>10} {'vs ref':>8} {'U drift':>8}")
for name, r in runs.items():
    sa, gb = sangani_acrivos(r["phi"]), gebart(r["phi"])
    ref = sa if r["phi"] < 0.35 else gb
    print(f"{name:<28} {r['phi']:6.3f} {r['Re']:9.3e} {r['K']:10.4e} {sa:10.4e} {gb:10.4e} {r['K'] / ref - 1:+8.1%} {r['drift']:8.1e}")
print("U drift: spread of U over the last 10% of the run, relative to its final value (steady state needs it small)")

if args.plot:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    phi = np.linspace(0.01, 0.78, 300)
    ax[0].semilogy(phi[phi < 0.4], sangani_acrivos(phi[phi < 0.4]), "k-", lw=1, label="Sangani & Acrivos (1982)")
    ax[0].semilogy(phi[phi > 0.3], gebart(phi[phi > 0.3]), "k--", lw=1, label="Gebart (1992)")
    for name, r in runs.items():
        ax[0].semilogy(r["phi"], r["K"], "o", mfc="none", label=f"{name} (Re = {r['Re']:.2g})")
        ax[1].plot(r["t"], r["U"] / r["U"][-1], label=name)
    ax[0].set_xlabel(r"solid fraction $\phi$")
    ax[0].set_ylabel(r"permeability $K/a^2$")
    ax[0].set_title("(a) permeability of a square cylinder array", loc="left")
    ax[1].set_xlabel(r"$t\ \nu/L^2$")
    ax[1].set_ylabel(r"$U/U_{final}$")
    ax[1].set_title("(b) approach to steady state", loc="left")
    for x in ax:
        x.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)
