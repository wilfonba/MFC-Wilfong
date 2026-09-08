#!/usr/bin/env pvpython
"""
Compute the Hysing et al. (2009) rising-bubble benchmark quantities from MFC's
Silo output with ParaView's pvpython:

  * center of mass    Y_c(t) = int(alpha2*y) / int(alpha2)
  * rise velocity     V_c(t) = int(alpha2*v)  / int(alpha2)
  * circularity       xi(t)  = 2*sqrt(pi*A)/P, A = int(alpha2), P = length of
                               the alpha2 = 0.5 contour
  * bubble shape at the final saved time (the paper reports it at t = 3)

The sharp-interface bubble region of the paper is replaced by alpha2 (bubble
volume fraction) weighting, the standard diffuse-interface analog. Reference
values for test case 1 (Hysing et al. 2009, Table 2ff): xi_min ~ 0.90 near
t ~ 1.9, V_c,max ~ 0.242 at t ~ 0.92, Y_c(3) ~ 1.081.

Usage:
  pvpython benchmark_quantities.py [case_dir] [--t-save 0.1] [--case 1|2]

--case overlays the finest-level results of the three benchmark groups
(TP2D, FreeLIFE, MooNMD) from reference/ on the plot; see reference/README.md.

Writes case_dir/benchmark_quantities.csv (t, y_c, v_c, circularity, area,
perimeter), case_dir/bubble_shape.csv (x0, y0, x1, y1 segments of the final-time contour), and
case_dir/benchmark_quantities.png (2x2 plot of the three time series and the
final bubble shape).
"""

import argparse
import csv
import glob
import math
import os
import sys

import matplotlib
import numpy as np
from paraview import servermanager

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from paraview.simple import (
    Calculator,
    CellDatatoPointData,
    Contour,
    Delete,
    IntegrateVariables,
    MergeBlocks,
    OpenDataFile,
)

parser = argparse.ArgumentParser(description="Hysing et al. (2009) benchmark quantities from MFC Silo output")
parser.add_argument("case_dir", nargs="?", default=os.path.dirname(os.path.abspath(__file__)), help="MFC case directory (default: this script's directory)")
parser.add_argument("--t-save", type=float, default=0.1, help="Save interval used by the case (default: 0.1)")
parser.add_argument("--case", type=int, choices=[1, 2], help="Hysing et al. (2009) test case whose reference data to overlay (default: none)")
args = parser.parse_args()

root_dir = os.path.join(args.case_dir, "root")
files = sorted(glob.glob(os.path.join(root_dir, "collection_*.silo")), key=lambda f: int(os.path.basename(f)[len("collection_") : -len(".silo")]))
if not files:
    sys.exit(f"No Silo files found under {root_dir} -- run the case (and post_process) first.")


def integrate(source):
    """Return {array name: integrated value} for all point/cell data of source."""
    integ = IntegrateVariables(Input=source)
    data = servermanager.Fetch(integ)
    result = {}
    for attr in (data.GetPointData(), data.GetCellData()):
        for i in range(attr.GetNumberOfArrays()):
            arr = attr.GetArray(i)
            result[arr.GetName()] = arr.GetTuple1(0)
    Delete(integ)
    return result


rows = []
shape_segments = []
for index, fname in enumerate(files):
    t = index * args.t_save

    reader = OpenDataFile(fname)
    reader.CellArrayStatus = ["alpha2", "vel2"]
    merged = MergeBlocks(Input=reader)
    pointed = CellDatatoPointData(Input=merged)

    weighted = Calculator(Input=pointed)
    weighted.ResultArrayName = "a2y"
    weighted.Function = "alpha2*coordsY"
    weighted2 = Calculator(Input=weighted)
    weighted2.ResultArrayName = "a2v"
    weighted2.Function = "alpha2*vel2"

    bulk = integrate(weighted2)
    area = bulk["alpha2"]
    y_c = bulk["a2y"] / area
    v_c = bulk["a2v"] / area

    interface = Contour(Input=pointed)
    interface.ContourBy = ["POINTS", "alpha2"]
    interface.Isosurfaces = [0.5]
    perimeter = integrate(interface).get("Length", float("nan"))
    circularity = 2.0 * math.sqrt(math.pi * area) / perimeter if perimeter > 0 else float("nan")

    if fname == files[-1]:
        contour_data = servermanager.Fetch(interface)
        for i in range(contour_data.GetNumberOfCells()):
            ids = contour_data.GetCell(i).GetPointIds()
            shape_segments.append(contour_data.GetPoint(ids.GetId(0))[:2] + contour_data.GetPoint(ids.GetId(1))[:2])

    rows.append((t, y_c, v_c, circularity, area, perimeter))
    print(f"t = {t:5.2f}  y_c = {y_c:.6f}  v_c = {v_c:+.6f}  circularity = {circularity:.6f}")

    for src in (interface, weighted2, weighted, pointed, merged, reader):
        Delete(src)

out_csv = os.path.join(args.case_dir, "benchmark_quantities.csv")
with open(out_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["t", "y_c", "v_c", "circularity", "area", "perimeter"])
    writer.writerows(rows)
print(f"Wrote {out_csv}")

# Bubble shape at the final saved time as line segments (same layout as reference/)
shape_csv = os.path.join(args.case_dir, "bubble_shape.csv")
with open(shape_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["x0", "y0", "x1", "y1"])
    writer.writerows(shape_segments)
print(f"Wrote {shape_csv}")

# 2x2 summary plot: y_c(t), v_c(t), circularity(t), final bubble shape,
# optionally against the benchmark groups' reference data
ref_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference")
ref_groups = ("tp2d", "freelife", "moonmd") if args.case else ()
ref = {g: np.loadtxt(os.path.join(ref_dir, f"c{args.case}_{g}.csv"), delimiter=",", skiprows=1) for g in ref_groups}
ref_shape = {g: np.loadtxt(os.path.join(ref_dir, f"c{args.case}_{g}_shape.csv"), delimiter=",", skiprows=1) for g in ref_groups}

data = np.array(rows)
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
# (MFC column, reference column, label)
for ax, col, ref_col, label in zip(axes.flat, (1, 2, 3), (3, 4, 2), ("center of mass $y_c$", "rise velocity $v_c$", "circularity")):
    for g, r in ref.items():
        ax.plot(r[:, 0], r[:, ref_col], "--", lw=1, label=g)
    ax.plot(data[:, 0], data[:, col], "k", label="MFC")
    ax.set_xlabel("$t$")
    ax.set_ylabel(label)
    ax.grid(True)
ax = axes[1, 1]


def plot_segments(ax, seg, *plot_args, **plot_kwargs):
    """Draw (x0, y0, x1, y1) rows as one NaN-separated polyline."""
    nan = np.full(len(seg), np.nan)
    ax.plot(np.column_stack([seg[:, 0], seg[:, 2], nan]).ravel(), np.column_stack([seg[:, 1], seg[:, 3], nan]).ravel(), *plot_args, **plot_kwargs)


for g, r in ref_shape.items():
    plot_segments(ax, r, "--", lw=1, label=g)
if shape_segments:
    plot_segments(ax, np.array(shape_segments), "k", label="MFC")
if ref:
    axes[0, 0].legend()
ax.set_aspect("equal")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_title(f"bubble shape at t = {rows[-1][0]:.2f}")
ax.grid(True)
fig.tight_layout()
out_png = os.path.join(args.case_dir, "benchmark_quantities.png")
fig.savefig(out_png, dpi=150)
print(f"Wrote {out_png}")
