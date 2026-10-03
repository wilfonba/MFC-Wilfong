#!/usr/bin/env python3
"""
Generate a random 2D porous medium in the unit square [0, 1] x [0, 1].

The solid phase is a set of non-overlapping grains placed by random sequential
addition. Grain shapes can be:

  * rounded  - irregular, lumpy, elongated grains (weathered sand / sediment)
  * angular  - random convex polygons with grid-scale rounded corners (crushed rock)
  * circle   - idealized discs

Grain sizes and spacings are expressed in grid cells of the target simulation
resolution, so the geometry is guaranteed to be resolved by the solver:

  * every part of every grain is >= --smallest-media cells thick. Each grain is
    morphologically opened by a disc of that diameter, i.e. it is exactly a
    union of such discs, so no tip, neck or lobe can be thinner.
  * every pore throat (grain-to-grain gap) is >= --min-gap cells

Outputs
  * <output>.stl  - closed, extruded solid (binary STL) for import into a CFD code, or
                    with --flat the z = 0 triangulation alone, the form MFC reads for a
                    2D immersed boundary (it takes the edges used once as the boundary)
  * <output>.png  - preview: exact geometry + the voxelized mask at grid resolution

Requires: numpy, matplotlib, shapely>=2.1

Example
  python generate_porous_media.py --grid-resolution 256 --smallest-media 6 \
      --shape rounded --porosity 0.6 --seed 42 --output media
"""

import argparse
import struct
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shapely
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as PolygonPatch
from shapely import affinity
from shapely.geometry import Point, Polygon


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate random 2D porous media in [0,1]x[0,1] and export an STL.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--grid-resolution", type=int, required=True, help="Simulation grid cells per unit length (domain is N x N cells).")
    p.add_argument("--smallest-media", type=int, required=True, help="Minimum grain feature size (thickness of any part of a grain), " "in grid cells.")
    p.add_argument("--largest-media", type=int, default=None, help="Maximum grain size (equivalent-area diameter), in grid cells " "(default: 4x smallest).")
    p.add_argument("--shape", choices=["rounded", "angular", "circle"], default="rounded", help="Grain shape family.")
    p.add_argument("--roughness", type=float, default=0.3, help="Lumpiness of 'rounded' grains (0 = smooth ellipse, ~0.5 = very lumpy).")
    p.add_argument("--max-aspect", type=float, default=1.8, help="Maximum grain elongation (length/width) for non-circular shapes.")
    p.add_argument("--porosity", type=float, default=0.6, help="Target porosity (void fraction). Packing stops when reached " "or when no more grains fit.")
    p.add_argument("--min-gap", type=int, default=2, help="Minimum gap between grains (pore throat width), in grid cells.")
    p.add_argument("--margin", type=int, default=None, help="Minimum clearance between grains and the domain edges, in grid " "cells (default: equal to --min-gap).")
    p.add_argument("--depth", type=float, default=None, help="Extrusion depth in z for the STL (default: one grid cell, 1/N).")
    p.add_argument("--max-attempts", type=int, default=200000, help="Maximum random placement attempts.")
    p.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    p.add_argument("--output", type=str, default="porous_media", help="Output file prefix (writes <prefix>.stl and <prefix>.png).")
    p.add_argument("--ascii", action="store_true", help="Write ASCII STL instead of binary.")
    p.add_argument("--flat", action="store_true", help="Write only the grains' z = 0 triangulation (a 2D model) instead of " "the extruded solid.")
    args = p.parse_args()

    if args.grid_resolution < 2:
        p.error("--grid-resolution must be >= 2")
    if args.smallest_media < 1:
        p.error("--smallest-media must be >= 1")
    if args.largest_media is None:
        args.largest_media = 4 * args.smallest_media
    if args.largest_media < args.smallest_media:
        p.error("--largest-media must be >= --smallest-media")
    if not 0.0 < args.porosity < 1.0:
        p.error("--porosity must be in (0, 1)")
    if args.min_gap < 0:
        p.error("--min-gap must be >= 0")
    if args.roughness < 0:
        p.error("--roughness must be >= 0")
    if args.max_aspect < 1:
        p.error("--max-aspect must be >= 1")
    if args.margin is None:
        args.margin = args.min_gap
    if args.largest_media + 2 * args.margin > args.grid_resolution:
        p.error("--largest-media plus margins does not fit in the domain at this resolution")
    return args


# Grain shapes
class GrainFactory:
    """Builds random grain polygons centered (by centroid) at the origin."""

    def __init__(self, rng, kind, r_min, roughness, max_aspect, tol):
        self.rng = rng
        self.kind = kind
        self.r_min = r_min
        self.roughness = roughness
        self.max_aspect = max_aspect
        self.tol = tol  # max deviation from the ideal curve (fraction of a cell)

    def _quad_segs(self, r):
        # Segments per quarter circle so the chord error is below tol.
        if r <= self.tol:
            return 2
        return int(np.clip(np.ceil(np.pi / (4 * np.arccos(1 - self.tol / r))), 2, 64))

    def _rounded(self):
        rng = self.rng
        th = np.linspace(0, 2 * np.pi, 180, endpoint=False)
        r = np.ones_like(th)
        for k in range(2, 9):
            amp = self.roughness * rng.normal() / k**1.3
            r += amp * np.cos(k * th + rng.uniform(0, 2 * np.pi))
        r = np.clip(r, 0.35, None)
        return Polygon(np.column_stack([r * np.cos(th), r * np.sin(th)]))

    def _angular(self):
        rng = self.rng
        n = rng.integers(5, 10)
        # Jittered angles avoid degenerate slivers; random radii give facets of unequal length.
        th = (np.arange(n) + rng.uniform(-0.35, 0.35, n)) * 2 * np.pi / n
        r = rng.uniform(0.65, 1.0, n)
        return Polygon(np.column_stack([r * np.cos(th), r * np.sin(th)])).convex_hull

    def make(self, diameter):
        """Grain whose area matches a disc of `diameter` (before opening)."""
        if self.kind == "circle":
            r = 0.5 * diameter
            return Point(0, 0).buffer(r, quad_segs=self._quad_segs(r))

        g = self._rounded() if self.kind == "rounded" else self._angular()
        aspect = self.rng.uniform(1.0, self.max_aspect)
        g = affinity.scale(g, np.sqrt(aspect), 1 / np.sqrt(aspect), origin=(0, 0))
        g = affinity.rotate(g, self.rng.uniform(0, 360), origin=(0, 0))
        target_area = np.pi * diameter**2 / 4
        s = np.sqrt(target_area / g.area)
        g = affinity.scale(g, s, s, origin=(0, 0))

        # Morphological opening: removes any feature thinner than 2*r_min and
        # rounds sharp corners to radius r_min, so the grain is exactly a union
        # of discs of diameter --smallest-media.
        q = self._quad_segs(self.r_min)
        opened = g.buffer(-self.r_min, quad_segs=q).buffer(self.r_min, quad_segs=q)
        if opened.is_empty:
            return Point(0, 0).buffer(0.5 * diameter, quad_segs=self._quad_segs(0.5 * diameter))
        if opened.geom_type == "MultiPolygon":
            opened = max(opened.geoms, key=lambda p: p.area)
        g = Polygon(opened.exterior).simplify(self.tol, preserve_topology=True)
        c = g.centroid
        return affinity.translate(g, -c.x, -c.y)


# Packing
def pack_grains(rng, factory, d_min, d_max, gap, margin, target_solid, max_attempts):
    """Random sequential addition of non-overlapping grains in the unit square.

    A uniform spatial hash on grain centroids (with per-grain bounding radii)
    limits the exact polygon-distance test to nearby grains. Grains are placed
    largest-first, which packs considerably denser than random insertion order.
    """
    bin_size = d_max + gap
    nbins = max(1, int(1.0 / bin_size))
    bins = [[[] for _ in range(nbins)] for _ in range(nbins)]
    placed, centers, brad = [], [], []
    max_brad = [0.0]

    def bin_of(v):
        return min(nbins - 1, max(0, int(v * nbins)))

    def draw_diameter():
        # Uniform in area between d_min and d_max.
        return float(np.sqrt(rng.uniform(d_min**2, d_max**2)))

    def try_place(shape, br):
        minx, miny, maxx, maxy = shape.bounds
        xlo, xhi = margin - minx, 1.0 - margin - maxx
        ylo, yhi = margin - miny, 1.0 - margin - maxy
        if xlo >= xhi or ylo >= yhi:
            return False
        x, y = rng.uniform(xlo, xhi), rng.uniform(ylo, yhi)

        reach = br + max_brad[0] + gap
        nb = int(np.ceil(reach * nbins))
        bx, by = bin_of(x), bin_of(y)
        near = []
        for i in range(max(0, bx - nb), min(nbins, bx + nb + 1)):
            for j in range(max(0, by - nb), min(nbins, by + nb + 1)):
                for k in bins[i][j]:
                    cx, cy = centers[k]
                    lim = br + brad[k] + gap
                    if (x - cx) ** 2 + (y - cy) ** 2 < lim * lim:
                        near.append(k)

        cand = affinity.translate(shape, x, y)
        if near and np.any(shapely.dwithin([placed[k] for k in near], cand, gap)):
            return False
        placed.append(cand)
        centers.append((x, y))
        brad.append(br)
        max_brad[0] = max(max_brad[0], br)
        bins[bx][by].append(len(placed) - 1)
        return True

    def new_grain(d):
        g = factory.make(d)
        br = float(np.max(np.hypot(*np.asarray(g.exterior.coords).T)))
        return g, br

    solid_area = 0.0  # domain area is 1, so this is the solid fraction

    # Pass 1: draw a set of sizes whose total area hits the target, then place
    # them largest-first (big grains are hardest to fit later).
    batch, area = [], 0.0
    while area < target_solid:
        d = draw_diameter()
        batch.append(d)
        area += np.pi * d * d / 4
    batch.sort(reverse=True)

    attempts, per_grain = 0, 500
    for d in batch:
        g, br = new_grain(d)
        for t in range(per_grain):
            attempts += 1
            if t and t % 100 == 0:  # a different shape/orientation may fit where this one won't
                g, br = new_grain(d)
            if try_place(g, br):
                solid_area += g.area
                break
        if solid_area >= target_solid or attempts >= max_attempts:
            break

    # Pass 2: top up with fresh random grains into the remaining space.
    while solid_area < target_solid and attempts < max_attempts:
        attempts += 1
        g, br = new_grain(draw_diameter())
        if try_place(g, br):
            solid_area += g.area

    return placed


def voxelize(grains, n):
    """Solid mask on the N x N simulation grid (cell-center sampling). mask[j, i] -> (x_i, y_j)."""
    xc = (np.arange(n) + 0.5) / n
    mask = np.zeros((n, n), dtype=bool)
    for g in grains:
        minx, miny, maxx, maxy = g.bounds
        i0, i1 = max(0, int(minx * n)), min(n, int(maxx * n) + 1)
        j0, j1 = max(0, int(miny * n)), min(n, int(maxy * n) + 1)
        X, Y = np.meshgrid(xc[i0:i1], xc[j0:j1])
        mask[j0:j1, i0:i1] |= shapely.contains_xy(g, X, Y)
    return mask


# STL export
def cap_triangles(poly, z):
    """CCW (+z facing) triangles of a polygon at height z, from a constrained Delaunay
    triangulation that uses the ring's own vertices, so they share its edges exactly."""
    tris = []
    for t in shapely.constrained_delaunay_triangles(poly).geoms:
        a, b, c = np.asarray(t.exterior.coords)[:3]
        if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) < 0:
            b, c = c, b  # make CCW
        tris.append(((a[0], a[1], z), (b[0], b[1], z), (c[0], c[1], z)))
    return tris


def extrude_polygon(poly, z0, z1):
    """Closed, outward-oriented triangulation of a polygon extruded in z."""
    poly = shapely.geometry.polygon.orient(poly, sign=1.0)  # CCW exterior
    ring = np.asarray(poly.exterior.coords)[:-1]
    nxt = np.roll(ring, -1, axis=0)
    tris = []
    # Side walls: for a CCW ring, (b0, b1, t1) has an outward normal.
    for (x0, y0), (x1, y1) in zip(ring, nxt):
        b0, b1, t0, t1 = (x0, y0, z0), (x1, y1, z0), (x0, y0, z1), (x1, y1, z1)
        tris.append((b0, b1, t1))
        tris.append((b0, t1, t0))
    # Caps share the ring's vertices, so cap and wall edges match exactly -> watertight.
    for a, b, c in cap_triangles(poly, z1):
        tris.append((a, b, c))  # +z
        tris.append(((a[0], a[1], z0), (c[0], c[1], z0), (b[0], b[1], z0)))  # -z
    return tris


def write_stl(path, triangles, ascii_mode=False, name="porous_media"):
    tri = np.asarray(triangles, dtype=np.float64)  # (M, 3, 3)
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)

    if ascii_mode:
        with open(path, "w") as f:
            f.write(f"solid {name}\n")
            for n, t in zip(nrm, tri):
                f.write(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n    outer loop\n")
                for v in t:
                    f.write(f"      vertex {v[0]:.9e} {v[1]:.9e} {v[2]:.9e}\n")
                f.write("    endloop\n  endfacet\n")
            f.write(f"endsolid {name}\n")
    else:
        rec = np.zeros(len(tri), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        rec["n"], rec["v"] = nrm, tri
        with open(path, "wb") as f:
            f.write(name.encode()[:80].ljust(80, b"\0"))
            f.write(struct.pack("<I", len(tri)))
            f.write(rec.tobytes())


# Preview
def save_preview(path, grains, mask, args, porosity_exact, porosity_grid):
    n = args.grid_resolution
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5), constrained_layout=True)

    ax = axes[0]
    patches = [PolygonPatch(np.asarray(g.exterior.coords)) for g in grains]
    ax.add_collection(PatchCollection(patches, facecolor="#3b3b3b", edgecolor="none"))
    ax.set_title(f"Geometry (STL)  -  {len(grains)} grains, porosity {porosity_exact:.3f}")

    ax = axes[1]
    ax.imshow(mask, origin="lower", extent=(0, 1, 0, 1), cmap="Greys", interpolation="nearest", vmin=0, vmax=1.4)
    ax.set_title(f"Voxelized at {n}x{n}  -  porosity {porosity_grid:.3f}")
    if n <= 128:
        ticks = np.linspace(0, 1, n + 1)
        ax.set_xticks(ticks, minor=True)
        ax.set_yticks(ticks, minor=True)
        ax.grid(which="minor", color="#cccccc", linewidth=0.3)

    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    fig.suptitle(f"grid={n}, shape={args.shape}, smallest={args.smallest_media} cells, " f"largest={args.largest_media} cells, min gap={args.min_gap} cells, " f"seed={args.seed}")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    n = args.grid_resolution
    dx = 1.0 / n
    rng = np.random.default_rng(args.seed)

    d_min = args.smallest_media * dx
    d_max = args.largest_media * dx
    factory = GrainFactory(rng, args.shape, r_min=0.5 * d_min, roughness=args.roughness, max_aspect=args.max_aspect, tol=0.02 * dx)
    grains = pack_grains(
        rng,
        factory,
        d_min,
        d_max,
        gap=args.min_gap * dx,
        margin=args.margin * dx,
        target_solid=1.0 - args.porosity,
        max_attempts=args.max_attempts,
    )
    if not grains:
        sys.exit("No grains could be placed; check your parameters.")

    porosity_exact = 1.0 - sum(g.area for g in grains)
    mask = voxelize(grains, n)
    porosity_grid = 1.0 - mask.mean()

    # STL: each grain extruded from z=0 to z=depth, or its z=0 face alone with --flat.
    depth = 0.0 if args.flat else (args.depth if args.depth is not None else dx)
    tris = []
    for g in grains:
        tris.extend(cap_triangles(g, 0.0) if args.flat else extrude_polygon(g, 0.0, depth))

    stl_path, png_path = f"{args.output}.stl", f"{args.output}.png"
    write_stl(stl_path, tris, ascii_mode=args.ascii)
    save_preview(png_path, grains, mask, args, porosity_exact, porosity_grid)

    d_eq = 2 * np.sqrt(np.array([g.area for g in grains]) / np.pi) * n
    print(f"Grains placed       : {len(grains)} ({args.shape})")
    print(f"Equiv. diam (cells) : min {d_eq.min():.2f}, max {d_eq.max():.2f}")
    print(f"Porosity (target)   : {args.porosity:.4f}")
    print(f"Porosity (geometry) : {porosity_exact:.4f}")
    print(f"Porosity ({n}x{n} grid): {porosity_grid:.4f}")
    if porosity_exact > args.porosity + 0.01:
        print("  note: could not pack to target porosity; try smaller --min-gap, " "a wider size range, or more --max-attempts.")
    print(f"STL                 : {stl_path} ({len(tris)} triangles, z in [0, {depth:g}])")
    print(f"Preview             : {png_path}")


if __name__ == "__main__":
    main()
