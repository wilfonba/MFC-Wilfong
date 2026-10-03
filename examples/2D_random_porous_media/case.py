#!/usr/bin/env python3
"""
Pressure-driven flow through a random 2D porous medium, with the all-Mach pressure projection.

The unit square holds the grains of porous_media.stl (from generate_porous_media.py --flat, see the README) as one
no-slip stationary immersed boundary. The domain is periodic and a uniform body force g stands in for the mean pressure
gradient -dp/dx = rho g; at steady state the superficial velocity U (the mean x velocity over the whole domain, solid
included) gives the permeability K = mu U/(rho g). Units are nondimensional with rho = mu = 1 and a unit domain; the
grain scale d = 4(1 - eps)/S, from the porosity eps and the solid perimeter per unit area S (a disc's diameter), sets
the Reynolds number and the time unit d^2/nu. The sound speed keeps the Mach number near 1e-3, so the projection steps
at the viscous limit; --explicit runs HLLC at the acoustic one as a control.

--Re sets the force from a permeability estimate, a square cylinder array of the same solid fraction and grain scale
(Sangani & Acrivos 1982 / Gebart 1992), to target Re = rho U d/mu; the actual U, and with it K, comes from the run.

The STL was built for a 2048^2 grid (grains >= 30 cells across, throats >= 5): a coarser --N under-resolves the throats.
"""

import argparse
import json
import math
import os
import struct
import sys

parser = argparse.ArgumentParser(description="2D random porous medium (STL immersed boundary), all-Mach pressure projection")
parser.add_argument("--stl", default="porous_media.stl", help="flat 2D STL of the grains in [0, 1]^2, relative to this file (default: %(default)s)")
parser.add_argument("--N", type=int, default=2048, help="cells along each side (default: %(default)s)")
parser.add_argument("--Re", type=float, default=0.1, help="target Reynolds number rho U d/mu (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="cfl_target, viscous or advective; the explicit viscous term is unstable past ~0.3 in 2D (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=2.0, help="end time in grain viscous times d^2/nu (default: %(default)s)")
parser.add_argument("--saves", type=int, default=100, help="number of outputs (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
args, _ = parser.parse_known_args()

stl = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.stl)


def read_flat_stl(path):
    """Solid area and perimeter of a flat binary STL: the boundary is the edges used by one triangle only"""
    with open(path, "rb") as f:
        f.read(80)
        (ntri,) = struct.unpack("<I", f.read(4))
        area, uses = 0.0, {}
        for _ in range(ntri):
            rec = struct.unpack("<12fH", f.read(50))
            v = [rec[3:6], rec[6:9], rec[9:12]]
            if any(abs(p[2] - v[0][2]) > 0 for p in v):
                sys.exit(f"{path} is not flat: regenerate it with generate_porous_media.py --flat")
            area += 0.5 * abs((v[1][0] - v[0][0]) * (v[2][1] - v[0][1]) - (v[1][1] - v[0][1]) * (v[2][0] - v[0][0]))
            for a, b in ((v[0], v[1]), (v[1], v[2]), (v[2], v[0])):
                key = (min(a[:2], b[:2]), max(a[:2], b[:2]))
                uses[key] = uses.get(key, 0) + 1
    perimeter = sum(math.dist(*e) for e, k in uses.items() if k == 1)
    return area, perimeter


def k_stokes(phi):
    """Stokes permeability K/a^2 of a square array: Sangani & Acrivos's expansion, Gebart's lubrication form near contact"""
    if phi < 0.4:
        return (-math.log(phi) - 1.476 + 2 * phi - 1.774 * phi**2 + 4.076 * phi**3) / (8 * phi)
    return 16 / (9 * math.pi * math.sqrt(2)) * (math.sqrt(math.pi / (4 * phi)) - 1) ** 2.5


L, rho, mu = 1.0, 1.0, 1.0
phi, perimeter = read_flat_stl(stl)
d = 4 * phi / perimeter
U = args.Re * mu / (rho * d)
g = mu * U / (rho * k_stokes(phi) * (d / 2) ** 2)
c = 1.0e3 * max(U, mu / (rho * d))  # Mach ~1e-3 against the flow and the viscous velocity scale
gamma, p0 = 4.4, 1.0
t_d = rho * d**2 / mu
print(f"porosity = {1 - phi:.4f}  d = {d:.4e} ({d * args.N:.1f} cells)  target U = {U:.4e}  g = {g:.4e}  d^2/nu = {t_d:.4e}", file=sys.stderr)
if args.N < 2048:
    print(f"warning: the STL was built for N = 2048; N = {args.N} gives its narrowest throats {5 * args.N / 2048:.1f} cells", file=sys.stderr)

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "m": args.N - 1,
            "n": args.N - 1,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": args.cfl,
            "n_start": 0,
            "t_stop": args.tstop * t_d,
            "t_save": args.tstop * t_d / args.saves,
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            "proj_method": "F" if args.explicit else "T",
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu,
            # The mean pressure gradient, as a uniform body force along x
            "bf_x": "T",
            "g_x": g,
            "k_x": 0.0,
            "w_x": 0.0,
            "p_x": 0.0,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * L,
            "patch_icpp(1)%y_centroid": 0.5 * L,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho,
            "patch_icpp(1)%alpha(1)": 1.0,
            # The grains: one no-slip stationary immersed boundary, the STL in its own (domain) coordinates
            "ib": "T",
            "num_ibs": 1,
            "fd_order": 2,
            "num_stl_models": 1,
            "patch_ib(1)%geometry": 5,
            "patch_ib(1)%model_id": 1,
            "patch_ib(1)%x_centroid": 0.0,
            "patch_ib(1)%y_centroid": 0.0,
            "patch_ib(1)%slip": "F",
            "stl_models(1)%model_filepath": stl,
            "stl_models(1)%model_threshold": 0.5,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": gamma * (rho * c**2 / gamma - p0) / (gamma - 1.0),
        }
    )
)
