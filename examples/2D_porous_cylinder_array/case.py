#!/usr/bin/env python3
"""
Pressure-driven flow through a square array of cylinders, the canonical model of a fibrous porous medium, with the
all-Mach pressure projection.

A periodic domain of nx x ny unit cells of side L = 1 holds a no-slip cylinder (a stationary immersed boundary) of radius
a in each, a solid fraction phi = pi a^2; --jitter moves each cylinder randomly within its cell for a disordered medium.
A uniform body force g stands in for the mean pressure gradient -dp/dx = rho g; at steady state the superficial velocity
U (the flow rate through the x = 0 face over its height) gives the permeability K = mu U/(rho g), for the ordered array
compared by analyze.py with Sangani & Acrivos (Int. J. Multiphase Flow 8:193-206, 1982) and Gebart (J. Compos. Mater.
26:1100-1133, 1992). Units are nondimensional with rho = mu = 1; the sound speed keeps the Mach number near 1e-3, so the
projection steps at the viscous limit; --explicit runs HLLC at the acoustic one, here only a few times smaller, as a control.

--Re sets the force from the Stokes estimate of K to target Re = rho U (2a)/mu; inertia lowers the Re actually reached
once it is of order one or more, and analyze.py reports it.
"""

import argparse
import json
import math
import random
import sys

parser = argparse.ArgumentParser(description="2D square cylinder array (porous medium), all-Mach pressure projection")
parser.add_argument("--phi", type=float, default=0.2, help="solid fraction pi a^2, below pi/4 (default: %(default)s)")
parser.add_argument("--Re", type=float, default=0.1, help="target Reynolds number rho U (2a)/mu (default: %(default)s)")
parser.add_argument("--nx", type=int, default=16, help="cylinders along x (default: %(default)s)")
parser.add_argument("--ny", type=int, default=16, help="cylinders along y (default: %(default)s)")
parser.add_argument("--ppc", type=int, default=64, help="cells across each unit cell (default: %(default)s)")
parser.add_argument("--jitter", type=float, default=0.0, help="random offset of each cylinder, as a fraction of its clearance to the cell edge (default: %(default)s)")
parser.add_argument("--seed", type=int, default=1, help="random seed for --jitter (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="cfl_target, viscous or advective; the explicit viscous term is unstable past ~0.3 in 2D (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=0.5, help="end time in viscous times L^2/nu (default: %(default)s)")
parser.add_argument("--saves", type=int, default=50, help="number of outputs (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
args, _ = parser.parse_known_args()
if not 0 < args.phi < math.pi / 4:
    parser.error("--phi must be in (0, pi/4): the cylinders touch at pi/4")
if not 0 <= args.jitter < 1:
    parser.error("--jitter must be in [0, 1) so each cylinder stays inside its cell")

L, rho, mu = 1.0, 1.0, 1.0
a = math.sqrt(args.phi / math.pi) * L


def k_stokes(phi):
    """Stokes permeability K/a^2 of a square array: Sangani & Acrivos's expansion, Gebart's lubrication form near contact"""
    if phi < 0.4:
        return (-math.log(phi) - 1.476 + 2 * phi - 1.774 * phi**2 + 4.076 * phi**3) / (8 * phi)
    return 16 / (9 * math.pi * math.sqrt(2)) * (math.sqrt(math.pi / (4 * phi)) - 1) ** 2.5


U = args.Re * mu / (rho * 2 * a)
g = mu * U / (rho * k_stokes(args.phi) * a**2)
c = 1.0e3 * max(U, mu / (rho * L))  # Mach ~1e-3 against the flow and the viscous velocity scale
gamma, p0 = 4.4, 1.0
print(f"a = {a:.4f}  target U = {U:.4e}  g = {g:.4e}  viscous time L^2/nu = {rho * L**2 / mu:g}  cylinders = {args.nx * args.ny}", file=sys.stderr)

# Each cylinder sits at its cell center, moved by up to jitter times its clearance to the cell edge
rng = random.Random(args.seed)
reach = args.jitter * (0.5 * L - a)
cylinders = {}
for k, (i, j) in enumerate(((i, j) for j in range(args.ny) for i in range(args.nx)), start=1):
    cylinders.update(
        {
            f"patch_ib({k})%geometry": 2,
            f"patch_ib({k})%x_centroid": (i + 0.5) * L + rng.uniform(-reach, reach),
            f"patch_ib({k})%y_centroid": (j + 0.5) * L + rng.uniform(-reach, reach),
            f"patch_ib({k})%radius": a,
            f"patch_ib({k})%slip": "F",
        }
    )

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": args.nx * L,
            "y_domain%beg": 0.0,
            "y_domain%end": args.ny * L,
            "m": args.nx * args.ppc - 1,
            "n": args.ny * args.ppc - 1,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": args.cfl,
            "n_start": 0,
            "t_stop": args.tstop * rho * L**2 / mu,
            "t_save": args.tstop * rho * L**2 / mu / args.saves,
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
            "patch_icpp(1)%x_centroid": 0.5 * args.nx * L,
            "patch_icpp(1)%y_centroid": 0.5 * args.ny * L,
            "patch_icpp(1)%length_x": args.nx * L,
            "patch_icpp(1)%length_y": args.ny * L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho,
            "patch_icpp(1)%alpha(1)": 1.0,
            # The cylinders: no-slip stationary immersed boundaries
            "ib": "T",
            "num_ibs": args.nx * args.ny,
            "fd_order": 2,
            **cylinders,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%pi_inf": gamma * (rho * c**2 / gamma - p0) / (gamma - 1.0),
        }
    )
)
