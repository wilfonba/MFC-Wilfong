#!/usr/bin/env python3
"""
2D bubble column: a swarm of air bubbles rising through a water column under an air headspace, with the all-Mach pressure
projection.

The bubbles start at rest on a jittered lattice near the bottom of a column of water of depth H, capped by an air
headspace at pressure p_top. Gravity drives them up; each expands as the hydrostatic pressure p_top + rho_l g (H - y)
falls, by p_bottom/p_top over the full depth, which changes its buoyancy, shape and wake on the way up and is what an
incompressible bubbly-flow solver leaves out. They burst into the headspace at the top. At atmospheric p_top that ratio
is 1.03 at the default 0.3 m depth and 2 at 10 m: deep columns, the case for scale, are where it matters.

Material properties are those of water and air at 20 C: water a stiffened gas (gamma = 6.12, pi_inf = 3.43e8 Pa, Le
Metayer et al. 2004; sound speed ~1450 m/s), air an ideal gas, so rho_l/rho_g ~ 830 at 1 atm. The flow Mach number in
the water is ~1e-4, so the projection steps at the flow's and the capillary pace; --explicit runs HLLC at the
acoustic one, as a control. SI units. The gas expands adiabatically (no heat conduction), where millimetre bubbles
rising this slowly are close to isothermal; in 2D each bubble is a cylinder.

The lattice holds --cols x --rows bubbles of mean diameter --d at spacing --spacing diameters, each moved and resized
randomly (--seed) by a hash of its lattice cell, so the initial condition is one analytic patch of fixed size.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="2D bubble column with a low-pressure headspace, all-Mach pressure projection")
parser.add_argument("--d", type=float, default=4.0e-3, help="mean bubble diameter [m] (default: %(default)s)")
parser.add_argument("--cols", type=int, default=5, help="bubbles across the column (default: %(default)s)")
parser.add_argument("--rows", type=int, default=6, help="rows of bubbles (default: %(default)s)")
parser.add_argument("--spacing", type=float, default=2.5, help="lattice spacing in diameters; sets the column width (default: %(default)s)")
parser.add_argument("--spread", type=float, default=0.3, help="relative spread of the bubble diameters (default: %(default)s)")
parser.add_argument("--depth", type=float, default=0.3, help="water depth H [m] (default: %(default)s)")
parser.add_argument("--p-top", type=float, default=101325.0, help="headspace pressure [Pa] (default: %(default)s)")
parser.add_argument("--ppd", type=int, default=16, help="cells per mean bubble diameter (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=1.5, help="end time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=75, help="number of outputs (default: %(default)s)")
parser.add_argument("--seed", type=int, default=1, help="seed of the bubble positions and sizes (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="cfl_target: advective with the projection, acoustic with --explicit (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
parser.add_argument("--single", action="store_true", help="one centred bubble of diameter --d in a column --spacing diameters wide, as an expansion check")
args, _ = parser.parse_known_args()
if args.single:
    args.cols, args.rows, args.spread = 1, 1, 0.0

# Water and air at 20 C
rho_l, mu_l, mu_g, sigma = 998.0, 1.0e-3, 1.8e-5, 0.072
gamma_l, pi_inf_l, gamma_g, R_air, T = 6.12, 3.43e8, 1.4, 287.0, 293.0
g = 9.81

R, dx = 0.5 * args.d, args.d / args.ppd
s = args.spacing * args.d
W, H, headspace = args.cols * s, args.depth, 0.1 * args.depth
y_b = s  # bottom of the bubble lattice
w = 0.75 * dx  # half-width of the bubbles' tanh interface
jitter = 0.0 if args.single else max(0.0, 0.5 * s - R * (1 + 0.5 * args.spread) - 3 * dx)  # keeps every bubble inside its lattice cell
if y_b + args.rows * s > 0.5 * H:
    parser.error("the bubble lattice must fit in the lower half of the column: fewer --rows or a larger --depth")


def num(v):
    return f"{v:.7g}"


# Per-cell draw in [0, 1]: sin() of a large multiple of the lattice cell's origin (x0, y0)
x0 = f"(x - mod(x, {num(s)}))"
yy = f"(y - {num(y_b)})"
y0 = f"({yy} - mod({yy}, {num(s)}))"
key = f"({x0}*{num(12.9898 / s)} + {y0}*{num(78.233 / s)} + {num(0.6180339887 * args.seed)})"


def draw(k):
    return f"(0.5 + 0.5*sin({key}*{num(437.585453 * k)}))"


cx = f"({x0} + {num(0.5 * s)} + {num(jitter)}*(2*{draw(1.0)} - 1))" if jitter else f"({x0} + {num(0.5 * s)})"
cy = f"({num(y_b + 0.5 * s)} + {y0} + {num(jitter)}*(2*{draw(1.37)} - 1))" if jitter else f"({num(y_b + 0.5 * s)} + {y0})"
rb = f"({num(R)}*(1 + {num(args.spread)}*({draw(1.91)} - 0.5)))" if args.spread else num(R)
inside = f"0.5*(1 - tanh((sqrt((x - {cx})**2 + (y - {cy})**2) - {rb})/{num(w)}))"  # 1 in a bubble, 0 outside

eps = 1.0e-8
p_l = f"({num(args.p_top)} + {num(rho_l * g)}*({num(H)} - y))"  # hydrostatic pressure of the water
alpha_g = f"({num(eps)} + {num(1 - 2 * eps)}*{inside})"
# Each bubble holds the water's pressure at its lattice cell's centre plus the Laplace jump sigma/R, uniform inside (the jitter
# offsets it by at most rho_l g jitter), and its air at that pressure and T
yl = f"({num(y_b + 0.5 * s)} + {y0})"  # height of the bubble's lattice cell centre, within the jitter of its own
p_c = f"({num(args.p_top + sigma / R)} + {num(rho_l * g)}*({num(H)} - {yl}))"  # the bubble's pressure
p_bubble = f"({p_l} + {inside}*({num(sigma / R)} + {num(rho_l * g)}*(y - {yl})))"
rho_g_bubble = f"{p_c}/{num(R_air * T)}"
rho_g_top = args.p_top / (R_air * T)

p_bottom = args.p_top + rho_l * g * H
print(
    f"{args.cols * args.rows} bubbles, column {W:.3g} x {H:.3g} m, {int(round(W / dx))} x {int(round((H + headspace) / dx))} cells, "
    f"p_bottom/p_top = {p_bottom / args.p_top:.3f}, rho_l/rho_g = {rho_l / rho_g_top:.0f} at the top, "
    f"Eo = {rho_l * g * args.d**2 / sigma:.2f}",
    file=__import__("sys").stderr,
)

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": W,
            "y_domain%beg": 0.0,
            "y_domain%end": H + headspace,
            "m": int(round(W / dx)) - 1,
            "n": int(round((H + headspace) / dx)) - 1,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": args.cfl,
            "n_start": 0,
            "t_stop": args.tstop,
            "t_save": args.tstop / args.saves,
            "num_patches": 3,
            "model_eqns": 2,
            "num_fluids": 2,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            # Slip side walls, no-slip floor and lid
            "bc_x%beg": -15,
            "bc_x%end": -15,
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "proj_method": "F" if args.explicit else "T",
            # MTHINC, softened, as in the rising-bubble benchmark
            "int_comp": 2,
            "ic_beta": 1.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu_l,
            "fluid_pp(2)%Re(1)": 1.0 / mu_g,
            "surface_tension": "T",
            "sigma": sigma,
            "surface_tension_model": "conservative" if args.explicit else "well_balanced",
            "bf_y": "T",
            "g_y": -g,
            "k_y": 0.0,
            "w_y": 0.0,
            "p_y": 0.0,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: the water, hydrostatic
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * W,
            "patch_icpp(1)%y_centroid": 0.5 * H,
            "patch_icpp(1)%length_x": W,
            "patch_icpp(1)%length_y": H,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_l,
            "patch_icpp(1)%alpha_rho(1)": (1 - eps) * rho_l,
            "patch_icpp(1)%alpha_rho(2)": eps * rho_g_top,
            "patch_icpp(1)%alpha(1)": 1 - eps,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(1)%cf_val": 0,
            # Patch 2: the air headspace at p_top
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%x_centroid": 0.5 * W,
            "patch_icpp(2)%y_centroid": H + 0.5 * headspace,
            "patch_icpp(2)%length_x": W,
            "patch_icpp(2)%length_y": headspace,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": args.p_top,
            "patch_icpp(2)%alpha_rho(1)": eps * rho_l,
            "patch_icpp(2)%alpha_rho(2)": (1 - eps) * rho_g_top,
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha(2)": 1 - eps,
            "patch_icpp(2)%cf_val": 1,
            # Patch 3: the bubble lattice, one analytic patch over its rows
            "patch_icpp(3)%geometry": 3,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%x_centroid": 0.5 * W,
            "patch_icpp(3)%y_centroid": y_b + 0.5 * args.rows * s,
            "patch_icpp(3)%length_x": W,
            "patch_icpp(3)%length_y": args.rows * s,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": p_bubble,
            "patch_icpp(3)%alpha_rho(1)": f"(1 - {alpha_g})*{num(rho_l)}",
            "patch_icpp(3)%alpha_rho(2)": f"{alpha_g}*{rho_g_bubble}",
            "patch_icpp(3)%alpha(1)": f"1 - {alpha_g}",
            "patch_icpp(3)%alpha(2)": alpha_g,
            "patch_icpp(3)%cf_val": inside,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
            "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
            "fluid_pp(2)%eos": "stiffened_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
        }
    )
)
