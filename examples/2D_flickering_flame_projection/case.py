#!/usr/bin/env python3
"""
2D flickering buoyant jet diffusion flame with the all-Mach pressure projection: pure hydrogen issues at --U from a --d wide slot
into a slow air coflow, burns, and the buoyant hot gas puffs at the flicker frequency ~ 0.5 sqrt(g/d) (Cetegen & Ahmed 1993), ~7 Hz
for the default 5 cm slot. Flicker needs buoyancy to dominate the jet's momentum, a Richardson number (drho/rho) g d/U^2 of order
one or more (~10 at the defaults), a weak coflow, which otherwise damps the shear layer, and enough heat release and width to
outgrow viscosity (Grashof number ~ g d^3/nu^2): a narrow or diluted flame burns steadily, as a candle does.

The bottom is a Dirichlet (-17) inlet throughout: fuel in the slot, air coflow elsewhere, both ramped from rest over --ramp. The
sides are slip walls and the top a pressure outlet at the air's hydrostatic pressure there; gravity acts in -y. The slot sits at
x = 0 in a 0.30 m by 0.30 m core of --dx cells, and the grid stretches (tanh, pre_process's stretch_x/y) beyond it to a --width
by --height domain, moving the walls and the outlet away from the flame at a fraction of a uniform grid's cells. The air starts
hydrostatic at rest, with a hot pocket above the slot that lights the arriving fuel. Unity-Lewis transport by default (--tm 1:
mixture-averaged), viscous, operator-split reactions (h2o2.yaml). SI units.

The flow is buoyancy-driven at ~1-3 m/s, a Mach number of ~3e-3 in the hot gas: the projection steps near the flow's pace,
~1e-4 s, where the explicit solver must resolve the 900 m/s sound speed of the products at ~3e-7 s.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="2D flickering buoyant hydrogen jet diffusion flame, all-Mach pressure projection")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--d", type=float, default=0.05, help="fuel slot width [m] (default: %(default)s)")
parser.add_argument("--U", type=float, default=0.2, help="fuel velocity [m/s] (default: %(default)s)")
parser.add_argument("--coflow", type=float, default=0.05, help="air coflow velocity [m/s] (default: %(default)s)")
parser.add_argument("--ramp", type=float, default=0.02, help="inflow ramp from rest [s] (default: %(default)s)")
parser.add_argument("--width", type=float, default=0.60, help="domain width [m]; the slip-walled sides must leave the plume air to entrain (default: %(default)s)")
parser.add_argument("--height", type=float, default=0.90, help="domain height [m] (default: %(default)s)")
parser.add_argument("--dx", type=float, default=2.0e-4, help="cell size in the core [m] (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=1.0, help="end time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=200, help="number of outputs (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.5, help="cfl_target (default: %(default)s)")
parser.add_argument("--tm", type=int, default=2, help="chem_params%%transport_model: 2 unity-Lewis, 1 mixture-averaged (default: %(default)s)")
args = parser.parse_args()

p0, Ru, g, T0, T_hot = 101325.0, 8314.46261815324, 9.81, 300.0, 1800.0
dx = args.dx
# Stretching: the --dx core's width and height, a_x/a_y (sharpness of its edge) and loops_x/loops_y (passes of the map). At the
# defaults the cells grow at most 1.5% per cell, to ~5 dx at the sides and ~15 dx at the top, in 1842 x 1776 cells
core_w, core_h, a, loops = 0.30, 0.30, 20.0, 2


def lc(v):
    return math.log(math.cosh(v))


def stretched_edges(beg, end, N, c_a, c_b):
    """Cell edges of N cells on [beg, end] after pre_process's tanh stretching about [c_a, c_b] (m_grid)"""
    L = end - beg
    x = [(beg + L * i / N) / L for i in range(N + 1)]
    ca, cb = c_a / L, c_b / L
    for _ in range(loops):
        x = [v / a * (a + lc(a * (v - ca)) + lc(a * (v - cb)) - 2 * lc(a * (cb - ca) / 2)) for v in x]
    return [v * L for v in x]


def axis(total, core, symmetric):
    """Domain bounds, cell count and stretching keys of one direction: a core of dx cells (about 0, or from 0 up) stretched to
    total. The map changes the extent, so the unstretched extent is found by bisection; returns the stretched edges too"""
    c_a, c_b = (-core / 2, core / 2) if symmetric else (-core, core)

    def bounds(Ln):
        return (-Ln / 2, Ln / 2) if symmetric else (0.0, Ln)

    if total <= core:
        beg, end = bounds(core)
        N = int(round(core / dx))
        return beg, end, N, {}, [beg + (end - beg) * i / N for i in range(N + 1)]
    lo, hi = core, total
    for _ in range(60):
        Ln = 0.5 * (lo + hi)
        e = stretched_edges(*bounds(Ln), int(round(Ln / dx)), c_a, c_b)
        lo, hi = (Ln, hi) if e[-1] - e[0] < total else (lo, Ln)
    beg, end = bounds(lo)
    N = int(round(lo / dx))
    return beg, end, N, {"a": a, "_a": c_a, "_b": c_b, "loops": loops}, stretched_edges(beg, end, N, c_a, c_b)


x_beg, x_end, Nx, sx, _ = axis(args.width, core_w, True)
y_beg, y_end, Ny, sy, ye = axis(args.height, core_h, False)
W, H = args.width, ye[-1]  # H: the outlet's height
# Species of h2o2.yaml: H2 (1), O2 (4), N2 (10); fuel pure H2, air O2:N2 = 0.233:0.767 by mass
Yf = {1: 1.0}
Ya = {4: 0.233, 10: 0.767}
Wf = 2.016
Wa = 1.0 / (Ya[4] / 31.998 + Ya[10] / 28.014)
rho_a = p0 * Wa / (Ru * T0)
pres = f"({p0!r} - {rho_a * g!r}*y)"  # the air's hydrostatic pressure


def gas(Y, W_mix, T, vel2):
    """Primitive state of a patch at temperature T with mass fractions Y"""
    out = {"pres": pres, "alpha(1)": 1.0, "alpha_rho(1)": f"{pres}*{W_mix / (Ru * T)!r}", "vel(1)": 0.0, "vel(2)": vel2}
    out.update({f"Y({i})": Y.get(i, 0.0) for i in range(1, 11)})
    return out


patches = [
    # 1: still air, hydrostatic
    ({"geometry": 3, "x_centroid": 0.0, "y_centroid": H / 2, "length_x": 2 * W, "length_y": 2 * H}, gas(Ya, Wa, T0, 0.0)),
    # 2: a hot air pocket from just above the inlet, past the slot edges, so the fuel issues into it and ignites
    ({"geometry": 3, "x_centroid": 0.0, "y_centroid": 0.5 * (dx + 0.03), "length_x": args.d + 0.02, "length_y": 0.03 - dx}, gas(Ya, Wa, T_hot, 0.0)),
    # 3, 4: the inflow states in the bottom cells, which the Dirichlet inlet holds: air coflow, then fuel in the slot
    ({"geometry": 3, "x_centroid": 0.0, "y_centroid": 0.0, "length_x": 2 * W, "length_y": 2 * dx}, gas(Ya, Wa, T0, args.coflow)),
    ({"geometry": 3, "x_centroid": 0.0, "y_centroid": 0.0, "length_x": args.d, "length_y": 2 * dx}, gas(Yf, Wf, T0, args.U)),
]

case = {
    "run_time_info": "T",
    "x_domain%beg": x_beg,
    "x_domain%end": x_end,
    "y_domain%beg": y_beg,
    "y_domain%end": y_end,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cfl_adap_dt": "T",
    "cfl_target": args.cfl,
    "n_start": 0,
    "t_stop": args.tstop,
    "t_save": args.tstop / args.saves,
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": len(patches),
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -15,
    "bc_x%end": -15,
    "bc_y%beg": -17,
    "bc_y%end": -3,
    "bc_y%vel_in_ramp": args.ramp,
    "bc_y%pres_out": p0 - rho_a * g * H,
    "proj_method": "T",
    "proj_max_acfl": 250,
    "bf_y": "T",
    "g_y": -g,
    "k_y": 0.0,
    "w_y": 0.0,
    "p_y": 0.0,
    "viscous": "T",
    "fluid_pp(1)%Re(1)": 1.0e5,  # unused: a reacting mixture takes its viscosity from the transport model
    "chemistry": "T",
    "cantera_file": "h2o2.yaml",
    "chem_params%diffusion": "T",
    "chem_params%transport_model": args.tm,
    "chem_params%reactions": "T",
    "chem_params%reaction_substeps": 4,
    "chem_params%adap_substeps": "T",
    "chem_params%reaction_substeps_max": 64,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "chem_wrt_T": "T",
    "diff_sts": "T",
    "diff_sts_max": 64,
}
for xy, st in [("x", sx), ("y", sy)]:
    if st:
        case.update({f"stretch_{xy}": "T", f"a_{xy}": st["a"], f"{xy}_a": st["_a"], f"{xy}_b": st["_b"], f"loops_{xy}": st["loops"]})
for k, (geo, state) in enumerate(patches, 1):
    case.update({f"patch_icpp({k})%{key}": val for key, val in {**geo, **state}.items()})
    if k > 1:
        case[f"patch_icpp({k})%alter_patch(1)"] = "T"
case["patch_icpp(4)%alter_patch(3)"] = "T"

print(json.dumps(case))
