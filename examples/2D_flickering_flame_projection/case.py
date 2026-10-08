#!/usr/bin/env python3
"""
2D flickering buoyant jet diffusion flame with the all-Mach pressure projection: diluted hydrogen (H2:N2 = 1:1 by volume) issues
at --U from a --d wide slot into a slow air coflow, burns, and the buoyant hot gas puffs at the flicker frequency ~ 0.5 sqrt(g/d)
(Cetegen & Ahmed 1993), ~11 Hz for the default 2 cm slot. Flicker needs buoyancy to dominate the jet's momentum, a Richardson
number (drho/rho) g d/U^2 of order one or more (~4 at the defaults), and a weak coflow, which otherwise damps the shear layer.

The bottom is a Dirichlet (-17) inlet throughout: fuel in the slot, air coflow elsewhere, both ramped from rest over --ramp. The
sides are slip walls and the top a pressure outlet at the air's hydrostatic pressure there; gravity acts in -y. The air starts
hydrostatic at rest, with a hot pocket above the slot that lights the arriving fuel. Unity-Lewis transport by default (--tm 1:
mixture-averaged), viscous, operator-split reactions (h2o2.yaml). SI units.

The flow is buoyancy-driven at ~1-3 m/s, a Mach number of ~3e-3 in the hot gas: the projection steps near the flow's pace,
~1e-4 s, where the explicit solver must resolve the 900 m/s sound speed of the products at ~3e-7 s.
"""

import argparse
import json

parser = argparse.ArgumentParser(description="2D flickering buoyant hydrogen jet diffusion flame, all-Mach pressure projection")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--d", type=float, default=0.02, help="fuel slot width [m] (default: %(default)s)")
parser.add_argument("--U", type=float, default=0.2, help="fuel velocity [m/s] (default: %(default)s)")
parser.add_argument("--coflow", type=float, default=0.05, help="air coflow velocity [m/s] (default: %(default)s)")
parser.add_argument("--ramp", type=float, default=0.02, help="inflow ramp from rest [s] (default: %(default)s)")
parser.add_argument("--width", type=float, default=0.30, help="domain width [m]; the slip-walled sides must leave the plume air to entrain (default: %(default)s)")
parser.add_argument("--height", type=float, default=0.30, help="domain height [m] (default: %(default)s)")
parser.add_argument("--dx", type=float, default=2.0e-4, help="cell size [m] (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=0.5, help="end time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=125, help="number of outputs (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.5, help="cfl_target (default: %(default)s)")
parser.add_argument("--tm", type=int, default=2, help="chem_params%%transport_model: 2 unity-Lewis, 1 mixture-averaged (default: %(default)s)")
args = parser.parse_args()

p0, Ru, g, T0, T_hot = 101325.0, 8314.46261815324, 9.81, 300.0, 1800.0
W, H, dx = args.width, args.height, args.dx
# Species of h2o2.yaml: H2 (1), O2 (4), N2 (10); fuel H2:N2 = 1:1 by volume, air O2:N2 = 0.233:0.767 by mass
Yf = {1: 2.016 / (2.016 + 28.014), 10: 28.014 / (2.016 + 28.014)}
Ya = {4: 0.233, 10: 0.767}
Wf = 1.0 / (Yf[1] / 2.016 + Yf[10] / 28.014)
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
    ({"geometry": 3, "x_centroid": W / 2, "y_centroid": H / 2, "length_x": W, "length_y": H}, gas(Ya, Wa, T0, 0.0)),
    # 2: a hot air pocket from just above the inlet, so the fuel issues into it and ignites
    ({"geometry": 3, "x_centroid": W / 2, "y_centroid": 0.5 * (dx + 0.03), "length_x": 0.03, "length_y": 0.03 - dx}, gas(Ya, Wa, T_hot, 0.0)),
    # 3, 4: the inflow states in the bottom cells, which the Dirichlet inlet holds: air coflow, then fuel in the slot
    ({"geometry": 3, "x_centroid": W / 2, "y_centroid": 0.0, "length_x": W, "length_y": 2 * dx}, gas(Ya, Wa, T0, args.coflow)),
    ({"geometry": 3, "x_centroid": W / 2, "y_centroid": 0.0, "length_x": args.d, "length_y": 2 * dx}, gas(Yf, Wf, T0, args.U)),
]

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": W,
    "y_domain%beg": 0.0,
    "y_domain%end": H,
    "m": int(round(W / dx)) - 1,
    "n": int(round(H / dx)) - 1,
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
    "proj_max_acfl": 1.0e4,
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
}
for k, (geo, state) in enumerate(patches, 1):
    case.update({f"patch_icpp({k})%{key}": val for key, val in {**geo, **state}.items()})
    if k > 1:
        case[f"patch_icpp({k})%alter_patch(1)"] = "T"
case["patch_icpp(4)%alter_patch(3)"] = "T"

print(json.dumps(case))
