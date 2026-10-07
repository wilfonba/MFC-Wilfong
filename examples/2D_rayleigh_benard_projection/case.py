#!/usr/bin/env python3
"""
2D Rayleigh-Benard convection: air between a hot lower plate and a cold upper plate, with the all-Mach pressure projection.

The plates are no-slip isothermal walls (bc_y%isothermal_in/out) at T0 +/- dT/2, periodic in x over --aspect plate gaps. --Ra and
--Pr set dT and the conductivity, with air's viscosity. Heat
crosses the gap by Fourier conduction (fluid_pp(1)%k_therm), which under the projection enters the pressure equation as
dilatation; buoyancy is the gravity body force on the density the heating lowers. The air starts in the conduction profile,
hydrostatic, with a small temperature perturbation that seeds the rolls. SI units.

Ra = g dT H^3 / (T0 nu kappa) (onset at 1708). At the defaults, Ra = 1e5 and air's Pr = 0.71, dT = 8.5 K and the free-fall
velocity is ~0.12 m/s, a Mach number of ~3e-4: the explicit solver would step ~3e7 times at the acoustic limit where the projection
takes ~1.5e4.
"""

import argparse
import json

parser = argparse.ArgumentParser(description="2D Rayleigh-Benard convection of air, all-Mach pressure projection")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--H", type=float, default=0.05, help="plate gap [m] (default: %(default)s)")
parser.add_argument("--Ra", type=float, default=1.0e5, help="Rayleigh number (default: %(default)s)")
parser.add_argument("--Pr", type=float, default=0.71, help="Prandtl number (default: %(default)s, air)")
parser.add_argument("--aspect", type=float, default=2.0, help="domain width in plate gaps (default: %(default)s)")
parser.add_argument("--ny", type=int, default=128, help="cells across the gap (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=20.0, help="end time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=40, help="number of outputs (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.5, help="cfl_target (default: %(default)s)")
args = parser.parse_args()

# Air at 300 K and 1 atm, its conductivity set by --Pr
T0, p0, gam, cv, mu = 300.0, 101325.0, 1.4, 718.0, 1.85e-5
R = cv * (gam - 1)
rho0 = p0 / (R * T0)
g = 9.81
H, L = args.H, args.aspect * args.H
nx = int(round(args.aspect * args.ny))
nu = mu / rho0
kappa = nu / args.Pr
k = rho0 * gam * cv * kappa
dT = args.Ra * T0 * nu * kappa / (g * H**3)
print(f"dT = {dT:.3g} K, k = {k:.3g} W/m/K, {nx} x {args.ny} cells", file=__import__("sys").stderr)

# Conduction profile with a two-mode perturbation of 1% dT; hydrostatic pressure of the mean density
T = f"({T0 + dT / 2!r} - {dT / H!r}*y + {0.01 * dT!r}*sin(pi*y/{H!r})*(cos(2*pi*x/{L!r}) + 0.5*cos(4*pi*x/{L!r} + 1.0)))"
pres = f"({p0!r} - {rho0 * g!r}*y)"

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": H,
            "m": nx - 1,
            "n": args.ny - 1,
            "p": 0,
            "cfl_adap_dt": "T",
            "cfl_target": args.cfl,
            "n_start": 0,
            "t_stop": args.tstop,
            "t_save": args.tstop / args.saves,
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
            "bc_y%beg": -16,
            "bc_y%end": -16,
            "bc_y%isothermal_in": "T",
            "bc_y%Twall_in": T0 + dT / 2,
            "bc_y%isothermal_out": "T",
            "bc_y%Twall_out": T0 - dT / 2,
            "proj_method": "T",
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu,
            "fluid_pp(1)%cv": cv,
            "fluid_pp(1)%k_therm": k,
            "bf_y": "T",
            "g_y": -g,
            "k_y": 0.0,
            "w_y": 0.0,
            "p_y": 0.0,
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "T_wrt": "T",
            "parallel_io": "T",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": L / 2,
            "patch_icpp(1)%y_centroid": H / 2,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": H,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": pres,
            "patch_icpp(1)%alpha_rho(1)": f"{pres}/({R!r}*{T})",
            "patch_icpp(1)%alpha(1)": 1.0,
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
        }
    )
)
