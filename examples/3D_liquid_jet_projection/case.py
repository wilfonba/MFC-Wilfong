#!/usr/bin/env python3
"""
3D round liquid jet: water issuing from a nozzle in a wall into still air, with the all-Mach pressure projection.

The nozzle is a circular Dirichlet (-17) boundary patch of diameter --d set into the no-slip floor; its inflow velocity ramps from rest to
--U over --ramp (bc_y%vel_in_ramp), so the start-up is part of the flow rather than an impulsive slug. The four sides and the top
are pressure outlets at the ambient pressure (extrapolation boundaries with bc_[x,y,z]%pres_out), through which the jet and the air it
entrains leave. Water and air at 20 C: water a stiffened gas (gamma = 6.12, pi_inf = 3.43e8 Pa), air an ideal gas, with
viscosity and surface tension. SI units.

The Mach number is ~3e-3 in the water and ~1.5e-2 in the air, so the projection steps at the flow's pace where the explicit
solver must resolve the water's 1450 m/s sound speed: at the defaults 1,960 steps against ~350,000 (the air around the jet's
head sets the projection's step). --explicit runs it as a control.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="3D round liquid jet into air, all-Mach pressure projection")
parser.add_argument("--d", type=float, default=2.0e-3, help="nozzle diameter [m] (default: %(default)s)")
parser.add_argument("--U", type=float, default=5.0, help="jet velocity [m/s] (default: %(default)s)")
parser.add_argument("--ramp", type=float, default=2.0e-3, help="duration of the inflow ramp from rest [s] (default: %(default)s)")
parser.add_argument("--width", type=float, default=12, help="domain width (x and z) in nozzle diameters (default: %(default)s)")
parser.add_argument("--height", type=float, default=30, help="domain height in nozzle diameters (default: %(default)s)")
parser.add_argument("--ppd", type=int, default=10, help="cells per nozzle diameter (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=0.012, help="end time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=40, help="number of outputs (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.25, help="cfl_target: advective with the projection, acoustic with --explicit (default: %(default)s)")
parser.add_argument("--profile", choices=["tophat", "turbulent"], default="tophat", help="nozzle exit profile: uniform, or the 1/7-power pipe profile of the same flux (default: %(default)s)")
parser.add_argument("--noise", type=float, default=0.0, help="simplex noise on the inflow velocity, as a fraction of the local speed (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
args, _ = parser.parse_known_args()

# Water and air at 20 C
rho_l, mu_l, rho_g, mu_g, sigma, p0 = 998.0, 1.0e-3, 1.2, 1.8e-5, 0.072, 101325.0
gamma_l, pi_inf_l, gamma_g = 6.12, 3.43e8, 1.4

D, U = args.d, args.U
W, H = args.width * D, args.height * D
dx = D / args.ppd
eps = 1.0e-8
# 1/7-power pipe profile; its mean is 98/120 of its centerline speed
r = f"sqrt((x - {0.5 * W!r})**2 + (z - {0.5 * W!r})**2) / {0.5 * D!r}"
vel_in = U if args.profile == "tophat" else f"{U * 120 / 98!r} * max(1.0 - {r}, 0.0)**{1 / 7!r}"
# Frozen noise on the nozzle's boundary cells, which the Dirichlet inflow is read from: features ~ d/2, decorrelated components
noise = {"simplex_perturb": "T"}
for i, o in enumerate([(12.3, -11.3, 34.6), (-70.3, 33.4, -34.6), (123.3, -654.3, -64.5)], 1):
    noise.update({f"simplex_params%perturb_vel({i})": "T", f"simplex_params%perturb_vel_freq({i})": 2 / D, f"simplex_params%perturb_vel_scale({i})": args.noise})
    noise.update({f"simplex_params%perturb_vel_offset({i},{j})": v * D for j, v in enumerate(o, 1)})
print(
    f"We = {rho_l * U**2 * D / sigma:.0f}, Re = {rho_l * U * D / mu_l:.0f}, Mach (water) = {U / math.sqrt(gamma_l * (p0 + pi_inf_l) / rho_l):.1e}, "
    f"{int(round(W / dx))} x {int(round(H / dx))} x {int(round(W / dx))} cells",
    file=__import__("sys").stderr,
)

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": W,
            "y_domain%beg": 0.0,
            "y_domain%end": H,
            "z_domain%beg": 0.0,
            "z_domain%end": W,
            "m": int(round(W / dx)) - 1,
            "n": int(round(H / dx)) - 1,
            "p": int(round(W / dx)) - 1,
            "cfl_adap_dt": "T",
            "cfl_target": args.cfl,
            "n_start": 0,
            "t_stop": args.tstop,
            "t_save": args.tstop / args.saves,
            "num_patches": 2,
            "model_eqns": 2,
            "num_fluids": 2,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            # Open sides and top at the ambient pressure; a no-slip floor holding the nozzle
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -16,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            "bc_x%pres_out": p0,
            "bc_y%pres_out": p0,
            "bc_z%pres_out": p0,
            "num_bc_patches": 1,
            "patch_bc(1)%dir": 2,
            "patch_bc(1)%loc": -1,
            "patch_bc(1)%geometry": 2,
            "patch_bc(1)%type": -17,
            "patch_bc(1)%centroid(1)": 0.5 * W,
            "patch_bc(1)%centroid(3)": 0.5 * W,
            "patch_bc(1)%radius": 0.5 * D,
            "bc_y%vel_in_ramp": args.ramp,
            "proj_method": "F" if args.explicit else "T",
            "proj_max_acfl": 0.0 if args.explicit else 1000.0,  # bounds only the start-up, while the air is at rest
            # MTHINC, softened, as in the rising-bubble benchmark
            "int_comp": 2,
            "ic_beta": 0.6,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 1.0 / mu_l,
            "fluid_pp(2)%Re(1)": 1.0 / mu_g,
            "surface_tension": "T",
            "sigma": sigma,
            "surface_tension_model": "conservative" if args.explicit else "well_balanced",
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            # Patch 1: still air
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.5 * W,
            "patch_icpp(1)%y_centroid": 0.5 * H,
            "patch_icpp(1)%z_centroid": 0.5 * W,
            "patch_icpp(1)%length_x": W,
            "patch_icpp(1)%length_y": H,
            "patch_icpp(1)%length_z": W,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
            "patch_icpp(1)%alpha_rho(2)": (1 - eps) * rho_g,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": 1 - eps,
            "patch_icpp(1)%cf_val": 0,
            # Patch 2: the water in the nozzle's boundary cells, which the Dirichlet patch takes its inflow state from
            "patch_icpp(2)%geometry": 10,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.5 * W,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.5 * W,
            "patch_icpp(2)%radius": 0.5 * D,
            "patch_icpp(2)%length_y": 2 * dx,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": vel_in,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%alpha_rho(1)": (1 - eps) * rho_l,
            "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
            "patch_icpp(2)%alpha(1)": 1 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            "patch_icpp(2)%cf_val": 1,
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
            "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
            "fluid_pp(2)%eos": "stiffened_gas",
            "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            **(noise if args.noise > 0 else {}),
        }
    )
)
