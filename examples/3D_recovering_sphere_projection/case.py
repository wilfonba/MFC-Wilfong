#!/usr/bin/env python3
"""
A cubic water droplet in air recovering a spherical shape under surface tension, with the all-Mach pressure projection.

The projection counterpart of examples/3D_recovering_sphere: one octant of a cube of side 0.15 in a box of side 0.375,
symmetry planes on the low faces and extrapolation on the high ones, sigma = 8. The flow is driven by surface tension
alone at speeds of order sqrt(sigma/(rho_w R)) ~ 0.3, against sound speeds of ~50 (water) and ~370 (air), so the
projection steps at the capillary limit, two to three orders of magnitude past the explicit acoustic one. --explicit runs
the HLLC solver instead, with the conservative surface tension model.
"""

import argparse
import json

parser = argparse.ArgumentParser(description="3D recovering droplet, all-Mach pressure projection")
parser.add_argument("--N", type=int, default=200, help="cells per direction (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.5, help="cfl_target: advective/capillary with the projection, acoustic with --explicit (default: %(default)s)")
parser.add_argument("--tstop", type=float, default=0.5, help="end time (default: %(default)s)")
parser.add_argument("--saves", type=int, default=50, help="number of outputs (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit instead, as a control")
args, _ = parser.parse_known_args()

L = 0.375
side = 0.15
eps = 1e-9
p0 = 1.0e5

print(
    json.dumps(
        {
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "z_domain%beg": 0.0,
            "z_domain%end": L,
            "m": args.N - 1,
            "n": args.N - 1,
            "p": args.N - 1,
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
            "bc_x%beg": -2,
            "bc_x%end": -3,
            "bc_y%beg": -2,
            "bc_y%end": -3,
            "bc_z%beg": -2,
            "bc_z%end": -3,
            "proj_method": "F" if args.explicit else "T",
            # The fluid starts at rest, so the step is capped until the capillary flow sets the advective limit
            "proj_max_acfl": 500.0,
            # MTHINC, softened as in the other projection examples
            "int_comp": 2,
            "ic_beta": 1.0,
            "surface_tension": "T",
            "sigma": 8.0,
            # The well-balanced model needs the projection's face pressure gradient
            "surface_tension_model": "conservative" if args.explicit else "well_balanced",
            "format": 1,
            "precision": 2,
            "alpha_wrt(1)": "T",
            "cf_wrt": "T",
            "parallel_io": "T",
            # Water
            "fluid_pp(1)%eos": "stiffened_gas",
            "fluid_pp(1)%gamma": 1.0 / (2.1 - 1.0),
            "fluid_pp(1)%pi_inf": 2.1 * 1.0e6 / (2.1 - 1.0),
            # Air
            "fluid_pp(2)%eos": "ideal_gas",
            "fluid_pp(2)%gamma": 1.0 / (1.4 - 1.0),
            # Air everywhere
            "patch_icpp(1)%geometry": 9,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_x": 2 * L,
            "patch_icpp(1)%length_y": 2 * L,
            "patch_icpp(1)%length_z": 2 * L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": eps * 1000,
            "patch_icpp(1)%alpha_rho(2)": (1 - eps) * 1,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": 1 - eps,
            "patch_icpp(1)%cf_val": 0,
            # The water cube, centered on the domain corner so this octant holds an eighth of it
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%geometry": 9,
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%length_x": side,
            "patch_icpp(2)%length_y": side,
            "patch_icpp(2)%length_z": side,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
            "patch_icpp(2)%pres": p0,
            "patch_icpp(2)%alpha_rho(1)": (1 - eps) * 1000,
            "patch_icpp(2)%alpha_rho(2)": eps * 1,
            "patch_icpp(2)%alpha(1)": 1 - eps,
            "patch_icpp(2)%alpha(2)": eps,
            "patch_icpp(2)%cf_val": 1,
        }
    )
)
