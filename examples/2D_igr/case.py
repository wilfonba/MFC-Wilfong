#!/usr/bin/env python3
# This case file demonstrates the Laplace pressure jump of a water droplet in air. The laplace pressure jump
# in 2D is given by delta = sigma / r where delta is the pressure jump, sigma is the surface tension coefficient,
# and r is the radius of the droplet. The results of this simulation agree with theory to well within 1%
# relative error.

import math
import json

l = 1

# Numerical setup
r0 = l/4
x0 = -l
x1 = l
y0 = -l
y1 = l

Nx = 199
Ny = 199

mydt = 2.5e-3

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "cfl_const_dt": "T",
    "cfl_target": 0.1,
    "n_start": 0,
    "t_stop": 2.5,
    "t_save": 2.5/100,
    # "dt": mydt,
    # "t_step_start": 0,
    # "t_step_stop": 1000,
    # "t_step_save": 10,
    # Simulation Algorithm
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "avg_state": 2,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "num_patches": 2,
    "num_fluids": 1,
    "igr": "T",
    "viscous": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Fluid Parameters (Gas)
    "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
    "fluid_pp(1)%pi_inf": 0.0e00,
    "fluid_pp(1)%Re(1)": 1e5,
    # Ambient pressure
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": 0,
    "patch_icpp(1)%length_x": 2,
    "patch_icpp(1)%length_y": 2,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": 1,
    "patch_icpp(1)%alpha_rho(1)": 1,
    "patch_icpp(1)%alpha(1)": 1,
    # High pressure
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.9,
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0,
    "patch_icpp(2)%y_centroid": 0,
    "patch_icpp(2)%radius": r0,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": 5,
    "patch_icpp(2)%alpha_rho(1)": 5,
    "patch_icpp(2)%alpha(1)": 1,
}

print(json.dumps(data))
