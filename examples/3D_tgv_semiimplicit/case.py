#!/usr/bin/env python3
"""
3D Taylor-Green vortex (Re = 1600, Mach 0.1) run with the semi-implicit
pressure projection method. The initial condition is the hard-coded TGV patch
(hcid 380, shared with 3D_TaylorGreenVortex), which pins p0 = 101325, L = 1,
and V0 = 0.1*c0 on a [-pi, pi]^3 periodic box. The adaptive advective time
step is ~1/Mach = 10x the acoustic limit the explicit twin runs at.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="3D Taylor-Green vortex with the semi-implicit projection method")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("--res", type=float, default=1.0, help="Resolution multiplier on the 256^3 base grid (default: 1)")
args = parser.parse_args()

N = round(256 * args.res)

Re = 1600
L = 1
P0 = 101325
C0 = math.sqrt(1.4 * P0)
V0 = 0.1 * C0  # matches Mach = 0.1 hard-coded in the hcid 380 patch
mu = V0 * L / Re

tC = L / V0
t_stop = 20 * tC

case = {
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": -math.pi * L,
    "x_domain%end": math.pi * L,
    "y_domain%beg": -math.pi * L,
    "y_domain%end": math.pi * L,
    "z_domain%beg": -math.pi * L,
    "z_domain%end": math.pi * L,
    "m": N - 1,
    "n": N - 1,
    "p": N - 1,
    "cfl_adap_dt": "T",
    "cfl_target": 0.5,
    "n_start": 0,
    "t_save": tC,
    "t_stop": t_stop,
    # Simulation Algorithm Parameters
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "mixture_err": "T",
    "time_stepper": "rk3",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    "viscous": "T",
    # Projection method
    "proj_method": "T",
    "proj_iter_solver": 4,
    "proj_tol_rel": 1.0e-8,
    "proj_max_iters": 500,
    "proj_cfl_ac": 10.0,
    # Formatted Database Files Structure Parameters
    "format": "silo",
    "precision": "double",
    "omega_wrt(1)": "T",
    "omega_wrt(2)": "T",
    "omega_wrt(3)": "T",
    "qm_wrt": "T",
    "fd_order": 4,
    "parallel_io": "T",
    # Patch 1: hard-coded Taylor-Green vortex initial condition
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0,
    "patch_icpp(1)%y_centroid": 0,
    "patch_icpp(1)%z_centroid": 0,
    "patch_icpp(1)%length_x": 2 * math.pi * L,
    "patch_icpp(1)%length_y": 2 * math.pi * L,
    "patch_icpp(1)%length_z": 2 * math.pi * L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": 0.0,
    "patch_icpp(1)%hcid": 380,
    "patch_icpp(1)%alpha_rho(1)": 1,
    "patch_icpp(1)%alpha(1)": 1,
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (1.4 - 1),
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1 / mu,
}

print(json.dumps(case))
