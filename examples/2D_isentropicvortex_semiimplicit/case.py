#!/usr/bin/env python3
"""
2D isentropic vortex (Shu 1998) advected with the semi-implicit pressure
projection method at Mach 0.1: background p = 2/(gamma*Ma^2) makes the sound
speed ~14 while the advection speed is ~sqrt(2), so the adaptive advective
time step is ~10x the acoustic limit an explicit scheme would need. One
diagonal period (t = 10) returns the vortex to its initial position.
Mirrors SemiImplicitFV cases/2D_isentropic_vortex.
"""

import json
import math

gamma = 1.4
beta = 5.0  # vortex strength
T_inf = 142.857142857143  # p_inf/(rho_inf*R) with p_inf = 2/(gamma*Ma^2), Ma = 0.1
xc, yc = 5.0, 5.0

N = 128
L = 10.0

# exp((1 - r^2)/2) with r^2 about the vortex center
f = f"exp((1.0 - (x - {xc})**2.0 - (y - {yc})**2.0)/2.0)"
# dT = -(gamma - 1)*beta^2/(8*gamma*pi^2) * f^2
dT = f"(-{(gamma - 1.0)*beta**2}/(8.0*{gamma}*pi**2.0)*{f}**2.0)"
T = f"({T_inf} + {dT})"
rho = f"({T}/{T_inf})**{1.0/(gamma - 1.0)}"

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "y_domain%beg": 0.0,
    "y_domain%end": L,
    "m": N - 1,
    "n": N - 1,
    "p": 0,
    "cfl_adap_dt": "T",
    "cfl_target": 0.5,
    "n_start": 0,
    "t_save": 1.0,
    "t_stop": 10.0,
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": "rk3",
    "recon_type": "weno",
    "weno_order": 3,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "proj_method": "T",
    "proj_iter_solver": 1,
    "proj_tol": 1.0e-8,
    "proj_max_iters": 500,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": xc,
    "patch_icpp(1)%y_centroid": yc,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": f"1.0 - {beta}/(2.0*pi)*(y - {yc})*{f}",
    "patch_icpp(1)%vel(2)": f"1.0 + {beta}/(2.0*pi)*(x - {xc})*{f}",
    "patch_icpp(1)%pres": f"{rho}*{T}",
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
}

print(json.dumps(case))
