#!/usr/bin/env python3
"""
1D hydrostatic water column with the semi-implicit projection method:
stiffened-gas water under gravity between no-slip walls, starting from a
uniform (non-equilibrium) state. The pressure solve balances gravity, so
the column settles to hydrostatic equilibrium with near-zero velocity at a
fixed dt ~80x the acoustic limit. Mirrors SemiImplicitFV
cases/1D_hydrostatic_water.
"""

import json

gamma = 4.4
p_inf = 6.0e8  # stiffened-gas reference pressure of water
rho0 = 1000.0
p0 = 1.0e5
g = -9.81

N = 200
L = 10.0
dt = 1.0e-3  # quiescent flow: fixed dt (advective CFL is unbounded at u = 0)
T_end = 1.6
Nt = int(T_end / dt)

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "m": N - 1,
    "n": 0,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Nt,
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "F",
    "time_stepper": "rk3",
    "recon_type": "weno",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "F",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    "bc_x%beg": -16,
    "bc_x%end": -16,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    "proj_method": "T",
    "proj_iter_solver": 1,
    "proj_tol": 1.0e-9,
    "proj_max_iters": 200,
    "bf_x": "T",
    "g_x": g,
    "k_x": 0.0,
    "w_x": 0.0,
    "p_x": 0.0,
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": L / 2,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": gamma * p_inf / (gamma - 1.0),
}

print(json.dumps(case))
