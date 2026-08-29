#!/usr/bin/env python3
"""2D Laplace pressure jump, semi-implicit projection: water drop (R=1mm) in
air, sigma=0.0728 -> dp = 72.8 Pa. Tanh-smoothed IC with the jump pre-imposed;
run 0.2 capillary times and check the jump holds. Mirrors SemiImplicitFV
cases/2D_laplace_pressure_jump (which is explicit there)."""

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
args = parser.parse_args()

L = 6.0e-3
N = 100
R = 1.0e-3
xc = yc = 3.0e-3
sigma = 0.0728
eps_i = 1.8e-4  # interface smoothing width
rho_w, rho_a = 1000.0, 1.225
p0, dp = 1.0e5, sigma / R
g_w, pi_w = 4.4, 6.0e8
g_a = 1.4
T_end = 7.41e-4  # 0.2 capillary times

H = f"(0.5*(1.0 - tanh((sqrt((x - {xc})**2.0 + (y - {yc})**2.0) - {R})/{eps_i})))"

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
    "t_save": T_end / 2,
    "t_stop": T_end,
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 2,
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
    "surface_tension": "T",
    "sigma": sigma,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    "proj_method": "T",
    "proj_iter_solver": 1,
    "proj_tol": 1.0e-9,
    "proj_max_iters": 300,
    "proj_cfl_ac": 20.0,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": xc,
    "patch_icpp(1)%y_centroid": yc,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": f"{p0} + {dp}*{H}",
    "patch_icpp(1)%alpha(1)": f"max(1.0e-8, {H})",
    "patch_icpp(1)%alpha(2)": f"max(1.0e-8, 1.0 - {H})",
    "patch_icpp(1)%alpha_rho(1)": f"max(1.0e-8, {H})*{rho_w}",
    "patch_icpp(1)%alpha_rho(2)": f"max(1.0e-8, 1.0 - {H})*{rho_a}",
    "patch_icpp(1)%cf_val": H,
    "fluid_pp(1)%gamma": 1.0 / (g_w - 1.0),
    "fluid_pp(1)%pi_inf": g_w * pi_w / (g_w - 1.0),
    "fluid_pp(2)%gamma": 1.0 / (g_a - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
}

print(json.dumps(case))
