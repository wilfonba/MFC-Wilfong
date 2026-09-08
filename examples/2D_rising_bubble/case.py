#!/usr/bin/env python3
"""
2D rising bubble example. --case selects the configuration:

  1 (default), 2: test cases 1 and 2 of the Hysing et al. (2009) benchmark, IJNMF
     60:1259-1288 (http://cds.iisc.ac.in/faculty/sashi/pub/BenchmarkBubble_IJNF2009.pdf).
     A circular bubble of radius 0.25 centered at (0.5, 0.5) rises through a denser
     liquid (rho = 1000, mu = 10) in a [0,1] x [0,2] box under gravity g = 0.98 over
     t in [0,3]. Case 1 has rho_b = 100, mu_b = 1, sigma = 24.5 (Re = 35, Eo = 10,
     the bubble stays connected); case 2 has rho_b = 1, mu_b = 0.1, sigma = 1.96
     (Re = 35, Eo = 125, thin skirts develop). Benchmark quantities are computed by
     benchmark_quantities.py.
  3: a physical air-water bubble of radius 0.5 mm in a 2 x 4 mm box
     (rho = 1000/1.2, mu = 1e-3/1.8e-5, sigma = 0.072, g = 9.81).

All configurations use free-slip side walls, no-slip top/bottom walls, an ideal-gas
EOS for both fluids with a hydrostatic initial pressure, and a background pressure
that keeps the liquid at low Mach number.

By default the case uses the semi-implicit pressure projection method; pass
--explicit for standard explicit time stepping.

Usage: ./mfc.sh run examples/2D_rising_bubble/case.py -- --case 2 --res 2 --explicit
"""

import argparse
import json

parser = argparse.ArgumentParser(description="2D rising bubble (Hysing et al. 2009 benchmark and air-water)")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT")
parser.add_argument("--res", type=float, default=1.0, help="Resolution multiplier on the 80x160 base grid (default: 1)")
parser.add_argument("--explicit", action="store_true", help="Use explicit time stepping instead of the projection method")
parser.add_argument("--case", type=int, choices=[1, 2, 3], default=1, help="Hysing et al. (2009) test case 1 or 2, or 3 for a physical air-water bubble (default: 1)")
args = parser.parse_args()

Nx = round(80 * args.res)
Ny = round(160 * args.res)

eps = 1e-6
gamma = 1.4
rho_l = 1000.0  # liquid (fluid 1)

if args.case == 3:
    Lx, Ly = 2e-3, 4e-3
    radius = 5e-4
    mu_l = 1e-3
    rho_b = 1.2  # bubble (fluid 2)
    mu_b = 1.8e-5
    sigma = 0.072
    g = 9.81
    p0 = 101325.0
    t_stop, t_save = 0.02, 0.0002
else:
    Lx, Ly = 1.0, 2.0
    radius = 0.25
    mu_l = 10.0
    g = 0.98
    p0 = 1.0e5
    # t_stop, t_save = 3.0, 0.1
    t_stop, t_save = 3.0, 0.03
    if args.case == 1:
        rho_b, mu_b, sigma = 100.0, 1.0, 24.5
    else:
        rho_b, mu_b, sigma = 1.0, 0.1, 1.96

# Hydrostatic initial pressure (liquid column; the bubble's deficit is negligible)
pres = f"{p0} + {rho_l*g}*({Ly} - y)"

case = {
    "run_time_info": "T",
    # Domain
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cfl_adap_dt": "T",
    "ramp_ratio": 2.0,
    "cfl_target": 0.8,
    "n_start": 0,
    "t_save": t_save,
    "t_stop": t_stop,
    "num_patches": 2,
    "model_eqns": "5eq",
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": "rk3",
    "recon_type": "weno",
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "mapped_weno": "T",
    "weno_avg": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    "avg_state": "arithmetic",
    # Free-slip side walls, no-slip top/bottom walls
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Gravity
    "bf_y": "T",
    "k_y": 0.0,
    "w_y": 0.0,
    "p_y": 0.0,
    "g_y": -g,
    # Viscosity
    "viscous": "T",
    "weno_Re_flux": "T",
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%Re(1)": 1.0 / mu_b,
    # Surface tension
    "surface_tension": "T",
    "sigma": sigma,
    "cf_wrt": "T",
    "patch_icpp(1)%cf_val": 0,
    "patch_icpp(2)%cf_val": 1,
    # Liquid background
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": pres,
    "patch_icpp(1)%alpha_rho(1)": (1 - eps) * rho_l,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_b,
    "patch_icpp(1)%alpha(1)": 1 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    # Bubble
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.8,
    "patch_icpp(2)%x_centroid": Lx / 2,
    "patch_icpp(2)%y_centroid": 2 * radius,
    "patch_icpp(2)%radius": radius,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": pres,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(2)%alpha_rho(2)": (1 - eps) * rho_b,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1 - eps,
    # Fluid properties
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
}

if not args.explicit:
    case.update(
        {
            "proj_method": "T",
            "proj_iter_solver": 3,
            "proj_tol_rel": 1.0e-8,
            "proj_max_iters": 500,
            "proj_cfl_ac": 25,
            "proj_normalization": "T",
        }
    )

print(json.dumps(case))
