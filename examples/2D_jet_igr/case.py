#!/usr/bin/env python3
import math
import json

pA = 1
rhoA = 1
gam = 1.4
c_l = math.sqrt(1.4 * pA / rhoA)

M = 10
pS = 5*pA
velS = M*c_l
rhoS = 2*rhoA

leng = 1
Ny = 499
Nx = Ny * 1
dx = leng / Nx

time_end = 5 * leng / (M * c_l)
cfl = 1.0

dt = cfl * dx / (M + 1)*c_l
Nt = int(time_end / dt)

eps = 1e-5

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -2 * leng,
            "x_domain%end": 2 * leng,
            "y_domain%beg": -2 * leng,
            "y_domain%end": 2 * leng,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            # "t_step_stop": Nt,
            # "t_step_save": int(Nt / 100.0),
            "t_step_stop": 20,
            "t_step_save": 20,#int(Nt / 20.0),
            # Simulation Algorithm Parameters
            "num_patches": 2,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "F",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "igr": "T",
            "elliptic_smoothing": "T",
            "elliptic_smoothing_iters": 100,

            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,

            # "num_bc_patches": 1,
            # Top jet inlet at left side of domain
            # "patch_bc(1)%dir": 1,
            # "patch_bc(1)%loc": -1,
            # "patch_bc(1)%geometry": 1,
            # "patch_bc(1)%type": -17,
            # "patch_bc(1)%centroid(2)": 0,
            # "patch_bc(1)%length(2)": 6*leng/8,

            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "T",
            # "file_per_process": "T",

            # Patch 1: Background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 8 * leng,
            "patch_icpp(1)%length_y": 4 * leng,
            "patch_icpp(1)%vel(1)": 0.0e00,
            "patch_icpp(1)%vel(2)": 0.0e00,
            "patch_icpp(1)%pres": pA,
            "patch_icpp(1)%alpha_rho(1)": rhoA,
            "patch_icpp(1)%alpha(1)": 1,

            # Patch 2: Shocked state
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -2 * leng,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%length_x": leng / 4,
            "patch_icpp(2)%length_y": leng / 2,
            "patch_icpp(2)%vel(1)": velS,
            "patch_icpp(2)%vel(2)": 0.0e00,
            "patch_icpp(2)%pres": pS,
            "patch_icpp(2)%alpha_rho(1)": rhoS,
            "patch_icpp(2)%alpha(1)": 1,

            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,

        }
    )
)
