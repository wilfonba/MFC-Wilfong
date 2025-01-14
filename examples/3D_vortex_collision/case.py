#!/usr/bin/env python3
import json

eps = 1e-9

r1 = 1
r2 = 0.25
v1 = 1
vr = 1
off1 = -1.1*r2
off2 = 1.1*r2

xb = -6*r2
xe = 6*r2
yb = -4*r1
ye = 4*r1
zb = -4*r1
ze = 4*r1

Ny = 699
Nz = 699
Nx = int((12*r2)/(8*r1)*Ny)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": xb,
            "x_domain%end": xe,
            "y_domain%beg": yb,
            "y_domain%end": ye,
            "z_domain%beg": zb,
            "z_domain%end": ze,
            "m": Nx,
            "n": Ny,
            "p": Nz,
            "cfl_adap_dt": "T",
            "cfl_target": 0.5,
            "n_start": 0,
            "t_stop": 10,
            "t_save": 2,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "num_fluids": 3,
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_z%beg": -3,
            "bc_z%end": -3,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            # "prim_vars_wrt": "T",
            "alpha_wrt(1)": "T",
            # "alpha_wrt(2)": "T",
            # "alpha_wrt(3)": "T",
            "vel_wrt(1)": "T",
            "vel_wrt(2)": "T",
            "vel_wrt(3)": "T",
            "parallel_io": "T",
            # Patch 1: Base
            "patch_icpp(1)%geometry": 13,
            "patch_icpp(1)%hcid": 302,
            "patch_icpp(1)%x_centroid": (xe + xb)/2,
            "patch_icpp(1)%y_centroid": (ye + yb)/2,
            "patch_icpp(1)%z_centroid": (ze + zb)/2,
            "patch_icpp(1)%length_x": (xe - xb),
            "patch_icpp(1)%length_y": (ye - yb),
            "patch_icpp(1)%length_z": (ze - zb),
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(1)%pres": 1,
            "patch_icpp(1)%alpha_rho(1)": 1 - eps,
            "patch_icpp(1)%alpha(1)": 1 - eps,
            "patch_icpp(1)%alpha_rho(2)": eps,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(1)%alpha_rho(3)": eps,
            "patch_icpp(1)%alpha(3)": eps,
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(1)%pi_inf": 0,
            "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(2)%pi_inf": 0,
            "fluid_pp(3)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
            "fluid_pp(3)%pi_inf": 0,
        }
    )
)
