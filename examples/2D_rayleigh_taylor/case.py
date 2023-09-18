#!/usr/bin/env python3

import math
import json

# Numerical setup
x0 = 0.
x1 = 0.2
y0 = 0.
y1 = 1.2

Nx = 119
Ny = 719

eps = 1e-8

mydt = 1e-6

#Configuration case dictionary
data = {
    # Logistics =============================
        #'case_dir'          : '\'.\'', 
        'run_time_info'     : 'T',
    # =======================================

    # Computational Domain ==================
        'x_domain%beg'      : x0,
        'x_domain%end'      : x1,
        'y_domain%beg'      : y0,
        'y_domain%end'      : y1,
        'm'                 : Nx,
        'n'                 : Ny,
        'p'                 : 0,
        'cyl_coord'        : 'F',
        'dt'                : mydt,
        't_step_start'      : 0,
        't_step_stop'       : 1000000,
        't_step_save'       : 10000,
    # =======================================

    # Simulation Algorithm ==================
        'model_eqns'        : 3,
        'alt_soundspeed'    : 'F',
        'adv_alphan'        : 'T',
        'mixture_err'       : 'T',
        'mpp_lim'           : 'T',
        'time_stepper'      : 3,
        'recon_type'        : 1,
        #'muscl_order'       : 2,
        #'muscl_lim'         : 2,
        #'int_comp'          : 'T',
        'avg_state'         : 2,
        'weno_order'        : 5,  
        'weno_eps'          : 1e-16,
        'mapped_weno'       : 'T',
        'null_weights'      : 'F',
        'mp_weno'           : 'F',
        'weno_Re_flux'      : 'T',
        'riemann_solver'    : 2,
        'wave_speeds'       : 1,
        'bc_x%beg'          : -2,
        'bc_x%end'          : -2,
        'bc_y%beg'          : -15,
        'bc_y%end'          : -15,
        'num_patches'       : 1,
        'num_fluids'        : 2,
    # =======================================

    # Database Structure Parameters =========
        'format'            : 1,
        'precision'         : 2,
        #'prim_vars_wrt'     : 'T',
        #'parallel_io'       : 'T',
        'alpha_wrt(1)'      : 'T',
        'pres_wrt'          : 'T',
    # =======================================

    # Fluid Parameters (Water) ==============
        'fluid_pp(1)%gamma'            : 1.E+00/(1.4E+00-1.E+00),
        'fluid_pp(1)%pi_inf'           : 0.E+00,
        'fluid_pp(1)%Re(1)'            : 1/0.1387, 
    # =======================================

    # Fluid Parameters (Gas) ================
        'fluid_pp(2)%gamma'            : 1.E+00/(1.4E+00-1.E+00),
        'fluid_pp(2)%pi_inf'           : 0.E+00,
        'fluid_pp(2)%Re(1)'            : 1/0.0073,
    # =======================================

    # Body Forces ===========================
        #'bf_x'                      : 1,
        #'k_x'                       : -98.1,
        #'w_x'                       : 0,
        #'p_x'                       : 1.57075,
        'bf_y'                      : 1,
        'k_y'                       : 9.81,
        'w_y'                       : 0,
        'p_y'                       : 1.57075,
    # ======================================

    # Water Patch ==========================
        'patch_icpp(1)%geometry'    : 7,
        'patch_icpp(1)%hcid'        : 201,
        'patch_icpp(1)%x_centroid'  : 0.1,
        'patch_icpp(1)%y_centroid'  : 0.6,
        'patch_icpp(1)%length_x'    : 0.2,
        'patch_icpp(1)%length_y'    : 1.2,
        'patch_icpp(1)%vel(1)'      : 0.0,
        'patch_icpp(1)%vel(2)'      : 0.0,
        'patch_icpp(1)%pres'        : 1e5,
        'patch_icpp(1)%alpha_rho(1)': (1-eps)*1000,
        'patch_icpp(1)%alpha_rho(2)': eps*1,
        'patch_icpp(1)%alpha(1)'    : 1-eps,
        'patch_icpp(1)%alpha(2)'    : eps,
    # ======================================

}

print(json.dumps(data))
