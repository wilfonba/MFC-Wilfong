#!/usr/bin/env python3

import math
import json

# Numerical setup
d0 = 0.
x0 = -0.1
x1 = 0.1
y0 = 0
y1 = 0.4

Nx = 399
Ny = 799

eps = 1e-6

mydt = 2e-7

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
        't_step_stop'       : 300000,
        't_step_save'       : 5000,
    # =======================================

    # Simulation Algorithm ==================
        'model_eqns'        : 3,
        'alt_soundspeed'    : 'F',
        'adv_alphan'        : 'T',
        'mixture_err'       : 'T',
        'mpp_lim'           : 'F',
        'time_stepper'      : 3,
        'recon_type'        : 1,
        'weno_order'        : 5,  
        'avg_state'         : 2,
        'weno_eps'          : 1e-16,
        'mapped_weno'       : 'F',
        'null_weights'      : 'F',
        'mp_weno'           : 'F',
        'weno_Re_flux'      : 'F',
        'riemann_solver'    : 2,
        'wave_speeds'       : 1,
        'bc_x%beg'          : -15,
        'bc_x%end'          : -15,
        'bc_y%beg'          : -15,
        'bc_y%end'          : -15,
        'num_patches'       : 2,
        'num_fluids'        : 2,
    # =======================================

    # Database Structure Parameters =========
        'format'            : 1,
        'precision'         : 2,
        'cons_vars_wrt'     : 'F',
        'prim_vars_wrt'     : 'T',
        'parallel_io'       : 'T',
    # =======================================

    # Fluid Parameters (Water) ==============
        'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),
        'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00-1.E+00),
        'fluid_pp(1)%Re(1)'         : 100,
    # =======================================

    # Fluid Parameters (Gas) ================
        'fluid_pp(2)%gamma'            : 1.E+00/(1.4E+00-1.E+00),
        'fluid_pp(2)%pi_inf'           : 0.E+00,
        'fluid_pp(2)%Re(1)'         : 10000,
    # =======================================

    # Body Forces ===========================
        #'bf_x'                      : 1,
        #'k_x'                       : 9.80,
        #'w_x'                       : 0,
        #'p_x'                       : 1.57075,
        'bf_y'                      : 1,
        'k_y'                       : 98.,
        'w_y'                       : 0,
        'p_y'                       : 1.57075,
    # ======================================

    # Water Patch ==========================
        'patch_icpp(1)%geometry'    : 3,
        'patch_icpp(1)%x_centroid'  : 0,
        'patch_icpp(1)%y_centroid'  : 1,
        'patch_icpp(1)%length_x'    : 1,
        'patch_icpp(1)%length_y'    : 2,
        'patch_icpp(1)%vel(1)'      : 0.0,
        'patch_icpp(1)%vel(2)'      : 0.0,
        'patch_icpp(1)%vel(3)'      : 0.0,
        'patch_icpp(1)%pres'        : 100000,
        'patch_icpp(1)%alpha_rho(1)': (1-eps)*1000,
        'patch_icpp(1)%alpha_rho(2)': (eps)*100,
        'patch_icpp(1)%alpha(1)'    : 1-eps,
        'patch_icpp(1)%alpha(2)'    : eps,
    # ======================================

    # Bubble Patch ========================
        'patch_icpp(2)%alter_patch(1)' : 'T',                       
        #'patch_icpp(2)%smoothen'    : 'T',
        #'patch_icpp(2)%smooth_patch_id' : 1,
        #'patch_icpp(2)%smooth_coeff'    : 1,
        'patch_icpp(2)%geometry'    : 2,
        'patch_icpp(2)%x_centroid'  : 0.0,
        'patch_icpp(2)%y_centroid'  : 0.1,
        'patch_icpp(2)%z_centroid'  : 0,
        'patch_icpp(2)%radius'      : 0.025,
        'patch_icpp(2)%vel(1)'      : 0.0,
        'patch_icpp(2)%vel(2)'      : 0.0,
        'patch_icpp(2)%vel(3)'      : 0.0,
        'patch_icpp(2)%pres'        : 100000,
        'patch_icpp(2)%alpha_rho(1)': (eps)*1000,
        'patch_icpp(2)%alpha_rho(2)': (1-eps)*100,
        'patch_icpp(2)%alpha(1)'    : eps,
        'patch_icpp(2)%alpha(2)'    : 1-eps,
    # ======================================

}

print(json.dumps(data))
