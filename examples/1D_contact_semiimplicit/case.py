#!/usr/bin/env python3
"""
1D air/water contact discontinuity with the semi-implicit projection method.

An exact-solution test for the projection method at a high density ratio with
NO body force and NO surface tension, so that the density ratio is the only
stressor.  A water slab sits in air at uniform pressure and uniform velocity;
both phases are inviscid and nothing drives the flow, so the exact solution is
that the initial profile translates at the imposed velocity forever, unchanged.

  --velocity 0    stationary contact: an exact steady state, nothing may move
  --velocity U    the profile advects at U and must return to itself after
                  L/U of physical time (the domain is periodic)

Any drift in the phase densities, the pressure or the velocity is therefore
error, and its size is the metric -- there is no judgement call about what the
answer should be.  --acfl sets the fixed step as a multiple of the acoustic
limit in water, which is the quantity the projection method exists to exceed.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D air/water contact discontinuity, semi-implicit projection method")
parser.add_argument(
    "--velocity",
    type=float,
    default=0.0,
    help="uniform advection velocity [m/s]; 0 => stationary contact (default: %(default)s)",
)
parser.add_argument(
    "--acfl",
    type=float,
    default=30.0,
    help="fixed dt as a multiple of the water acoustic limit (default: %(default)s)",
)
parser.add_argument(
    "--explicit",
    action="store_true",
    help="run the explicit solver instead, as a control",
)
args, _ = parser.parse_known_args()

# stiffened-gas water and ideal air, ~833:1 density ratio
gamma_w, p_inf_w, rho_w = 4.4, 6.0e8, 1000.0
gamma_a, p_inf_a, rho_a = 1.4, 0.0, 1.2
p0 = 1.0e5

N = 200
L = 1.0
dx = L / N

# the slab occupies the middle half of a periodic domain, so there are two
# contacts and neither sits on a boundary
slab_beg, slab_end = 0.25 * L, 0.75 * L

c_w = math.sqrt(gamma_w * (p0 + p_inf_w) / rho_w)
dt_acoustic = dx / c_w
dt = args.acfl * dt_acoustic

# long enough for one full traversal at the reference velocity of 5 m/s, so a
# sweep in --velocity compares equal amounts of interface transport
T_end = L / 5.0
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
    "num_patches": 3,
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
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "parallel_io": "F",
    # water slab
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.5 * (slab_beg + slab_end),
    "patch_icpp(1)%length_x": slab_end - slab_beg,
    "patch_icpp(1)%vel(1)": args.velocity,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho_w,
    "patch_icpp(1)%alpha_rho(2)": 0.0,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    # the air surrounding it, as the two remaining segments
    "patch_icpp(2)%geometry": 1,
    "patch_icpp(2)%x_centroid": 0.5 * slab_beg,
    "patch_icpp(2)%length_x": slab_beg,
    "patch_icpp(2)%vel(1)": args.velocity,
    "patch_icpp(2)%pres": p0,
    "patch_icpp(2)%alpha_rho(1)": 0.0,
    "patch_icpp(2)%alpha_rho(2)": rho_a,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(3)%geometry": 1,
    "patch_icpp(3)%x_centroid": 0.5 * (slab_end + L),
    "patch_icpp(3)%length_x": L - slab_end,
    "patch_icpp(3)%vel(1)": args.velocity,
    "patch_icpp(3)%pres": p0,
    "patch_icpp(3)%alpha_rho(1)": 0.0,
    "patch_icpp(3)%alpha_rho(2)": rho_a,
    "patch_icpp(3)%alpha(1)": 0.0,
    "patch_icpp(3)%alpha(2)": 1.0,
    "fluid_pp(1)%gamma": 1.0 / (gamma_w - 1.0),
    "fluid_pp(1)%pi_inf": gamma_w * p_inf_w / (gamma_w - 1.0),
    "fluid_pp(2)%gamma": 1.0 / (gamma_a - 1.0),
    "fluid_pp(2)%pi_inf": gamma_a * p_inf_a / (gamma_a - 1.0),
}

if not args.explicit:
    case.update(
        {
            "proj_method": "T",
            "proj_iter_solver": 1,
            "proj_tol": 1.0e-9,
            "proj_max_iters": 200,
        }
    )

print(json.dumps(case))
