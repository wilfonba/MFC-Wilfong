#!/usr/bin/env python3
"""
1D air/water contact discontinuity on the staggered (MAC) grid.

The same exact-solution gate the projection method uses, run with the staggered
discretization instead: scalars at cell centres, velocity on the faces.  A water
slab sits in air at uniform pressure and uniform velocity, so the exact solution
is that the initial profile translates unchanged.  Any drift in the phase
densities, the pressure or the velocity is therefore error.

Running with `run_time_info` on -- which this case sets -- prints three
acceptance tests before the first step:

  * the operator gate, which checks that the discrete divergence and gradient
    are exact negative adjoints and that composing them gives the compact
    Laplacian rather than the collocated grid's checkerboard-blind wide one;
  * scalar transport, which reports free-stream preservation, conservation and
    the accuracy of one full period of advection;
  * the interface test, which checks Abgrall's condition -- that a state at
    uniform pressure and velocity stays that way across an 833:1 jump.

All three pass.  The time integration does NOT yet: the scheme carries no
acoustic dissipation and diverges around step 32 at an acoustic CFL of 0.5.  See
docs/staggered_handoff.md for where that stands.  Use --acfl well below one when
looking at the time stepping at all.
"""

import argparse
import json
import math

parser = argparse.ArgumentParser(description="1D air/water contact discontinuity, staggered (MAC) grid")
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
parser.add_argument(
    "--accel",
    type=float,
    default=0.0,
    help="background body-force acceleration [m/s^2]; 0 => no body force (default: %(default)s)",
)
parser.add_argument(
    "--osc-ratio",
    type=float,
    default=0.0,
    help="oscillatory acceleration amplitude as a multiple of --accel (default: %(default)s)",
)
parser.add_argument(
    "--osc-freq",
    type=float,
    default=0.0,
    help="angular frequency of the oscillatory acceleration [rad/s] (default: %(default)s)",
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

if args.accel != 0.0:
    # A body force uniform in space accelerates both phases equally, so in a
    # periodic domain it drives no relative motion and the exact solution is
    # still an unchanged density profile at uniform pressure
    case.update(
        {
            "bf_x": "T",
            "g_x": args.accel,
            "k_x": args.accel * args.osc_ratio,
            "w_x": args.osc_freq,
            "p_x": 0.0,
        }
    )

if not args.explicit:
    case.update({"stagger": "T"})

print(json.dumps(case))
