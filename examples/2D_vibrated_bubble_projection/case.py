#!/usr/bin/env python3
"""
2D bubbles beneath a Faraday-unstable free surface in a vertically vibrated water column, with the all-Mach pressure projection:
a proof of concept for 3D simulations of vibrated bubbly liquids.

A water column of depth H in a closed container, under a partly evacuated gas headspace at p_amb, is vibrated vertically: in
its frame the effective gravity is g + Gamma g sin(w t), through MFC's body force accel_y = g_y + k_y sin(w_y t - p_y). The
hydrostatic pressure at depth h then oscillates by rho Gamma g h, and a gas bubble's volume oscillates against it, smallest
when the effective gravity is largest. For a 3D bubble in a large tank, the time-averaged buoyancy rho <V (g + a)> therefore
falls with depth and reverses below the critical depth (Bleich 1956; Jameson 1966)

    h_c = 2 gamma p_amb / (rho g Gamma^2)      (quasi-static, linear in the volume oscillation)

~29 m / Gamma^2 for water and air at 1 atm, so the headspace is evacuated to 0.1 bar (as in experiments), which brings h_c to
~4.5 cm at Gamma = 8. Two bubbles start at rest side by side, placed by that 3D estimate at h_c/2 and 2 h_c.

In 2D this estimate does not hold. A 2D bubble is a cylinder, which holds far more gas per unit width of liquid than a
sphere, so the water above it bounces on its gas at omega_0^2 ~ gamma p w/(rho h A), with w the bubble's share of the
periodic width, h the water above it and A its area: ~63 t_c^-1 for the deep bubble, next to the forcing's 69.5, against far
above it in 3D. The deep bubble is driven near that resonance and the column's pressure lags and overshoots the hydrostatic
(the same column without bubbles follows it to 0.2%).

The forcing must also leave the bubble intact. Its velocity relative to the water, of amplitude Gamma g/w, sets a Weber
number rho (Gamma g/w)^2 D/sigma that must stay of order one (~0.9 at the defaults, 100 Hz). Whenever the effective gravity
points up, water lies above the bubble's lower surface, which is Rayleigh-Taylor unstable at wavelengths above
2 pi sqrt(sigma/(rho Gamma g)), ~6 mm here, longer than the bubble. The column is short enough that the pressure at its
bottom stays positive through the forcing's trough.

The free surface is far past its Faraday threshold (roughly Gamma ~ 0.3 for water at 100 Hz, from the viscous damping of
its ~6 mm, 50 Hz capillary-gravity waves): its subharmonic waves, and how they interact with the bubbles beneath, are part
of the problem. The headspace (--headspace, a fraction of H) is large enough that the bubbles' breathing hardly moves its
pressure, and closing it gives both solvers the same boundaries: periodic sides, a slip wall at the top, no-slip below.

Nondimensionalization (code units), as in examples/2D_interface_breakup: rho_w = H = g0 = 1, t_c = sqrt(H/g0),
U_c = sqrt(g0 H), p_c = rho_w g0 H. The control number is the bubble's Bond number Bo = rho_w g D^2/sigma, which sets the
length scale (D = sqrt(Bo sigma/(rho_w g)), H = D/(D/H)) for the real water-air pair; Re = rho_w U_c H/mu_i follows.
Water keeps its physical stiffness, and the projection steps at the flow's pace, not water's sound speed.
--explicit runs the explicit HLLC solver instead (low_Mach = 2, conservative surface tension), with water softened to the
gas's sound speed, which it must still resolve.

Limitation: the five-equation model compresses mixed interface cells with the advected mixture's stiff modulus rather than
Wood's, so a bubble's volume oscillation is under-predicted, and with it the Bjerknes force: the simulated h_c lies deeper
than the formula's by the inverse of the captured fraction, which depends on resolution and on the size of the oscillation.
"""

import argparse
import json
import math
import os
import sys

parser = argparse.ArgumentParser(description="2D bubbles under a Faraday-unstable free surface in a vibrated water column, all-Mach pressure projection")
parser.add_argument("--Bo", type=float, default=2.2, help="bubble Bond number rho_w g D^2 / sigma; sets the length scale (default: %(default)s, D ~ 4 mm)")
parser.add_argument("--forcing-ratio", type=float, default=8.0, help="Gamma = k_y / g0: vibration acceleration amplitude / gravity (default: %(default)s)")
parser.add_argument("--forcing-freq", type=float, default=69.5, help="w_y * t_c: dimensionless angular forcing frequency (default: %(default)s, ~100 Hz)")
parser.add_argument("--p-amb", type=float, default=0.1, help="headspace pressure [bar] (default: %(default)s)")
parser.add_argument("--headspace", type=float, default=0.25, help="gas headspace height above the free surface, in H (default: %(default)s)")
parser.add_argument("--cells-per-D", type=int, default=100, help="cells across a bubble diameter (default: %(default)s)")
parser.add_argument("--periods", type=float, default=10.0, help="forcing periods to simulate (default: %(default)s)")
parser.add_argument("--frames-per-period", type=float, default=20.0, help="outputs per forcing period (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit HLLC at the acoustic limit (low_Mach = 2) instead of the pressure projection")
args, _ = parser.parse_known_args()

# 1. Physical fluid properties (SI): fluid 1 water (stiffened gas), fluid 2 air (ideal gas)
rho_w_SI = 1000.0
mu_w_SI, mu_g_SI = 1.0e-3, 1.8e-5
sigma_SI, g_SI, p_amb_SI = 0.072, 9.81, 1.0e5 * args.p_amb
rho_g_SI = p_amb_SI / (287.0 * 293.0)  # air at the headspace pressure, 20 C
gam_w, pi_inf_w_SI, gam_g = 6.12, 3.43e8, 1.4

# 2. Length scale from the Bond number; the column depth is a fixed multiple of the bubble diameter
D_over_H = 1.0 / 30.0
D_SI = math.sqrt(args.Bo * sigma_SI / (rho_w_SI * g_SI))
H_SI = D_SI / D_over_H
U_SI, t_SI = math.sqrt(g_SI * H_SI), math.sqrt(H_SI / g_SI)
p_SI = rho_w_SI * g_SI * H_SI

# 3. Code units
rho_w, H, g0 = 1.0, 1.0, 1.0
rho_g = rho_g_SI / rho_w_SI
D, R = D_over_H * H, 0.5 * D_over_H * H
sigma = rho_w * g0 * D**2 / args.Bo
p_amb = p_amb_SI / p_SI
Re_w = rho_w_SI * U_SI * H_SI / mu_w_SI
Re_g = rho_w_SI * U_SI * H_SI / mu_g_SI
Gamma = args.forcing_ratio
w_y = args.forcing_freq / math.sqrt(H / g0)
T_osc = 2.0 * math.pi / w_y
h_c = 2.0 * gam_g * p_amb / (rho_w * g0 * Gamma**2)  # critical depth, code units

# 4. Domain: the column of depth H, 0.4 H (12 D) wide and periodic, so the two bubbles, 0.2 H apart either way, are evenly
# spaced, under the headspace; the free surface is at y = H
Lx, Ly = 0.4 * H, (1.0 + args.headspace) * H
dx = D / args.cells_per_D
Nx, Ny = int(round(Lx / dx)), int(round(Ly / dx))
# Bubbles: at h_c/2 and 2 h_c, by the 3D estimate
bubbles = [(-0.25 * Lx, H - 0.5 * h_c), (0.25 * Lx, H - 2.0 * h_c)]

# Hydrostatic under the background gravity from p_amb at the free surface; a bubble adds its Laplace jump (2D: sigma/R)
p_water = f"({p_amb!r} + {rho_w * g0!r}*({H!r} - y))"
p_bubble = f"({p_amb + sigma / R!r} + {rho_w * g0!r}*({H!r} - y))"
p_head = f"({p_amb!r} + {rho_g * g0!r}*({H!r} - y))"

# Stiffened-gas coefficients (MFC): Gamma = 1/(gam - 1), Pi = gam pi_inf/(gam - 1)
pi_inf_w = pi_inf_w_SI / p_SI
if args.explicit:
    # Soften water so its sound speed at p_amb matches the gas's (2D_interface_breakup's "match_gas"): the explicit step
    # grows by c_water/c_gas ~ 4.4 while water stays at Mach ~1e-2, and the column's quarter-wave resonance, c/(4H), stays
    # ~7x the forcing frequency
    pi_inf_w = gam_g * p_amb / rho_g * rho_w / gam_w - p_amb

print(f"  D = {D_SI * 1e3:.3g} mm, H = {H_SI * 100:.3g} cm, t_c = {t_SI * 1e3:.3g} ms, p_c = {p_SI:.4g} Pa", file=sys.stderr)
print(f"  Bo = {args.Bo:.3g}, Re_water = {Re_w:.3g}, Re_gas = {Re_g:.3g}, sigma = {sigma:.3g}, p_amb = {p_amb:.4g}", file=sys.stderr)
print(f"  forcing: Gamma = {Gamma:.3g}, f = {w_y / t_SI / (2 * math.pi):.3g} Hz, amplitude = {Gamma * g_SI / (w_y / t_SI) ** 2 * 1e3:.3g} mm", file=sys.stderr)
print(f"  relative-oscillation Weber number {rho_w_SI * (Gamma * g_SI / (w_y / t_SI)) ** 2 * D_SI / sigma_SI:.3g}", file=sys.stderr)
print(f"  Rayleigh-Taylor cutoff wavelength {2 * math.pi * math.sqrt(sigma_SI / (rho_w_SI * Gamma * g_SI)) * 1e3:.3g} mm (D = {D_SI * 1e3:.3g} mm)", file=sys.stderr)
print(f"  lowest pressure, at the bottom in the trough: {(p_amb + rho_w * g0 * H * (1 - Gamma)) * p_SI:.4g} Pa", file=sys.stderr)
print(f"  critical depth h_c = {h_c * H_SI * 100:.3g} cm ({h_c:.3g} H); bubble depths {', '.join(f'{(H - y) * H_SI * 100:.3g} cm' for _, y in bubbles)}", file=sys.stderr)
print(f"  sound speeds: water {math.sqrt(gam_w * (p_amb + pi_inf_w) / rho_w) * U_SI:.4g} m/s, gas {math.sqrt(gam_g * p_amb / rho_g) * U_SI:.4g} m/s", file=sys.stderr)
print(f"  grid {Nx} x {Ny}, {args.cells_per_D} cells per D; {args.periods:.3g} periods = {args.periods * T_osc * t_SI:.3g} s", file=sys.stderr)
if 2.0 * h_c > H - D or p_amb + rho_w * g0 * H * (1 - Gamma) <= 0:
    print("  warning: the deep bubble does not fit above the bottom, or the bottom goes into tension; adjust --forcing-ratio or --p-amb", file=sys.stderr)

eps = 1.0e-8
water = {"alpha_rho(1)": (1 - eps) * rho_w, "alpha_rho(2)": eps * rho_g, "alpha(1)": 1 - eps, "alpha(2)": eps, "cf_val": 0}
gas = {"alpha_rho(1)": eps * rho_w, "alpha_rho(2)": (1 - eps) * rho_g, "alpha(1)": eps, "alpha(2)": 1 - eps, "cf_val": 1}

case = {
    "run_time_info": "T",
    "x_domain%beg": -0.5 * Lx,
    "x_domain%end": 0.5 * Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cfl_adap_dt": "T",
    "cfl_target": 0.4,
    "n_start": 0,
    "t_stop": args.periods * T_osc,
    "t_save": T_osc / args.frames_per_period,
    "model_eqns": 2,
    "num_fluids": 2,
    "num_patches": 2 + len(bubbles),
    "time_stepper": 3,
    "weno_order": 3,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Periodic sides; a slip wall at the headspace's top, no-slip below
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -16,
    "bc_y%end": -15,
    "int_comp": 2,
    "ic_beta": 0.6,
    "viscous": "T",
    "weno_Re_flux": "F",
    "fluid_pp(1)%Re(1)": Re_w,
    "fluid_pp(2)%Re(1)": Re_g,
    "surface_tension": "T",
    "sigma": sigma,
    # Effective gravity in the vibrated frame: g_y + k_y sin(w_y t)
    "bf_y": "T",
    "g_y": -g0,
    "k_y": Gamma * g0,
    "w_y": w_y,
    "p_y": 0.0,
    "fluid_pp(1)%gamma": 1.0 / (gam_w - 1.0),
    "fluid_pp(1)%pi_inf": gam_w * pi_inf_w / (gam_w - 1.0),
    "fluid_pp(2)%gamma": 1.0 / (gam_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
}
if args.explicit:
    case.update({"low_Mach": 2, "surface_tension_model": "conservative"})
else:
    # The well-balanced surface tension model needs the projection's face pressure gradient
    case.update({"proj_method": "T", "proj_max_acfl": 200, "diff_sts": "T", "diff_sts_max": 64, "surface_tension_model": "well_balanced"})
patches = [
    ({"geometry": 3, "x_centroid": 0.0, "y_centroid": 0.5 * Ly, "length_x": Lx, "length_y": Ly, "pres": p_water}, water),
    ({"geometry": 3, "x_centroid": 0.0, "y_centroid": H + 0.5 * (Ly - H), "length_x": Lx, "length_y": Ly - H, "pres": p_head, "alter_patch(1)": "T"}, gas),
]
for xb, yb in bubbles:
    geo = {"geometry": 2, "x_centroid": xb, "y_centroid": yb, "radius": R, "pres": p_bubble, "alter_patch(1)": "T", "smoothen": "T", "smooth_patch_id": 1, "smooth_coeff": 0.95}
    patches.append((geo, gas))
for k, (geo, state) in enumerate(patches, 1):
    case.update({f"patch_icpp({k})%{key}": val for key, val in {**geo, **state, "vel(1)": 0.0, "vel(2)": 0.0}.items()})

# For analyze.py: the free surface's height, p_amb and h_c, in code units
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vibrated_bubble.json"), "w") as f:
    json.dump({"H": H, "p_amb": p_amb, "h_c": h_c}, f)

print(json.dumps(case))
