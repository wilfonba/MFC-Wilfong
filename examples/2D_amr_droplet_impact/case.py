#!/usr/bin/env python3
# 2D off-centre collision of two water droplets in air, with surface tension, viscosity, THINC interface
# compression and two levels of block-structured AMR (see docs/documentation/amr.md).
#
# Two droplets of diameter D move toward each other at +-U with their centres offset by B*D across the
# impact axis. They merge, the merged drop spins and stretches into a long ligament, and the ligament
# pinches off near both end bulbs (stretching separation). Dynamic regridding keeps both 2:1 levels
# (finest cell = dx0/4, D/dx = 128) on the interfaces.
#
# Densities, pressure and the air EOS are those of water and air. Surface tension and viscosities are
# raised so the dimensionless numbers (We, Re) are reached at a droplet speed ~1/35 of the air sound speed;
# the liquid stiffness pi_inf is lowered so the water sound speed matches air's rather than setting a
# 4x smaller time step. The flow stays nearly incompressible (Mach < 0.2, the peak being the air squeezed
# out of the gap before contact).
import json
import math

# Physical setup (SI)
D = 1.0e-3  # droplet diameter
R = D / 2
U = 10.0  # speed of each droplet (closing speed 2U)
B = 0.5  # impact parameter: centre offset across the impact axis / D
rho_l, rho_g = 1000.0, 1.2
p0 = 1.0e5
We = 60.0  # rho_l (2U)^2 D / sigma
Re = 600.0  # rho_l (2U) D / mu_l
sigma = rho_l * (2 * U) ** 2 * D / We
mu_l = rho_l * (2 * U) * D / Re
mu_g = mu_l / 55.0  # water/air viscosity ratio
gam_l, pi_l = 6.12, 2.0e7  # water stiffened gas (pi_inf softened from 3.43e8)
gam_g = 1.4

# Grid: level 0 is N x N over a 12D square (D/dx0 = 32); two 2:1 levels give D/dx = 128 on the finest level.
# The wide domain lets one level-1 box (at most N/2 cells wide) hold the whole stretched ligament.
N = 384
L = 12 * D
dx_f = L / N / 4

# Lock-step AMR: every level advances at dt, set from the finest cell (acoustic limit, CFL 0.4).
c_max = max(math.sqrt(gam_l * (p0 + pi_l) / rho_l), math.sqrt(gam_g * p0 / rho_g)) + 2 * U
dt = 0.4 * dx_f / c_max
t_end = 4.5e-4  # through the ligament pinch-off
Nt = int(round(t_end / dt / 300)) * 300

xc = 1.25 * R  # droplet centres at (-+xc, -+B*R): a 0.5R gap along the impact axis
yc = B * R
eps = 1e-9

data = {
    # Logistics
    "run_time_info": "T",
    # Computational domain
    "x_domain%beg": -L / 2,
    "x_domain%end": L / 2,
    "y_domain%beg": -L / 2,
    "y_domain%end": L / 2,
    "m": N - 1,
    "n": N - 1,
    "p": 0,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": Nt // 30,
    # Simulation algorithm
    "num_patches": 3,
    "model_eqns": 2,
    "num_fluids": 2,
    "mpp_lim": "T",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "F",
    "weno_avg": "T",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "int_comp": 1,  # THINC; MTHINC (2) is stable here only while each AMR level is a single block
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma,
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # Patch 1: air at rest
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 2 * L,
    "patch_icpp(1)%length_y": 2 * L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1 - eps,
}

# Patches 2 and 3: the droplets, at the Laplace pressure (2D: sigma/R), moving toward each other
for i, s in [(2, -1), (3, 1)]:
    data.update(
        {
            f"patch_icpp({i})%geometry": 2,
            f"patch_icpp({i})%x_centroid": s * xc,
            f"patch_icpp({i})%y_centroid": s * yc,
            f"patch_icpp({i})%radius": R,
            f"patch_icpp({i})%alter_patch(1)": "T",
            f"patch_icpp({i})%smoothen": "T",
            f"patch_icpp({i})%smooth_patch_id": 1,
            f"patch_icpp({i})%smooth_coeff": 2.0,
            f"patch_icpp({i})%vel(1)": -s * U,
            f"patch_icpp({i})%vel(2)": 0.0,
            f"patch_icpp({i})%pres": p0 + sigma / R,
            f"patch_icpp({i})%alpha_rho(1)": (1 - eps) * rho_l,
            f"patch_icpp({i})%alpha_rho(2)": eps * rho_g,
            f"patch_icpp({i})%alpha(1)": 1 - eps,
            f"patch_icpp({i})%alpha(2)": eps,
        }
    )

data.update(
    {
        # Fluid 1: water (stiffened gas), fluid 2: air (ideal gas). Re(1) = 1/mu.
        "fluid_pp(1)%gamma": 1.0 / (gam_l - 1.0),
        "fluid_pp(1)%pi_inf": gam_l * pi_l / (gam_l - 1.0),
        "fluid_pp(1)%Re(1)": 1.0 / mu_l,
        "fluid_pp(2)%gamma": 1.0 / (gam_g - 1.0),
        "fluid_pp(2)%pi_inf": 0.0,
        "fluid_pp(2)%Re(1)": 1.0 / mu_g,
        # --- AMR ---
        # Two 2:1 levels; regrid every 8 steps (an interface moves < 1 finest cell in that time).
        # A coarse/fine face must never cut an interface, and a level-2 block cannot straddle two level-1
        # blocks, so the level-1 box must not be split across the liquid: a wide pad (amr_buf, level-0 cells,
        # also the level-2 pad inside its parent), a low clustering efficiency so the liquid gets one compact
        # level-1 box, and a box cap (N/2) large enough for the whole ligament.
        "amr": "T",
        "amr_max_level": 2,
        "amr_ref_ratio": 2,
        "amr_block_beg(1)": N // 2 - 46,
        "amr_block_end(1)": N // 2 + 45,
        "amr_block_beg(2)": N // 2 - 34,
        "amr_block_end(2)": N // 2 + 33,
        "amr_regrid_int": 8,
        "amr_tag_eps": 0.1,
        "amr_buf": 8,
        "amr_cluster_eff": 0.1,
        "amr_max_grid_size": N // 2,
        "amr_max_blocks": 512,
        "amr_prolong": 1,
    }
)

print(json.dumps(data))
