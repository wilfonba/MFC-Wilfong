#!/usr/bin/env python3
import math
import json
import argparse

parser = argparse.ArgumentParser(prog="3D IGR Performance Test", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    "--mfc",
    type=json.loads,
    default="{}",
    metavar="DICT",
    help="MFC's toolchain's internal state.",
)
parser.add_argument("-dim", "--dim", metavar="dim", type=int, help="Dimension", default=1)
parser.add_argument("-nt", "--nt", metavar="nt", type=int, help="number of time-steps", default=100)
parser.add_argument("-ns", "--ns", metavar="ns", type=int, help="number of saves", default=1)

args = parser.parse_args()

pA = 1.0
rhoA = 1
gam = 1.4
c_l = math.sqrt(1.4 * pA / rhoA)

M = 10
pS = 10*pA
velS = M*c_l
rhoS = 1.5*rhoA

leng = 1

Ny = 399
Nx = Ny
Nz = Ny

dx = 6*leng/Nx

cfl = 1.0
dt = cfl * dx / (M*c_l + velS)

eps = 1e-5

# Configuring case dictionary
case = {
        "rdma_mpi": "T",
        # Logistics
        # "run_time_info": "T",
        # Computational Domain Parameters
        "x_domain%beg": 0,
        "x_domain%end": 8*leng,
        "y_domain%beg": 0*leng,
        "y_domain%end": 8*leng,
        "z_domain%beg": 0*leng,
        "z_domain%end": 8*leng,
        "m": Nx,
        "n": Ny,
        "p": Nz,

        "dt": dt,
        "t_step_start": 0,
        "t_step_stop": args.nt,
        "t_step_save": int(args.nt/args.ns),

        # Simulation Algorithm Parameters
        "num_patches": 2,
        "model_eqns": 2,
        "alt_soundspeed": "F",
        "num_fluids": 2,
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
        "alf_factor": 15,
        "num_igr_iters": 2,

        "bc_x%beg": -3,
        "bc_x%end": -3,
        "bc_y%beg": -3,
        "bc_y%end": -3,
        "bc_z%beg": -3,
        "bc_z%end": -3,

        "num_bc_patches": 1,

        # Formatted Database Files Structure Parameters
        "format": 1,
        "prim_vars_wrt": "T",
        "parallel_io": "T",

        # Patch 1: Background
        "patch_icpp(1)%geometry": 9,
        "patch_icpp(1)%x_centroid": 0.0,
        "patch_icpp(1)%y_centroid": 0.0,
        "patch_icpp(1)%z_centroid": 0.0,
        "patch_icpp(1)%length_x": 100*leng,
        "patch_icpp(1)%length_y": 100*leng,
        "patch_icpp(1)%length_z": 100*leng,
        "patch_icpp(1)%vel(1)": 0.0*velS,
        "patch_icpp(1)%vel(2)": 0.0e00,
        "patch_icpp(1)%vel(3)": 0.0e00,
        "patch_icpp(1)%pres": pA,
        "patch_icpp(1)%alpha_rho(1)": eps*rhoA,
        "patch_icpp(1)%alpha(1)": eps,
        "patch_icpp(1)%alpha_rho(2)": (1-eps)*rhoA,
        "patch_icpp(1)%alpha(2)": 1-eps,

        # Patch 2: Shocked state
        "patch_icpp(2)%geometry": 8,
        "patch_icpp(2)%alter_patch(1)": "T",
        "patch_icpp(2)%smoothen": "T",
        "patch_icpp(2)%smooth_patch_id": 1,
        "patch_icpp(2)%smooth_coeff": 0.1,
        "patch_icpp(2)%pres": pS,
        "patch_icpp(2)%alpha_rho(1)": (1-eps)*rhoS,
        "patch_icpp(2)%alpha(1)": 1-eps,
        "patch_icpp(2)%alpha_rho(2)": eps*rhoS,
        "patch_icpp(2)%alpha(2)": eps,

        # Fluids Physical Parameters
        "fluid_pp(1)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
        "fluid_pp(1)%pi_inf": 0.0,
        "fluid_pp(2)%gamma": 1.0e00 / (1.4e00 - 1.0e00),
        "fluid_pp(2)%pi_inf": 0.0,
        "viscous": "T",
        "fluid_pp(1)%Re(1)": 1.81e5,
        "fluid_pp(2)%Re(1)": 1.81e5,
    }

if args.dim == 1:
    case.update(
        {
        "patch_bc(1)%dir": 1,
        "patch_bc(1)%loc": -1,
        "patch_bc(1)%geometry": 2,
        "patch_bc(1)%type": -17,
        "patch_bc(1)%centroid(2)": 4*leng,
        "patch_bc(1)%centroid(3)": 4*leng,
        "patch_bc(1)%radius": 1.5*leng/2,

        "patch_icpp(2)%x_centroid": 0.0,
        "patch_icpp(2)%y_centroid": 4*leng,
        "patch_icpp(2)%z_centroid": 4*leng,
        "patch_icpp(2)%radius": leng/2,
        "patch_icpp(2)%vel(1)": velS,
        "patch_icpp(2)%vel(2)": 0.0e00,
        "patch_icpp(2)%vel(3)": 0.0e00,
        }
    )
elif args.dim == 2:
    case.update(
        {
        "patch_bc(1)%dir": 2,
        "patch_bc(1)%loc": -1,
        "patch_bc(1)%geometry": 2,
        "patch_bc(1)%type": -17,
        "patch_bc(1)%centroid(1)": 4*leng,
        "patch_bc(1)%centroid(3)": 4*leng,
        "patch_bc(1)%radius": 1.5*leng/2,

        "patch_icpp(2)%x_centroid": 4*leng,
        "patch_icpp(2)%y_centroid": 0,
        "patch_icpp(2)%z_centroid": 4*leng,
        "patch_icpp(2)%radius": leng/2,
        "patch_icpp(2)%vel(1)": 0.0e00,
        "patch_icpp(2)%vel(2)": velS,
        "patch_icpp(2)%vel(3)": 0.0e00,
        }
    )
elif args.dim == 3:
    case.update(
        {
        "patch_bc(1)%dir": 2,
        "patch_bc(1)%loc": -1,
        "patch_bc(1)%geometry": 2,
        "patch_bc(1)%type": -17,
        "patch_bc(1)%centroid(1)": 4*leng,
        "patch_bc(1)%centroid(2)": 4*leng,
        "patch_bc(1)%radius": 1.5*leng/2,

        "patch_icpp(2)%x_centroid": 4*leng,
        "patch_icpp(2)%y_centroid": 4*leng,
        "patch_icpp(2)%z_centroid": 0.0,
        "patch_icpp(2)%radius": leng/2,
        "patch_icpp(2)%vel(1)": 0.0e00,
        "patch_icpp(2)%vel(2)": 0.0e00,
        "patch_icpp(2)%vel(3)": velS,
        }
    )

if 'single' in args.mfc:
    if args.mfc['single']:
        case.update({"precision": 1})
    else:
        case.update({"precision": 2})
else:
    case.update({"precision": 2})


if __name__ == "__main__":
    print(json.dumps(case,indent=4))
