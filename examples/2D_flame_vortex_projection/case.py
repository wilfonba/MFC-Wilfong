#!/usr/bin/env python3
"""
2D premixed flame-vortex interaction with the all-Mach pressure projection: a lean hydrogen-air flame (h2o2.yaml), steady in a
0.77 m/s inflow, is struck by a counter-rotating vortex pair that drives into it from upstream.

The flame and the vortices are those of examples/2D_premixed_flame_vortex (hcid 271), at half its resolution: IC/ holds that
example's steady 1D flame interpolated to this grid. The fresh gas enters through a Dirichlet (-17) inlet that holds the
initial inflow state and leaves through a pressure outlet at 1 atm; y is periodic. Diffusion is mixture-averaged and the
reactions are operator-split (the projection's only reaction mode). SI units.

The vortices reach ~1.8 m/s, a Mach number of ~5e-3, and the projection steps at the explicit diffusion limit of the burned gas
rather than at its sound speed. --explicit runs the same case with the explicit solver, as a control.
"""

import argparse
import json
import os

parser = argparse.ArgumentParser(description="2D premixed flame-vortex interaction, all-Mach pressure projection")
parser.add_argument("--mfc", type=json.loads, default="{}", metavar="DICT", help="MFC's toolchain's internal state.")
parser.add_argument("--tstop", type=float, default=2.0e-3, help="end time [s] (default: %(default)s)")
parser.add_argument("--saves", type=int, default=20, help="number of outputs (default: %(default)s)")
parser.add_argument("--cfl", type=float, default=0.5, help="cfl_target (default: %(default)s)")
parser.add_argument("--explicit", action="store_true", help="explicit solver at the acoustic limit instead, as a control")
args = parser.parse_args()

L, N, p0 = 0.016, 512, 101325.0  # IC/ is tabulated on this x grid
current_dir = os.path.dirname(os.path.abspath(__file__))

case = {
    "run_time_info": "T",
    "x_domain%beg": 0.0,
    "x_domain%end": L,
    "y_domain%beg": 0.0,
    "y_domain%end": L / 2,
    "m": N - 1,
    "n": N // 2 - 1,
    "p": 0,
    "cfl_adap_dt": "T",
    "cfl_target": args.cfl,
    "n_start": 0,
    "t_stop": args.tstop,
    "t_save": args.tstop / args.saves,
    "model_eqns": 2,
    "num_fluids": 1,
    "num_patches": 1,
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 1,
    # Fresh gas in through the inlet, out at ambient pressure; periodic across the flame
    "bc_x%beg": -17,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "viscous": "T",
    "fluid_pp(1)%Re(1)": 200000,  # unused: a reacting mixture takes its viscosity from the transport model
    "chemistry": "T",
    "cantera_file": "h2o2.yaml",
    "chem_params%diffusion": "T",
    "chem_params%transport_model": 1,
    "chem_params%reactions": "T",
    "chem_params%reaction_substeps": 4,
    "chem_params%adap_substeps": "T",
    "chem_params%reaction_substeps_max": 64,
    "files_dir": os.path.join(current_dir, "IC"),
    "file_extension": "000000",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "chem_wrt_T": "T",
    "omega_wrt(3)": "T",
    "fd_order": 2,
    "fluid_pp(1)%gamma": 1.0 / (1.4 - 1.0),
    "fluid_pp(1)%eos": "stiffened_gas",
    "fluid_pp(1)%pi_inf": 0.0,
    # The steady 1D flame of IC/ along x, with the vortex pair added (hcid 271)
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%hcid": 271,
    "patch_icpp(1)%x_centroid": L / 2,
    "patch_icpp(1)%y_centroid": L / 4,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L / 2,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": 1.0,
    "patch_icpp(1)%alpha(1)": 1.0,
}
if not args.explicit:
    case.update({"proj_method": "T", "proj_max_acfl": 1.0e4, "bc_x%pres_out": p0})

print(json.dumps(case))
