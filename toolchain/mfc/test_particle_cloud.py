import pytest

from mfc.case import Case
from mfc.case_validator import CaseConstraintError, CaseValidator
from mfc.common import MFCException


def _valid_cloud_params():
    return {
        "m": 31,
        "n": 31,
        "p": 31,
        "model_eqns": 2,
        "num_fluids": 1,
        "num_patches": 1,
        "t_step_start": 0,
        "t_step_stop": 1,
        "t_step_save": 1,
        "dt": 1.0e-6,
        "x_domain%beg": -1.0,
        "x_domain%end": 1.0,
        "y_domain%beg": -1.0,
        "y_domain%end": 1.0,
        "z_domain%beg": -1.0,
        "z_domain%end": 1.0,
        "bc_x%beg": -1,
        "bc_x%end": -1,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "bc_z%beg": -1,
        "bc_z%end": -1,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "avg_state": 2,
        "fd_order": 2,
        "ib": "T",
        "num_ibs": 0,
        "num_particle_clouds": 1,
        "particle_cloud(1)%cloud_geometry": 2,
        "particle_cloud(1)%packing_method": 1,
        "particle_cloud(1)%x_centroid": 0.0,
        "particle_cloud(1)%y_centroid": 0.0,
        "particle_cloud(1)%z_centroid": 0.0,
        "particle_cloud(1)%num_particles": 8,
        "particle_cloud(1)%radius": 0.04,
        "particle_cloud(1)%mass": 1.0,
        "particle_cloud(1)%min_spacing": 0.01,
        "particle_cloud(1)%shell_inner_radius": 0.1,
        "particle_cloud(1)%shell_outer_radius": 0.3,
    }


def test_hemi_shell_validator_accepts_feasible_shell():
    CaseValidator(_valid_cloud_params()).validate("simulation")


def test_hemi_shell_validator_rejects_negative_inner_radius():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_inner_radius": -0.01}
    with pytest.raises(CaseConstraintError, match="shell_inner_radius >= 0"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_shell_thinner_than_particle_diameter():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_outer_radius": 0.17}
    with pytest.raises(CaseConstraintError, match="shell_outer_radius > shell_inner_radius"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_shell_outside_domain():
    params = {**_valid_cloud_params(), "particle_cloud(1)%x_centroid": 0.9, "particle_cloud(1)%shell_outer_radius": 0.3}
    with pytest.raises(CaseConstraintError, match="x-extent must lie within x_domain"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_invalid_shell_axis():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_axis": 0}
    with pytest.raises(CaseConstraintError, match="shell_axis must be 1"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_applies_standoff_to_the_open_axis():
    # shell_axis = 1 opens the shell toward +x, so the flat face at x_centroid only needs one
    # particle radius of clearance from x_domain%beg, not the full outer radius.
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_axis": 1, "particle_cloud(1)%x_centroid": -0.95}
    CaseValidator(params).validate("simulation")


def test_void_fraction_resolves_to_a_box_particle_count():
    # A 3D box of 1.0 x 1.0 x 1.0 filled to 20% by spheres of radius 0.1 holds
    # 0.2 / ((4/3)*pi*0.1**3) = 47.7 -> 48 particles.
    params = {k: v for k, v in _valid_cloud_params().items() if k != "particle_cloud(1)%num_particles"}
    params.update(
        {
            "particle_cloud(1)%cloud_geometry": 1,
            "particle_cloud(1)%length_x": 1.0,
            "particle_cloud(1)%length_y": 1.0,
            "particle_cloud(1)%length_z": 1.0,
            "particle_cloud(1)%radius": 0.1,
            "particle_cloud(1)%void_fraction": 0.2,
        }
    )
    resolved = Case(params).get_parameters()
    assert resolved["particle_cloud(1)%num_particles"] == 48
    assert "particle_cloud(1)%void_fraction" not in resolved
    CaseValidator(resolved).validate("simulation")


def test_void_fraction_and_num_particles_together_are_rejected():
    params = {**_valid_cloud_params(), "particle_cloud(1)%void_fraction": 0.2}
    with pytest.raises(CaseConstraintError, match="not both"):
        CaseValidator(Case(params).get_parameters()).validate("simulation")


def test_neither_void_fraction_nor_num_particles_is_rejected():
    params = {k: v for k, v in _valid_cloud_params().items() if k != "particle_cloud(1)%num_particles"}
    with pytest.raises(CaseConstraintError, match="either num_particles or void_fraction"):
        CaseValidator(params).validate("simulation")


def test_void_fraction_outside_the_unit_interval_is_rejected():
    params = {k: v for k, v in _valid_cloud_params().items() if k != "particle_cloud(1)%num_particles"}
    params["particle_cloud(1)%void_fraction"] = 1.4
    with pytest.raises(MFCException, match="must be between 0 and 1"):
        Case(params)


def test_void_fraction_uses_the_shell_volume_for_hemisphere_clouds():
    # The 3D shell holds (2/3)*pi*(0.3**3 - 0.1**3) = 0.05445; filling 5% of it with spheres of
    # radius 0.04 (volume 2.681e-4) takes 10.2 -> 10 particles.
    params = {k: v for k, v in _valid_cloud_params().items() if k != "particle_cloud(1)%num_particles"}
    params["particle_cloud(1)%void_fraction"] = 0.05
    assert Case(params).get_parameters()["particle_cloud(1)%num_particles"] == 10


def test_validator_accepts_relaxation_packing_on_a_shell():
    # Relaxation handles both cloud geometries, unlike lattice packing.
    params = {**_valid_cloud_params(), "particle_cloud(1)%packing_method": 3}
    CaseValidator(params).validate("simulation")


def test_validator_rejects_unknown_packing_method():
    params = {**_valid_cloud_params(), "particle_cloud(1)%packing_method": 4}
    with pytest.raises(CaseConstraintError, match="packing_method must be 1"):
        CaseValidator(params).validate("simulation")
