#:include 'macros.fpp'

!> @brief This module is used to compute source terms for hypoelastic model
module m_surface_tension

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! ==========================================================================

    implicit none

    private; public :: s_compute_capillary_stress_tensor

contains

    subroutine s_compute_capillary_stress_tensor(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

    end subroutine


end module m_surface_tension