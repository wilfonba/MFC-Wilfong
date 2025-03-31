module m_signal_handler

    use iso_c_binding

    use m_global_parameters

    implicit none

    ! Make only the necessary procedures public
    private
    public :: s_setup_signal_handlers

    logical, public :: checkpoint_now

    ! Declare constants for signals
    integer(c_int), parameter :: SIGTERM = 15
    integer(c_int), parameter :: SIGUSR1 = 10

    interface
        ! Declare an interface to the C signal function
        function c_signal(sig, handler) bind(C, name="signal")
            import :: c_int, c_funptr
            implicit none
            integer(c_int) :: c_signal
            integer(c_int), value :: sig
            type(c_funptr), value :: handler
        end function c_signal
    end interface

contains

    subroutine s_setup_signal_handlers()
        implicit none
        integer(c_int) :: status

        ! Register signal handlers for SIGTERM and SIGUSR1
        status = c_signal(SIGTERM, c_funloc(my_signal_handler))
        if (status == 0 .and. proc_rank == 0) print *, "SIGTERM handler set."

        status = c_signal(SIGUSR1, c_funloc(my_signal_handler))
        if (status == 0 .and. proc_rank == 0) print *, "SIGUSR1 handler set."

    end subroutine s_setup_signal_handlers

    ! Signal handler function
    subroutine my_signal_handler(sig) bind(C)
        use iso_c_binding
        implicit none
        integer(c_int), value :: sig

        if (sig == SIGTERM) then
            if (proc_rank == 0) then
                print *, "Received SIGTERM! Cleaning up..."
            end if
            ! Add cleanup actions here
            checkpoint_now = .true.
        else if (sig == SIGUSR1) then
            if (proc_rank == 0) then
                print *, "Received SIGUSR1! Checkpointing..."
            end if
            checkpoint_now = .true.
        end if
    end subroutine my_signal_handler

end module m_signal_handler
