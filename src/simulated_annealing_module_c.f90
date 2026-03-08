module simulated_annealing_module_c

    use iso_c_binding, only: c_double, c_int, c_char, c_null_char, &
                             c_intptr_t, c_ptr, c_loc, c_f_pointer, &
                             c_null_ptr, c_associated, c_bool, c_funptr, &
                             c_f_procpointer
    use simulated_annealing_module, only: simulated_annealing_type, wp => simann_wp

    implicit none

    public :: initialize_simulated_annealing
    public :: destroy_simulated_annealing
    public :: solve_simulated_annealing

    ! C function pointer types for callbacks
    abstract interface
        subroutine c_sa_func(ipointer, x, n, f, istat) bind(c)
            import :: c_double, c_int, c_intptr_t
            implicit none
            integer(c_intptr_t), value :: ipointer
            integer(c_int), value :: n
            real(c_double), dimension(n) :: x
            real(c_double) :: f
            integer(c_int) :: istat
        end subroutine c_sa_func

        subroutine c_sa_func_parallel_inputs(ipointer, n_inputs) bind(c)
            import :: c_int, c_intptr_t
            implicit none
            integer(c_intptr_t), value :: ipointer
            integer(c_int) :: n_inputs
        end subroutine c_sa_func_parallel_inputs

        subroutine c_sa_func_parallel_inputs_func(ipointer, x, n_inputs, n) bind(c)
            import :: c_double, c_int, c_intptr_t
            implicit none
            integer(c_intptr_t), value :: ipointer
            integer(c_int), value :: n_inputs
            integer(c_int), value :: n
            real(c_double), dimension(n_inputs, n) :: x
        end subroutine c_sa_func_parallel_inputs_func

        subroutine c_sa_func_parallel_output_func(ipointer, x, n, f, istat) bind(c)
            import :: c_double, c_int, c_intptr_t
            implicit none
            integer(c_intptr_t), value :: ipointer
            integer(c_int), value :: n
            real(c_double), dimension(n) :: x
            real(c_double) :: f
            integer(c_int) :: istat
        end subroutine c_sa_func_parallel_output_func
    end interface

    type,extends(simulated_annealing_type) :: c_sa_wrapper_type
        !! Wrapper type that encapsulates the C interface.
        !! extend the fortran type to hold the C function pointers.
        procedure(c_sa_func), pointer, nopass :: c_fcn_ptr => null()
        procedure(c_sa_func_parallel_inputs), pointer, nopass :: c_n_inputs_ptr => null()
        procedure(c_sa_func_parallel_inputs_func), pointer, nopass :: c_fcn_parallel_input_ptr => null()
        procedure(c_sa_func_parallel_output_func), pointer, nopass :: c_fcn_parallel_output_ptr => null()
        integer(c_intptr_t) :: ipointer = 0  !! pointer to this wrapper (for callbacks)
    end type c_sa_wrapper_type

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert an integer pointer to a [[c_sa_wrapper_type]] pointer.

    subroutine int_pointer_to_f_pointer(ipointer, p)

        integer(c_intptr_t), intent(in) :: ipointer !! integer pointer from C
        type(c_sa_wrapper_type), pointer :: p !! fortran pointer

        type(c_ptr) :: cp

        cp = transfer(ipointer, c_null_ptr)
        if (c_associated(cp)) then
            call c_f_pointer(cp, p)
        else
            p => null()
        end if

    end subroutine int_pointer_to_f_pointer

!*****************************************************************************************
!>
!  Wrapper procedures for C callbacks (module-level, work with wrapper type)

    subroutine fcn_wrapper(me, x, f, istat)
        class(simulated_annealing_type), intent(inout) :: me
        real(wp), dimension(:), intent(in) :: x
        real(wp), intent(out) :: f
        integer, intent(out) :: istat

        select type (me)
         type is (c_sa_wrapper_type)
          call me%c_fcn_ptr(me%ipointer, x, size(x), f, istat)
        end select

    end subroutine fcn_wrapper

    subroutine n_inputs_wrapper(me, n_inputs)
        class(simulated_annealing_type), intent(inout) :: me
        integer, intent(out) :: n_inputs

        select type (me)
         type is (c_sa_wrapper_type)
          call me%c_n_inputs_ptr(me%ipointer, n_inputs)
        end select

    end subroutine n_inputs_wrapper

    subroutine fcn_parallel_input_wrapper(me, x)
        class(simulated_annealing_type), intent(inout) :: me
        real(wp), dimension(:, :), intent(in) :: x

        select type (me)
         type is (c_sa_wrapper_type)
            call me%c_fcn_parallel_input_ptr(me%ipointer, x, size(x, 1), size(x, 2))
        end select

    end subroutine fcn_parallel_input_wrapper

    subroutine fcn_parallel_output_wrapper(me, x, f, istat)
        class(simulated_annealing_type), intent(inout) :: me
        real(wp), dimension(:), intent(out) :: x
        real(wp), intent(out) :: f
        integer, intent(out) :: istat

        select type (me)
         type is (c_sa_wrapper_type)
            call me%c_fcn_parallel_output_ptr(me%ipointer, x, size(x), f, istat)
        end select

    end subroutine fcn_parallel_output_wrapper

!*****************************************************************************************
!>
!  create a [[simulated_annealing_type]] from C

    subroutine initialize_simulated_annealing(ipointer, n, lb, ub, c, &
                                              maximize, eps, ns, nt, neps, maxevl, &
                                              iprint, iseed1, iseed2, step_mode, vms, iunit, &
                                              use_initial_guess, n_resets, &
                                              cooling_schedule, cooling_param, &
                                              optimal_f_specified, optimal_f, optimal_f_tol, &
                                              distribution_mode, dist_std_dev, &
                                              dist_scale, dist_shape, &
                                              fcn, n_inputs_to_send, fcn_parallel_input, fcn_parallel_output) &
      bind(C, name="initialize_simulated_annealing")

        integer(c_intptr_t), intent(out) :: ipointer
        integer(c_int), intent(in), value :: n
        real(c_double), dimension(n), intent(in) :: lb
        real(c_double), dimension(n), intent(in) :: ub
        real(c_double), dimension(n), intent(in) :: c
        logical(c_bool), intent(in), value :: maximize
        real(c_double), intent(in), value :: eps
        integer(c_int), intent(in), value :: ns
        integer(c_int), intent(in), value :: nt
        integer(c_int), intent(in), value :: neps
        integer(c_int), intent(in), value :: maxevl
        integer(c_int), intent(in), value :: iprint
        integer(c_int), intent(in), value :: iseed1
        integer(c_int), intent(in), value :: iseed2
        integer(c_int), intent(in), value :: step_mode
        real(c_double), intent(in), value :: vms
        integer(c_int), intent(in), value :: iunit
        logical(c_bool), intent(in), value :: use_initial_guess
        integer(c_int), intent(in), value :: n_resets
        integer(c_int), intent(in), value :: cooling_schedule
        real(c_double), intent(in), value :: cooling_param
        logical(c_bool), intent(in), value :: optimal_f_specified
        real(c_double), intent(in), value :: optimal_f
        real(c_double), intent(in), value :: optimal_f_tol
        integer(c_int), dimension(n), intent(in) :: distribution_mode
        real(c_double), dimension(n), intent(in) :: dist_std_dev
        real(c_double), dimension(n), intent(in) :: dist_scale
        real(c_double), dimension(n), intent(in) :: dist_shape
        type(c_funptr), intent(in), value :: fcn  !! C function pointer (can be C_NULL_FUNPTR)
        type(c_funptr), intent(in), value :: n_inputs_to_send  !! C function pointer (can be C_NULL_FUNPTR)
        type(c_funptr), intent(in), value :: fcn_parallel_input  !! C function pointer (can be C_NULL_FUNPTR)
        type(c_funptr), intent(in), value :: fcn_parallel_output  !! C function pointer (can be C_NULL_FUNPTR)

        type(c_sa_wrapper_type), pointer :: wrapper
        type(c_ptr) :: cp
        logical :: use_serial_mode, use_parallel_mode

        ! Allocate the wrapper
        allocate (wrapper)

        ! Store the wrapper pointer for use in callbacks
        cp = c_loc(wrapper)
        wrapper%ipointer = transfer(cp, 0_c_intptr_t)

        ! Convert C function pointers to Fortran procedure pointers
        use_serial_mode = .false.
        use_parallel_mode = .false.

        if (c_associated(fcn)) then
            call c_f_procpointer(fcn, wrapper%c_fcn_ptr)
            use_serial_mode = .true.
        end if

        if (c_associated(n_inputs_to_send) .and. &
            c_associated(fcn_parallel_input) .and. &
            c_associated(fcn_parallel_output)) then
            call c_f_procpointer(n_inputs_to_send, wrapper%c_n_inputs_ptr)
            call c_f_procpointer(fcn_parallel_input, wrapper%c_fcn_parallel_input_ptr)
            call c_f_procpointer(fcn_parallel_output, wrapper%c_fcn_parallel_output_ptr)
            use_parallel_mode = .true.
        end if

        ! Initialize the class with appropriate function pointers
        if (use_serial_mode) then
            call wrapper%initialize(n=n, lb=lb, ub=ub, c=c, &
                             maximize=logical(maximize), &
                             eps=eps, ns=ns, nt=nt, neps=neps, maxevl=maxevl, &
                             iprint=iprint, iseed1=iseed1, iseed2=iseed2, &
                             step_mode=step_mode, vms=vms, iunit=iunit, &
                             use_initial_guess=logical(use_initial_guess), &
                             n_resets=n_resets, &
                             cooling_schedule=cooling_schedule, &
                             cooling_param=cooling_param, &
                             optimal_f_specified=logical(optimal_f_specified), &
                             optimal_f=optimal_f, &
                             optimal_f_tol=optimal_f_tol, &
                             distribution_mode=distribution_mode, &
                             dist_std_dev=dist_std_dev, &
                             dist_scale=dist_scale, &
                             dist_shape=dist_shape, &
                             fcn=fcn_wrapper)
        else if (use_parallel_mode) then
            call wrapper%initialize(n=n, lb=lb, ub=ub, c=c, &
                             maximize=logical(maximize), &
                             eps=eps, ns=ns, nt=nt, neps=neps, maxevl=maxevl, &
                             iprint=iprint, iseed1=iseed1, iseed2=iseed2, &
                             step_mode=step_mode, vms=vms, iunit=iunit, &
                             use_initial_guess=logical(use_initial_guess), &
                             n_resets=n_resets, &
                             cooling_schedule=cooling_schedule, &
                             cooling_param=cooling_param, &
                             optimal_f_specified=logical(optimal_f_specified), &
                             optimal_f=optimal_f, &
                             optimal_f_tol=optimal_f_tol, &
                             distribution_mode=distribution_mode, &
                             dist_std_dev=dist_std_dev, &
                             dist_scale=dist_scale, &
                             dist_shape=dist_shape, &
                             n_inputs_to_send=n_inputs_wrapper, &
                             fcn_parallel_input=fcn_parallel_input_wrapper, &
                             fcn_parallel_output=fcn_parallel_output_wrapper)
        else
            error stop 'Error: either fcn (serial mode) or all of n_inputs_to_send, '//&
                       'fcn_parallel_input and fcn_parallel_output (parallel mode) must be provided.'
        end if

        ! Return converted pointer to C (pointer to the wrapper)
        ipointer = wrapper%ipointer

    end subroutine initialize_simulated_annealing

!*****************************************************************************************
!>
!  destroy a [[simulated_annealing_type]] from C

    subroutine destroy_simulated_annealing(ipointer) &
      bind(C, name="destroy_simulated_annealing")

        integer(c_intptr_t), intent(in) :: ipointer

        type(c_sa_wrapper_type), pointer :: wrapper

        call int_pointer_to_f_pointer(ipointer, wrapper)

        if (associated(wrapper)) then
         call wrapper%destroy()
            deallocate (wrapper)
        end if

    end subroutine destroy_simulated_annealing

!*****************************************************************************************
!>
!  solve optimization problem using simulated annealing from C

    subroutine solve_simulated_annealing(ipointer, n, x, rt, t, vm, xopt, fopt, nacc, nfcnev, ier) &
      bind(C, name="solve_simulated_annealing")

        integer(c_intptr_t), intent(in) :: ipointer
        integer(c_int), intent(in), value :: n
        real(c_double), dimension(n), intent(inout) :: x
        real(c_double), intent(in), value :: rt
        real(c_double), intent(inout) :: t
        real(c_double), dimension(n), intent(inout) :: vm
        real(c_double), dimension(n), intent(out) :: xopt
        real(c_double), intent(out) :: fopt
        integer(c_int), intent(out) :: nacc
        integer(c_int), intent(out) :: nfcnev
        integer(c_int), intent(out) :: ier

        type(c_sa_wrapper_type), pointer :: wrapper

        call int_pointer_to_f_pointer(ipointer, wrapper)

        if (associated(wrapper)) then
            call wrapper%optimize(x=x, rt=rt, t=t, vm=vm, &
                                  xopt=xopt, fopt=fopt, nacc=nacc, &
                                  nfcnev=nfcnev, ier=ier)
        else
            error stop 'error in solve_simulated_annealing: wrapper or SA class is not associated'
        end if

    end subroutine solve_simulated_annealing

end module simulated_annealing_module_c