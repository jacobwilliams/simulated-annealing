


!TODO : - add a restart option. if the current point is > best found so far after so many iterations,
!         then reset back to the best one and continue (increase temp and/or vm when this happens?)
!       - input rand distributation option
!       - monotonic option - only accept a step if solution is better


!********************************************************************************
!>
!   simulated annealing is a global optimization method that distinguishes
!   between different local optima.  starting from an initial point, the
!   algorithm takes a step and the function is evaluated. when minimizing a
!   function, any downhill step is accepted and the process repeats from this
!   new point. an uphill step may be accepted. thus, it can escape from local
!   optima. this uphill decision is made by the metropolis criteria. as the
!   optimization process proceeds, the length of the steps decline and the
!   algorithm closes in on the global optimum. since the algorithm makes very
!   few assumptions regarding the function to be optimized, it is quite
!   robust with respect to non-quadratic surfaces. the degree of robustness
!   can be adjusted by the user. in fact, simulated annealing can be used as
!   a local optimizer for difficult functions.
!
!### Reference
!
!  * corana et al., "minimizing multimodal functions of continuous variables
!    with the "simulated annealing" algorithm", september 1987(vol. 13, no. 3, pp. 262-280),
!    acm transactions on mathematical software.
!  * goffe, ferrier and rogers, "global optimization of statistical functions
!    with simulated annealing", journal of econometrics, vol. 60, no. 1/2,
!    jan./feb. 1994, pp. 65-100.
!
!### History
!
!  * based on reference by bill goffe : 1/22/94, version: 3.2
!    See: https://www.netlib.org/opt/simann.f
!  * modifications by alan miller : fortran 90 version - 2 october 2013
!  * Jacob Williams, 8/26/2019 : modernized Fortran

    module simulated_anneal

    use iso_fortran_env, only: dp => real64, output_unit

    implicit none

    private

    !*******************************************************
    type,public :: simulated_annealing_type

        private

        integer :: n = 0              !! number of variables in the function to be optimized.
        logical :: maximize = .false. !! denotes whether the function should be maximized or minimized.
                                      !! a true value denotes maximization while a false value denotes
                                      !! minimization.  intermediate output (see iprint) takes this into
                                      !! account.
        real(dp) :: eps = 1.0e-9_dp   !! error tolerance for termination. if the final function
                                      !! values from the last neps temperatures differ from the
                                      !! corresponding value at the current temperature by less than
                                      !! eps and the final function value at the current temperature
                                      !! differs from the current optimal function value by less than
                                      !! eps, execution terminates and ier = 0 is returned.
        integer  :: ns = 20           !! number of cycles.  after ns function evaluations, each element of
                                      !! vm is adjusted according to the input `step_mode`
                                      !! the suggested value is 20.
        integer  :: nt = 100          !! number of iterations before temperature reduction. after
                                      !! nt*ns function evaluations, temperature (t) is changed
                                      !! by the factor rt.  value suggested by corana et al. is
                                      !! max(100, 5*n).  see goffe et al. for further advice.
        integer  :: neps = 4          !! number of final function values used to decide upon
                                      !! termination.  see eps.  suggested value is 4.
        integer  :: maxevl = 10000    !! the maximum number of function evaluations.  if it is
                                      !! exceeded, ier = 1.
        real(dp), dimension(:),allocatable :: lb        !! the lower bound for the allowable solution variables.
        real(dp), dimension(:),allocatable :: ub        !! the upper bound for the allowable solution variables.
                                                        !! if the algorithm chooses x(i) < lb(i) or x(i) > ub(i),
                                                        !! i = 1, n, a point is from inside is randomly selected. this
                                                        !! this focuses the algorithm on the region inside ub and lb.
                                                        !! unless the user wishes to concentrate the search to a particular
                                                        !! region, ub and lb should be set to very large positive
                                                        !! and negative values, respectively.  note that the starting
                                                        !! vector x should be inside this region.  also note that lb and
                                                        !! ub are fixed in position, while vm is centered on the last
                                                        !! accepted trial set of variables that optimizes the function.
        real(dp), dimension(:),allocatable :: c         !! vector that controls the step length adjustment.  the suggested
                                                        !! value for all elements is 2.0.
        integer :: iprint = 1           !! controls printing inside sa:
                                        !!
                                        !!  * 0 - nothing printed.
                                        !!  * 1 - function value for the starting value and
                                        !!    summary results before each temperature
                                        !!    reduction. this includes the optimal
                                        !!    function value found so far, the total
                                        !!    number of moves (broken up into uphill,
                                        !!    downhill, accepted and rejected), the
                                        !!    number of out of bounds trials, the
                                        !!    number of new optima found at this
                                        !!    temperature, the current optimal x and
                                        !!    the step length vm. note that there are
                                        !!    n*ns*nt function evalutations before each
                                        !!    temperature reduction. finally, notice is
                                        !!    is also given upon achieveing the termination
                                        !!    criteria.
                                        !!  * 2 - each new step length (vm), the current optimal
                                        !!    x (xopt) and the current trial x (x). this
                                        !!    gives the user some idea about how far x
                                        !!    strays from xopt as well as how vm is adapting
                                        !!    to the function.
                                        !!  * 3 - each function evaluation, its acceptance or
                                        !!    rejection and new optima. for many problems,
                                        !!    this option will likely require a small tree
                                        !!    if hard copy is used. this option is best
                                        !!    used to learn about the algorithm. a small
                                        !!    value for maxevl is thus recommended when
                                        !!    using iprint = 3.
                                        !!
                                        !! suggested value: 1
                                        !! note: for a given value of iprint, the lower valued
                                        !! options (other than 0) are utilized.
        integer :: iseed1 = 1234    !! the first seed for the random number generator.
        integer :: iseed2 = 5678    !! the second seed for the random number generator.
                                    !! different values for iseed1 and iseed2 will lead
                                    !! to an entirely different sequence of trial points
                                    !! and decisions on downhill moves (when maximizing).
                                    !! see goffe et al. on how this can be used to test
                                    !! the results of sa.
        integer   :: step_mode = 1  !! how to variable `vm` after `ns` cycles.
                                    !!
                                    !! * 1 : original method: adjust vm so that approximately
                                    !!   half of all evaluations are accepted
                                    !! * 2 : keep vm constant
                                    !! * 3 : adjust by a constant `vms` factor
        real(dp)  :: vms = 0.1_dp         !! for `step_mode=3`, the factor to adjust `vm`
        integer   :: iunit = output_unit  !! unit number for prints.

        procedure(sa_func),pointer :: fcn => null()





    contains

        private

        procedure,public :: initialize => initialize_sa
        procedure,public :: optimize   => sa
        procedure,public :: destroy    => destroy_sa

        procedure :: func
        procedure :: perturb_and_evaluate

    end type simulated_annealing_type
    !*******************************************************

    !*******************************************************
    abstract interface
        subroutine sa_func(me, x, f, istat)
            !! interface to function to be maximized/minimized
            import
            implicit none
            class(simulated_annealing_type),intent(inout) :: me
            real(dp),dimension(:),intent(in) :: x
            real(dp), intent(out)            :: f
            integer, intent(out)             :: istat !! status flag:
                                                      !!
                                                      !! * 0 : valid function evaluation
                                                      !! * -1 : invalid function evaluation.
                                                      !!   try a different input vector.
                                                      !! * -2 : stop the optimization process

        end subroutine sa_func
    end interface
    !*******************************************************

    public :: prtvec
    public :: dp

    contains
!********************************************************************************

!********************************************************************************
!>
!  Destructor.

    subroutine destroy_sa(me)

    implicit none

    class(simulated_annealing_type),intent(out) :: me

    end subroutine destroy_sa
!********************************************************************************

!********************************************************************************
!>
!  Initialize the class.
!
!  See the definition of [[simulated_annealing_type]] for the options.
!  Any options not set here will use the default values in the type.

    subroutine initialize_sa(me,fcn,n,lb,ub,c,&
                                maximize,eps,ns,nt,neps,maxevl,&
                                iprint,iseed1,iseed2,step_mode,vms,iunit)

    implicit none

    class(simulated_annealing_type),intent(inout) :: me
    procedure(sa_func)                          :: fcn
    integer, intent(in)                         :: n
    real(dp), dimension(n), intent(in)          :: lb
    real(dp), dimension(n), intent(in)          :: ub
    real(dp), dimension(n), intent(in),optional :: c
    logical, intent(in),optional                :: maximize
    real(dp), intent(in),optional               :: eps
    integer, intent(in),optional                :: ns
    integer, intent(in),optional                :: nt
    integer, intent(in),optional                :: neps
    integer, intent(in),optional                :: maxevl
    integer, intent(in),optional                :: iprint
    integer, intent(in),optional                :: iseed1
    integer, intent(in),optional                :: iseed2
    integer,intent(in),optional                 :: step_mode
    real(dp), intent(in), optional              :: vms
    integer, intent(in), optional               :: iunit

    call me%destroy()

    me%fcn => fcn
    me%n  = n
    me%lb = lb
    me%ub = ub

    ! check validity of bounds:
    if (any(me%lb > me%ub)) then
        error stop 'Error: at least one LB > UB.'
    end if

    allocate(me%c(n))
    if (present(c)) then
        me%c = c
    else
        me%c = 2.0_dp
    end if

    if (present(maximize) ) me%maximize = maximize
    if (present(eps)      ) me%eps = eps
    if (present(ns)       ) me%ns = ns
    if (present(nt)       ) me%nt = nt
    if (present(neps)     ) me%neps = neps
    if (present(maxevl)   ) me%maxevl = maxevl
    if (present(iprint)   ) me%iprint = iprint
    if (present(iseed1)   ) me%iseed1 = iseed1
    if (present(iseed2)   ) me%iseed2 = iseed2
    if (present(step_mode)) me%step_mode = step_mode
    if (present(vms)      ) me%vms = vms
    if (present(iunit)    ) me%iunit = iunit

    end subroutine initialize_sa
!********************************************************************************

!********************************************************************************
!>
!  Continuous simulated annealing global optimization algorithm
!
!### Brief Overview
!
!  sa tries to find the global optimum of an n dimensional function.
!  it moves both up and downhill and as the optimization process
!  proceeds, it focuses on the most promising area.
!
!  to start, it randomly chooses a trial point within the step length
!  vm (a vector of length n) of the user selected starting point. the
!  function is evaluated at this trial point and its value is compared
!  to its value at the initial point.
!
!  in a maximization problem, all uphill moves are accepted and the
!  algorithm continues from that trial point. downhill moves may be
!  accepted; the decision is made by the metropolis criteria. it uses t
!  (temperature) and the size of the downhill move in a probabilistic
!  manner. the smaller t and the size of the downhill move are, the more
!  likely that move will be accepted. if the trial is accepted, the
!  algorithm moves on from that point. if it is rejected, another point
!  is chosen instead for a trial evaluation.
!
!  each element of vm periodically adjusted so that half of all
!  function evaluations in that direction are accepted.
!
!  a fall in t is imposed upon the system with the rt variable by
!  t(i+1) = rt*t(i) where i is the ith iteration. thus, as t declines,
!  downhill moves are less likely to be accepted and the percentage of
!  rejections rise. given the scheme for the selection for vm, vm falls.
!  thus, as t declines, vm falls and sa focuses upon the most promising
!  area for optimization.
!
!  the parameter t is crucial in using sa successfully. it influences
!  vm, the step length over which the algorithm searches for optima. for
!  a small initial t, the step length may be too small; thus not enough
!  of the function might be evaluated to find the global optima. the user
!  should carefully examine vm in the intermediate output (set iprint =
!  1) to make sure that vm is appropriate. the relationship between the
!  initial temperature and the resulting step length is function
!  dependent.
!
!  to determine the starting temperature that is consistent with
!  optimizing a function, it is worthwhile to run a trial run first. set
!  rt = 1.5 and t = 1.0. with rt > 1.0, the temperature increases and vm
!  rises as well. then select the t that produces a large enough vm.
!
!### Notes
!  * as far as possible, the parameters here have the same name as in
!    the description of the algorithm on pp. 266-8 of corana et al.
!  * the suggested input values generally come from corana et al. to
!    drastically reduce runtime, see goffe et al., pp. 90-1 for
!    suggestions on choosing the appropriate rt and nt.
!
!### History
!  * differences compared to version 2.0:
!     1. if a trial is out of bounds, a point is randomly selected
!        from lb(i) to ub(i). unlike in version 2.0, this trial is
!        evaluated and is counted in acceptances and rejections.
!        all corresponding documentation was changed as well.
!  * differences compared to version 3.0:
!     1. if vm(i) > (ub(i) - lb(i)), vm is set to ub(i) - lb(i).
!        the idea is that if t is high relative to lb & ub, most
!        points will be accepted, causing vm to rise. but, in this
!        situation, vm has little meaning; particularly if vm is
!        larger than the acceptable region. setting vm to this size
!        still allows all parts of the allowable region to be selected.
!  * differences compared to version 3.1:
!     1. test made to see if the initial temperature is positive.
!     2. write statements prettied up.
!     3. references to paper updated.

    subroutine sa(me, x, rt, t, vm, xopt, fopt, nacc, nfcnev, ier)

    implicit none

    class(simulated_annealing_type),intent(inout) :: me
    real(dp), dimension(me%n), intent(inout)  :: x       !! on input: the starting values for the variables of the
                                                         !! function to be optimized. [Will be replaced by final point]
    real(dp), intent(in)                      :: rt      !! the temperature reduction factor.  the value suggested by
                                                         !! corana et al. is .85. see goffe et al. for more advice.
    real(dp), intent(inout)                   :: t       !! on input, the initial temperature. see goffe et al. for advice.
                                                         !! on output, the final temperature.
    real(dp), dimension(me%n), intent(inout)  :: vm      !! the step length vector. on input it should encompass the region of
                                                         !! interest given the starting value x.  for point x(i), the next
                                                         !! trial point is selected is from x(i) - vm(i)  to  x(i) + vm(i).
                                                         !! since vm is adjusted so that about half of all points are accepted,
                                                         !! the input value is not very important (i.e. is the value is off,
                                                         !! sa adjusts vm to the correct value).
    real(dp), dimension(me%n), intent(out)  :: xopt      !! the variables that optimize the function.
    real(dp), intent(out)                   :: fopt      !! the optimal value of the function.
    integer, intent(out)                    :: nacc      !! the number of accepted function evaluations.
    integer, intent(out)                    :: nfcnev    !! the total number of function evaluations. in a minor
                                                         !! point, note that the first evaluation is not used in the
                                                         !! core of the algorithm; it simply initializes the
                                                         !! algorithm.
    integer, intent(out)                    :: ier       !! the error return number:
                                                         !!
                                                         !!  * 0 - normal return; termination criteria achieved.
                                                         !!  * 1 - number of function evaluations (nfcnev) is
                                                         !!    greater than the maximum number (maxevl).
                                                         !!  * 3 - the initial temperature is not positive.
                                                         !!  * 4 - user stop in problem function.
                                                         !!  * 99 - should not be seen; only used internally.

    real(dp),dimension(me%n)    :: xp       !! perturbed `x` vector
    real(dp),dimension(me%neps) :: fstar    !! array of optimals found so far
    real(dp)                    :: f        !! function value
    real(dp)                    :: fp       !! function value
    real(dp)                    :: p        !! for metropolis criteria
    real(dp)                    :: pp       !! for metropolis criteria
    real(dp)                    :: ratio    !! ratio of `nacp` to `ns`
    integer                     :: nup      !! number of uphill steps
    integer                     :: ndown    !! number of accepted downhill steps
    integer                     :: nrej     !! number of rejected downhill steps
    integer                     :: nnew     !! new optimal for this temperature
    integer                     :: i, j, m  !! counter
    integer                     :: totmov   !! total moves
    integer                     :: nacp     !! number of accepts steps in a complete `ns` cycle
    logical                     :: quit     !! if the termination criteria was achieved
    integer                     :: unit     !! unit for printing

    ! initialize:
    unit   = me%iunit
    nacc   = 0
    nacp   = 0
    nfcnev = 0
    ier    = 99
    xopt   = x
    fstar  = huge(1.0_dp)
    fopt   = huge(1.0_dp)
    vm     = abs(vm) ! JW : just in case

    !  if the initial temperature is not positive, notify the user and
    !  return to the calling routine.
    if (t <= 0.0_dp) then
        if (me%iprint>0) write(unit,'(A)') 'Error: The initial temperature is not positive.'
        ier = 3
        return
    end if

    ! if the initial value is out of bounds, then just put
    ! the violated ones on the nearest bound.
    where (x > me%ub)
        x = me%ub
    else where (x < me%lb)
        x = me%lb
    end where

    !  initialize the random number generator
    call rand_init(me%iseed1,me%iseed2)

    ! evaluate the function with input x and return value as f.
    ! keep trying if necessary until there is a valid function
    ! evaluation.
    call me%perturb_and_evaluate(x,vm,xp,f,nfcnev,ier,first=.true.)
    select case (ier)
    case(1,4)
        return
    end select

    f = me%func(f)
    xopt = xp ! it may have been perturbed above
    fopt = f
    fstar(1) = f
    if (me%iprint >= 1) then
        write(unit, '(A)') '  '
        call prtvec(unit,x,me%n,'initial x')
        write(unit, '(A,1X,G25.18)')  '  initial f: ', me%func(f)
    end if

    do
        !  start the main loop. note that it terminates if (i) the algorithm
        !  succesfully optimizes the function or (ii) there are too many
        !  function evaluations (more than maxevl).
        nup = 0
        nrej = 0
        nnew = 0
        ndown = 0

        do m = 1, me%nt
            do j = 1, me%ns

                !  evaluate the function with the trial point xp and return as fp.
                call me%perturb_and_evaluate(x,vm,xp,fp,nfcnev,ier)
                select case (ier)
                case(1,4)
                    ! ...
                    ! ... need to set outputs to something ...
                    ! ...
                    fopt = me%func(fopt)
                    return
                end select

                fp = me%func(fp)
                if (me%iprint >= 3) then
                    write(unit,'(A)') ' '
                    call prtvec(unit,x,me%n,'current x')
                    write(unit,'(A,G25.18)') ' current f: ', me%func(f)
                    call prtvec(unit,xp,me%n,'trial x')
                    write(unit,'(A,G25.18)') ' resulting f: ', me%func(fp)
                end if

                !  accept the new point if the function value increases.
                if (fp >= f) then
                    if (me%iprint >= 3) write(unit,'(A)') '  point accepted'
                    x = xp
                    f = fp
                    nacc = nacc + 1
                    nacp = nacp + 1
                    nup = nup + 1

                    !  if greater than any other point, record as new optimum.
                    if (fp > fopt) then
                        if (me%iprint >= 3) write(unit,'(A)') '  new optimum'
                        xopt = xp
                        fopt = fp
                        nnew = nnew + 1
                    end if

                !  if the point is lower, use the metropolis criteria to decide on
                !  acceptance or rejection.
                else
                    p = exprep((fp - f)/t)
                    pp = rand()
                    if (pp < p) then
                        if (me%iprint >= 3) then
                            if (me%maximize) then
                                write(unit,'(A)')  '  though lower, point accepted'
                            else
                                write(unit,'(A)')  '  though higher, point accepted'
                            end if
                        end if
                        x = xp
                        f = fp
                        nacc = nacc + 1
                        nacp = nacp + 1
                        ndown = ndown + 1
                    else
                        nrej = nrej + 1
                        if (me%iprint >= 3) then
                            if (me%maximize) then
                                write(unit,'(A)') '  lower point rejected'
                            else
                                write(unit,'(A)') '  higher point rejected'
                            end if
                        end if
                    end if
                end if

            end do ! j - ns loop

            select case (me%step_mode)
            case(1)
                !  adjust vm so that approximately half of all evaluations are accepted.
                ratio = real(nacp,dp) /real(me%ns,dp)
                do i = 1, me%n
                    if (ratio > 0.6_dp) then
                        vm(i) = vm(i)*(1.0_dp + me%c(i)*(ratio - 0.6_dp)/0.4_dp)
                    else if (ratio < 0.4_dp) then
                        vm(i) = vm(i)/(1.0_dp + me%c(i)*((0.4_dp - ratio)/0.4_dp))
                    end if
                    if (vm(i) > (me%ub(i)-me%lb(i))) then
                        vm(i) = me%ub(i) - me%lb(i)
                    end if
                end do
                if (me%iprint >= 2) then
                    write(unit,'(/A)') '---------------------------------------------------'
                    write(unit,'(A)')  ' intermediate results after step length adjustment '
                    write(unit,'(A/)') '---------------------------------------------------'
                    call prtvec(unit, vm, me%n, 'new step length (vm)')
                    call prtvec(unit, xopt, me%n, 'current optimal x')
                    call prtvec(unit, x, me%n, 'current x')
                    write(unit,'(A)') ' '
                end if
            case (2)
                ! keep vm as is
            case (3)
                ! use the constant factor:
                vm = abs(me%vms) * vm
            case default
                error stop ' invalid value of step_mode'
            end select

            nacp = 0

        end do ! m - nt loop

        if (me%iprint >= 1) then
            totmov = nup + ndown + nrej
            write(unit,'(/A)')      '--------------------------------------------------------'
            write(unit,'(A)')       ' intermediate results before next temperature reduction '
            write(unit,'(A/)')      '--------------------------------------------------------'
            write(unit,'(A,G12.5)') '  current temperature:            ', t
            if (me%maximize) then
                write(unit,'(A,G25.18)') '  max function value so far:  ', fopt
                write(unit,'(A,I8)')     '  total moves:                ', totmov
                write(unit,'(A,I8)')     '     uphill:                  ', nup
                write(unit,'(A,I8)')     '     accepted downhill:       ', ndown
                write(unit,'(A,I8)')     '     rejected downhill:       ', nrej
                write(unit,'(A,I8)')     '  new maxima this temperature:', nnew
            else
                write(unit,'(A,G25.18)') '  min function value so far:  ', -fopt
                write(unit,'(A,I8)')     '  total moves:                ', totmov
                write(unit,'(A,I8)')     '     downhill:                ', nup
                write(unit,'(A,I8)')     '     accepted uphill:         ', ndown
                write(unit,'(A,I8)')     '     rejected uphill:         ', nrej
                write(unit,'(A,I8)')     '  new minima this temperature:', nnew
            end if
            call prtvec(unit,xopt, me%n, 'current optimal x')
            call prtvec(unit,vm,   me%n, 'step length (vm)')
            write(unit, '(A)') ' '
        end if

        !  check termination criteria.
        quit = .false.
        fstar(1) = f
        if ((fopt - fstar(1)) <= me%eps) quit = .true.
        do i = 1, me%neps
            if (abs(f - fstar(i)) > me%eps) quit = .false.
        end do

        !  terminate sa if appropriate.
        if (quit) then
            x = xopt
            ier = 0
            fopt = me%func(fopt)
            if (me%iprint >= 1) then
                write(unit,'(/A)') '----------------------------------------------'
                write(unit,'(A)')  '  sa achieved termination criteria. ier = 0.  '
                write(unit,'(A/)') '----------------------------------------------'
            end if
            return
        end if

        !  if termination criteria is not met, prepare for another loop.
        t = abs(rt)*t
        do i = me%neps, 2, -1
            fstar(i) = fstar(i-1)
        end do
        f = fopt
        x = xopt

    end do

    end subroutine sa
!********************************************************************************

!********************************************************************************
!>
!  Perturb the `x` vector and evaluate the function.
!
!  If the function evaluation is not valid, it will perturb
!  and try again. Until a valid function is obtained or the
!  maximum number of function evaluations is reached.

    subroutine perturb_and_evaluate(me,x,vm,xp,fp,nfcnev,ier,first)

    implicit none

    class(simulated_annealing_type),intent(inout) :: me
    real(dp),dimension(:),intent(in)  :: x       !! input optimization variable vector to perburb
    real(dp),dimension(:),intent(in)  :: vm      !! step length vector
    real(dp),dimension(:),intent(out) :: xp      !! the perturbed `x` value
    real(dp),intent(out)              :: fp      !! the value of the user function at `xp`
    integer,intent(inout)             :: nfcnev  !! total number of function evaluations
    integer,intent(inout)             :: ier     !! status output code
    logical,intent(in),optional       :: first   !! to use the input `x` the first time

    integer :: i         !! counter
    integer :: istat     !! user function status code
    logical :: first_try !! local copy of `first`
    real(dp) :: lower    !! lower bound to use for random interval
    real(dp) :: upper    !! upper bound to use for random interval

    if (present(first)) then
        first_try = first
    else
        first_try = .false.
    end if

    do

        if (first_try) then
            xp = x
            first_try = .false.
        else
            ! perturb all of them:
            do i = 1, me%n
                lower = max( me%lb(i), x(i) - vm(i) )
                upper = min( me%ub(i), x(i) + vm(i) )
                xp(i) = lower + (upper-lower)*rand()
            end do
        end if

        ! evaluate the function with the trial
        ! point xp and return as fp.
        call me%fcn(xp, fp, istat)

        ! function eval counter:
        nfcnev = nfcnev + 1

        !  if too many function evaluations occur, terminate the algorithm.
        if (nfcnev > me%maxevl) then
            if (me%iprint>0) then
                write(me%iunit, '(A)') ' too many function evaluations; consider'
                write(me%iunit, '(A)') ' increasing maxevl or eps, or decreasing'
                write(me%iunit, '(A)') ' nt or rt. these results are likely to be'
                write(me%iunit, '(A)') ' poor.'
            end if
            ier = 1
            return
        end if

        select case (istat)
        case(-2)
            ! user stop
            write(me%iunit, '(A)') ' user stopped in function.'
            ier = 4
            exit
        case(-1)
            ! try another one until we get a valid evaluation
            cycle
        case default
            ! continue
        end select

        exit ! finished

    end do

    end subroutine perturb_and_evaluate
!********************************************************************************

!********************************************************************************
!>
!  if the function is to be minimized, switch the sign of the function.
!  note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.

    pure function func(me,f)

    implicit none

    class(simulated_annealing_type),intent(in) :: me
    real(dp),intent(in) :: f
    real(dp) :: func

    if (me%maximize) then
        func = f
    else
        func = -f
    end if

    end function func
!********************************************************************************

!********************************************************************************
!>
!  this function replaces `exp()` to avoid underflow and overflow.
!
!  note that the maximum and minimum values of
!  exprep are such that they has no effect on the algorithm.
!
!### History
!  * Jacob Williams, 8/30/2019 : new version of this routine that uses IEEE flags.

    pure function exprep(x) result(f)

    use, intrinsic :: ieee_exceptions

    implicit none

    logical,dimension(2) :: flags
    type(ieee_flag_type),parameter,dimension(2) :: out_of_range = &
                                                    [ieee_overflow,ieee_underflow]

    real(dp), intent(in) :: x
    real(dp) :: f

    call ieee_set_halting_mode(out_of_range,.false.)

    f = exp(x)

    call ieee_get_flag(out_of_range,flags)
    if (any(flags)) then
        call ieee_set_flag(out_of_range,.false.)
        if (flags(1)) then
            f = huge(1.0_dp)
        else
            f = 0.0_dp
        end if
    end if

    end function exprep
!********************************************************************************

!********************************************************************************
!>
!  Initialize the random number generator.
!
!### Author
!  * Jacob Williams, 8/30/2019

    subroutine rand_init(seed1,seed2)

    implicit none

    integer,intent(in) :: seed1  !! the first seed for the random number generator.
    integer,intent(in) :: seed2  !! the second seed for the random number generator.

    integer,dimension(:),allocatable :: s
    integer :: n
    integer :: i

    if (seed1==0 .and. seed2==0) then
        call random_seed()
    else
        call random_seed(size = n)
        allocate(s(n))
        s(1) = seed1
        if (n>1) then
            s(2) = seed2
            do i = 3, n ! just in case
                s(i) = seed1 + seed2 + i
            end do
        end if
        call random_seed(put=s)
    end if

    end subroutine rand_init
!********************************************************************************

!********************************************************************************
!>
!  Get a new uniform random number.
!
!### Author
!  * Jacob Williams, 8/30/2019

    function rand() result(f)

    implicit none

    real :: f

    call random_number(f)

    end function rand
!********************************************************************************

!********************************************************************************
!>
!  this subroutine prints the double precision vector named vector.
!  elements 1 thru ncols will be printed. name is a character variable
!  that describes vector. note that if name is given in the call to
!  prtvec, it must be enclosed in quotes. if there are more than 10
!  elements in vector, 10 elements will be printed on each line.

    subroutine prtvec(iunit, vector, ncols, name)

    implicit none

    integer, intent(in)                     :: iunit
    integer, intent(in)                     :: ncols
    real(dp), dimension(ncols), intent(in)  :: vector
    character(len=*), intent(in)            :: name

    integer :: i, lines, ll

    write(iunit,'(/,25X,A)') trim(name)

    if (ncols > 10) then
        lines = int(ncols/10.0_dp)
        do i = 1, lines
            ll = 10*(i - 1)
            write(iunit,'(10(G12.5,1X))') vector(1+ll:10+ll)
        end do
        write(iunit,'(10(G12.5,1X))') vector(11+ll:ncols)
    else
        write(iunit,'(10(G12.5,1X))') vector(1:ncols)
    end if

    end subroutine prtvec
!********************************************************************************

!********************************************************************************
    end module simulated_anneal
!********************************************************************************