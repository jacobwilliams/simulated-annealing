


!TODO : - maybe just force all values between the bounds all the time... rather than having special logic for that
!       - add a function flag to indicate a bad function eval, continue with best so far ...
!          -- has to count toward fcn eval counter
!       - maybe add a DMAX per variable per step ? or is that what the vm is doing ?
!       - if initial point is outside the bounds, just move it inside and continue ...


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

    abstract interface
        subroutine sa_func(n, x, f)
            !! interface to function to be maximized/minimized
            import :: dp
            implicit none
            integer, intent(in)              :: n
            real(dp),dimension(:),intent(in) :: x
            real(dp), intent(out)            :: f
        end subroutine sa_func
    end interface

    public :: sa  ! main routine

    public :: prtvec
    public :: dp

    contains
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

    subroutine sa(fcn, n, x, max, rt, eps, ns, nt, neps, maxevl, lb, ub, c, iprint,  &
                  iseed1, iseed2, t, vm, xopt, fopt, nacc, nfcnev, nobds, ier, iunit)

    implicit none

    procedure(sa_func)                      :: fcn    !! the function to maximize or minimize.
    integer, intent(in)                     :: n      !! number of variables in the function to be optimized.
    real(dp), dimension(n), intent(inout)   :: x      !! on input: the starting values for the variables of the
                                                      !! function to be optimized. [Will be replaced by final point]
    logical, intent(in)                     :: max    !! denotes whether the function should be maximized or minimized.
                                                      !! a true value denotes maximization while a false value denotes
                                                      !! minimization.  intermediate output (see iprint) takes this into
                                                      !! account.
    real(dp), intent(in)                    :: rt     !! the temperature reduction factor.  the value suggested by
                                                      !! corana et al. is .85. see goffe et al. for more advice.
    real(dp), intent(in)                    :: eps    !! error tolerance for termination. if the final function
                                                      !! values from the last neps temperatures differ from the
                                                      !! corresponding value at the current temperature by less than
                                                      !! eps and the final function value at the current temperature
                                                      !! differs from the current optimal function value by less than
                                                      !! eps, execution terminates and ier = 0 is returned.
    integer, intent(in)                     :: ns     !! number of cycles.  after ns*n function evaluations, each element of
                                                      !! vm is adjusted so that approximately half of all function evaluations
                                                      !! are accepted.  the suggested value is 20.
    integer, intent(in)                     :: nt     !! number of iterations before temperature reduction. after
                                                      !! nt*ns*n function evaluations, temperature (t) is changed
                                                      !! by the factor rt.  value suggested by corana et al. is
                                                      !! max(100, 5*n).  see goffe et al. for further advice.
    integer, intent(in)                     :: neps   !! number of final function values used to decide upon
                                                      !! termination.  see eps.  suggested value is 4.
    integer, intent(in)                     :: maxevl !! the maximum number of function evaluations.  if it is
                                                      !! exceeded, ier = 1.
    real(dp), dimension(n), intent(in)      :: lb     !! the lower bound for the allowable solution variables.
    real(dp), dimension(n), intent(in)      :: ub     !! the upper bound for the allowable solution variables.
                                                      !! if the algorithm chooses x(i) < lb(i) or x(i) > ub(i),
                                                      !! i = 1, n, a point is from inside is randomly selected. this
                                                      !! this focuses the algorithm on the region inside ub and lb.
                                                      !! unless the user wishes to concentrate the search to a particular
                                                      !! region, ub and lb should be set to very large positive
                                                      !! and negative values, respectively.  note that the starting
                                                      !! vector x should be inside this region.  also note that lb and
                                                      !! ub are fixed in position, while vm is centered on the last
                                                      !! accepted trial set of variables that optimizes the function.
    real(dp), dimension(n), intent(in)      :: c      !! vector that controls the step length adjustment.  the suggested
                                                      !! value for all elements is 2.0.
    integer, intent(in)                     :: iprint !! controls printing inside sa:
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
    integer, intent(in)                     :: iseed1 !! the first seed for the random number generator.
    integer, intent(in)                     :: iseed2 !! the second seed for the random number generator.
                                                      !! different values for iseed1 and iseed2 will lead
                                                      !! to an entirely different sequence of trial points
                                                      !! and decisions on downhill moves (when maximizing).
                                                      !! see goffe et al. on how this can be used to test
                                                      !! the results of sa.
    real(dp), intent(inout)                 :: t      !! on input, the initial temperature. see goffe et al. for advice.
                                                      !! on output, the final temperature.
    real(dp), dimension(n), intent(inout)   :: vm     !! the step length vector. on input it should encompass the region of
                                                      !! interest given the starting value x.  for point x(i), the next
                                                      !! trial point is selected is from x(i) - vm(i)  to  x(i) + vm(i).
                                                      !! since vm is adjusted so that about half of all points are accepted,
                                                      !! the input value is not very important (i.e. is the value is off,
                                                      !! sa adjusts vm to the correct value).
    real(dp), dimension(n), intent(out)     :: xopt   !! the variables that optimize the function.
    real(dp), intent(out)                   :: fopt   !! the optimal value of the function.
    integer, intent(out)                    :: nacc   !! the number of accepted function evaluations.
    integer, intent(out)                    :: nfcnev !! the total number of function evaluations. in a minor
                                                      !! point, note that the first evaluation is not used in the
                                                      !! core of the algorithm; it simply initializes the
                                                      !! algorithm.
    integer, intent(out)                    :: nobds  !! the total number of trial function evaluations that
                                                      !! would have been out of bounds of lb and ub. note that
                                                      !! a trial point is randomly selected between lb and ub.
    integer, intent(out)                    :: ier    !! the error return number:
                                                      !!
                                                      !!  * 0 - normal return; termination criteria achieved.
                                                      !!  * 1 - number of function evaluations (nfcnev) is
                                                      !!    greater than the maximum number (maxevl).
                                                      !!  * 2 - the starting value (x) is not inside the
                                                      !!    bounds (lb and ub).
                                                      !!  * 3 - the initial temperature is not positive.
                                                      !!  * 99 - should not be seen; only used internally.
    integer, intent(in), optional           :: iunit  !! unit number for prints. If not present, then
                                                      !! standard `output_unit` is used.

    real(dp)                 :: f, fp, p, pp, ratio
    real(dp),dimension(n)    :: xp
    real(dp),dimension(neps) :: fstar
    integer                  :: nup, ndown, nrej, nnew, lnobds, h, i, j, m, unit, totmov
    integer,dimension(n)     :: nacp
    logical                  :: quit

    if (present(iunit)) then
        unit = iunit
    else
        unit = output_unit
    end if

    !  initialize the random number generator
    call rmarin(iseed1,iseed2)

    !  set initial values.
    nacc = 0
    nobds = 0
    nfcnev = 0
    ier = 99
    xopt = x
    nacp = 0
    fstar = huge(1.0_dp)

    !  if the initial temperature is not positive, notify the user and
    !  return to the calling routine.
    if (t <= 0.0_dp) then
        if (iprint>0) write(unit,'(A)') 'The initial temperature is not positive. reset the variable t.'
        ier = 3
        return
    end if

    !  if the initial value is out of bounds, notify the user and return
    !  to the calling routine.
    do i = 1, n
        if ((x(i) > ub(i)) .or. (x(i) < lb(i))) then
            if (iprint>0) then
                write(unit, '(A)') '  the starting value (x) is outside the bounds'
                write(unit, '(A)') '  (lb and ub). execution terminated without any'
                write(unit, '(A)') '  optimization. respecify x, ub or lb so that'
                write(unit, '(A)') '  lb(i) < x(i) < ub(i), i = 1, n.'
            end if
            ier = 2
            return
        end if
    end do

    !  evaluate the function with input x and return value as f.
    call fcn(n, x, f)
    f = func(f,max)

    nfcnev = nfcnev + 1
    fopt = f
    fstar(1) = f
    if (iprint >= 1) then
        write(*, '(A)') '  '
        call prtvec(unit,x,n,'initial x')
        write(unit, '(A,1X,G25.18)')  '  initial f: ', func(f,max)
    end if

    do
        !  start the main loop. note that it terminates if (i) the algorithm
        !  succesfully optimizes the function or (ii) there are too many
        !  function evaluations (more than maxevl).
        nup = 0
        nrej = 0
        nnew = 0
        ndown = 0
        lnobds = 0

        do m = 1, nt
            do j = 1, ns
                do h = 1, n

                    !  generate xp, the trial value of x. note use of vm to choose xp.
                    do i = 1, n
                        if (i == h) then
                            xp(i) = x(i) + (ranmar()*2.0_dp - 1.0_dp) * vm(i)
                        else
                            xp(i) = x(i)
                        end if

                        !  if xp is out of bounds, select a point in bounds for the trial.
                        if ((xp(i) < lb(i)) .or. (xp(i) > ub(i))) then
                            xp(i) = lb(i) + (ub(i) - lb(i))*ranmar()
                            lnobds = lnobds + 1
                            nobds = nobds + 1
                            if (iprint >= 3) then
                                write(unit, '(A)') '  '
                                call prtvec(unit,x, n, 'current x')
                                write(unit, '(A,1X,G25.18)') '  current f: ', func(f,max)
                                call prtvec(unit,xp, n, 'trial x')
                                write(unit, '(A)') '  point rejected since out of bounds'
                            end if
                        end if
                    end do

                    !  evaluate the function with the trial point xp and return as fp.
                    call fcn(n, xp, fp)
                    fp = func(fp,max)
                    nfcnev = nfcnev + 1
                    if (iprint >= 3) then
                        write(*,'(A)') ' '
                        call prtvec(unit,x,n,'current x')
                        write(*,'(A,G25.18)') ' current f: ', func(f,max)
                        call prtvec(unit,xp,n,'trial x')
                        write(*,'(A,G25.18)') ' resulting f: ', func(fp,max)
                    end if

                    !  if too many function evaluations occur, terminate the algorithm.
                    if (nfcnev >= maxevl) then
                        if (iprint>0) then
                            write(unit, '(A)') '  too many function evaluations; consider'
                            write(unit, '(A)') '  increasing maxevl or eps, or decreasing'
                            write(unit, '(A)') '  nt or rt. these results are likely to be'
                            write(unit, '(A)') '  poor.'
                        end if
                        fopt = func(fopt,max)
                        ier = 1
                        return
                    end if

                    !  accept the new point if the function value increases.
                    if (fp >= f) then
                        if (iprint >= 3) then
                            write(unit,'(A)') '  point accepted'
                        end if
                        x = xp
                        f = fp
                        nacc = nacc + 1
                        nacp(h) = nacp(h) + 1
                        nup = nup + 1

                        !  if greater than any other point, record as new optimum.
                        if (fp > fopt) then
                            if (iprint >= 3) then
                                write(unit,'(A)') '  new optimum'
                            end if
                            xopt = xp
                            fopt = fp
                            nnew = nnew + 1
                        end if

                    !  if the point is lower, use the metropolis criteria to decide on
                    !  acceptance or rejection.
                    else
                        p = exprep((fp - f)/t)
                        pp = ranmar()
                        if (pp < p) then
                            if (iprint >= 3) then
                                if (max) then
                                    write(unit,'(A)')  '  though lower, point accepted'
                                else
                                    write(unit,'(A)')  '  though higher, point accepted'
                                end if
                            end if
                            x = xp
                            f = fp
                            nacc = nacc + 1
                            nacp(h) = nacp(h) + 1
                            ndown = ndown + 1
                        else
                            nrej = nrej + 1
                            if (iprint >= 3) then
                                if (max) then
                                    write(unit,'(A)') '  lower point rejected'
                                else
                                    write(unit,'(A)') '  higher point rejected'
                                end if
                            end if
                        end if
                    end if

                end do
            end do

            !  adjust vm so that approximately half of all evaluations are accepted.
            do i = 1, n
                ratio = real(nacp(i),dp) /real(ns,dp)
                if (ratio > 0.6_dp) then
                    vm(i) = vm(i)*(1.0_dp + c(i)*(ratio - 0.6_dp)/0.4_dp)
                else if (ratio < 0.4_dp) then
                    vm(i) = vm(i)/(1.0_dp + c(i)*((0.4_dp - ratio)/0.4_dp))
                end if
                if (vm(i) > (ub(i)-lb(i))) then
                    vm(i) = ub(i) - lb(i)
                end if
            end do

            if (iprint >= 2) then
                write(unit,'(/A)') '---------------------------------------------------'
                write(unit,'(A)')  ' intermediate results after step length adjustment '
                write(unit,'(A/)') '---------------------------------------------------'
                call prtvec(unit, vm, n, 'new step length (vm)')
                call prtvec(unit, xopt, n, 'current optimal x')
                call prtvec(unit, x, n, 'current x')
                write(unit,'(A)') ' '
            end if

            nacp = 0

        end do

        if (iprint >= 1) then
            totmov = nup + ndown + nrej
            write(unit,'(/A)') '--------------------------------------------------------'
            write(*,'(A)')     ' intermediate results before next temperature reduction '
            write(unit,'(A/)') '--------------------------------------------------------'
            write(*,'(A,G12.5)') '  current temperature:            ', t
            if (max) then
                write(*,'(A,G25.18)') '  max function value so far:  ', fopt
                write(*,'(A,I8)')     '  total moves:                ', totmov
                write(*,'(A,I8)')     '     uphill:                  ', nup
                write(*,'(A,I8)')     '     accepted downhill:       ', ndown
                write(*,'(A,I8)')     '     rejected downhill:       ', nrej
                write(*,'(A,I8)')     '  out of bounds trials:       ', lnobds
                write(*,'(A,I8)')     '  new maxima this temperature:', nnew
            else
                write(*,'(A,G25.18)') '  min function value so far:  ', -fopt
                write(*,'(A,I8)')     '  total moves:                ', totmov
                write(*,'(A,I8)')     '     downhill:                ', nup
                write(*,'(A,I8)')     '     accepted uphill:         ', ndown
                write(*,'(A,I8)')     '     rejected uphill:         ', nrej
                write(*,'(A,I8)')     '  trials out of bounds:       ', lnobds
                write(*,'(A,I8)')     '  new minima this temperature:', nnew
            end if
            call prtvec(unit,xopt, n, 'current optimal x')
            call prtvec(unit,vm, n, 'step length (vm)')
            write(*, '(A)') ' '
        end if

        !  check termination criteria.
        quit = .false.
        fstar(1) = f
        if ((fopt - fstar(1)) <= eps) quit = .true.
        do i = 1, neps
            if (abs(f - fstar(i)) > eps) quit = .false.
        end do

        !  terminate sa if appropriate.
        if (quit) then
            x = xopt
            ier = 0
            fopt = func(fopt,max)
            if (iprint >= 1) then
                write(unit,'(/A)') '----------------------------------------------'
                write(unit, '(A)') '  sa achieved termination criteria. ier = 0.  '
                write(unit,'(A/)') '----------------------------------------------'
            end if
            return
        end if

        !  if termination criteria is not met, prepare for another loop.
        t = rt*t
        do i = neps, 2, -1
            fstar(i) = fstar(i-1)
        end do
        f = fopt
        x = xopt

    end do

    end subroutine sa
!********************************************************************************

!********************************************************************************
!>
!  if the function is to be minimized, switch the sign of the function.
!  note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.

    pure function func(f,max)

    implicit none

    real(dp),intent(in) :: f
    logical,intent(in)  :: max
    real(dp)            :: func

    if (max) then
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

    subroutine rmarin(seed1,seed2)

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

    end subroutine rmarin
!********************************************************************************

!********************************************************************************
!>
!  Get a new uniform random number.
!
!### Author
!  * Jacob Williams, 8/30/2019

    function ranmar() result(f)

    implicit none

    real :: f

    call random_number(f)

    end function ranmar
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