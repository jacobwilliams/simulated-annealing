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
!    with the "simulated annealing" algorithm", september 1987
!    (vol. 13, no. 3, pp. 262-280),
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
!
!### TODO
!  * input rand distributation option
!  * a way to specify that some variables are not to be changed... this could be part of the distributation
!    selection (a constant option). each variable can have a different distribution.
!  * get ideas from: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html

module simulated_annealing_module

   use iso_fortran_env

   implicit none

   private

#ifdef REAL32
   integer,parameter :: wp = real32   !! real kind used by this module [4 bytes]
#elif REAL64
   integer,parameter :: wp = real64   !! real kind used by this module [8 bytes]
#elif REAL128
   integer,parameter :: wp = real128  !! real kind used by this module [16 bytes]
#else
   integer,parameter :: wp = real64   !! real kind used by this module [8 bytes]
#endif
   integer,parameter,public :: simann_wp = wp   !! for exporting from the module

   real(wp),parameter :: pi = acos(-1.0_wp)  !! value of pi for this module's real kind
   real(wp),parameter :: twopi = 2.0_wp * pi

   integer,parameter,public :: n_distribution_modes = 5  !! number of modes
   integer,parameter,public :: sa_mode_uniform = 1
   integer,parameter,public :: sa_mode_normal = 2
   integer,parameter,public :: sa_mode_cauchy = 3
   integer,parameter,public :: sa_mode_triangular = 4
   integer,parameter,public :: sa_mode_bipareto = 5

   !*******************************************************
   type,public :: simulated_annealing_type

      private

      integer :: n = 0              !! number of variables in the function to be optimized.
      logical :: maximize = .false. !! denotes whether the function should be maximized or minimized.
                                    !! a true value denotes maximization while a false value denotes
                                    !! minimization.  intermediate output (see iprint) takes this into
                                    !! account.
      real(wp) :: eps = 1.0e-9_wp   !! error tolerance for termination. if the final function
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

      logical :: use_initial_guess = .true. !! if false, the initial guess is ignored and a
                                            !! random point in the bounds is used for the first
                                            !! function evaluation

      integer :: n_resets = 2 !! number of times to run the main loop (must be >=1)

      real(wp), dimension(:),allocatable :: lb   !! the lower bound for the allowable solution variables.
      real(wp), dimension(:),allocatable :: ub   !! the upper bound for the allowable solution variables.
                                                 !! if the algorithm chooses x(i) < lb(i) or x(i) > ub(i),
                                                 !! i = 1, n, a point is from inside is randomly selected. this
                                                 !! this focuses the algorithm on the region inside ub and lb.
                                                 !! unless the user wishes to concentrate the search to a particular
                                                 !! region, ub and lb should be set to very large positive
                                                 !! and negative values, respectively.  note that the starting
                                                 !! vector x should be inside this region.  also note that lb and
                                                 !! ub are fixed in position, while vm is centered on the last
                                                 !! accepted trial set of variables that optimizes the function.
      real(wp), dimension(:),allocatable :: c  !! vector that controls the step length adjustment.  the suggested
                                               !! value for all elements is 2.0.
      integer :: iprint = 1     !! controls printing inside [[sa]]:
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
                                  !! the results of [[sa]].
      integer   :: step_mode = 1  !! how to vary `vm` after `ns` cycles.
                                  !!
                                  !! * 1 : original method: adjust vm so that approximately
                                  !!   half of all evaluations are accepted
                                  !! * 2 : keep vm constant
                                  !! * 3 : adjust by a constant `vms` factor
      real(wp)  :: vms = 0.1_wp         !! for `step_mode=3`, the factor to adjust `vm`
      integer   :: iunit = output_unit  !! unit number for prints.

      ! cooling schedule options:
      integer   :: cooling_schedule = 1 !! temperature reduction schedule:
                                        !!
                                        !! * 1 : geometric (default): T(k+1) = rt * T(k)
                                        !! * 2 : fast annealing (Cauchy): T(k) = T0 / (1 + k)
                                        !! * 3 : huang: T(k) = T0 / (1 + cooling_param * k)^n
                                        !! * 4 : boltzmann: T(k) = T0 / log(1 + k + exp(1))
                                        !! * 5 : logarithmic: T(k) = T0 / (1 + cooling_param * log(1 + k))
      real(wp)  :: cooling_param = 1.0_wp !! parameter for cooling schedules 3 and 5. suggested value: 1.0

      ! option to specify the known answer, so the solver will stop if it finds it:
      logical :: optimal_f_specified = .false. !! if the optional `f` value is known,
                                               !! it can be specified by `optimal_f`.
      real(wp) :: optimal_f = 0.0_wp !! if `optimal_f_specified=True` the solver will stop
                                     !! if this value is achieved.
      real(wp) :: optimal_f_tol = 0.0_wp !! absolute tolerance for the `optimal_f` check

      ! distribution selection for perturbations (per-variable):
      integer, dimension(:), allocatable :: distribution_mode !! distribution to use for perturbations for each variable:
                                                              !!
                                                              !! * `sa_mode_uniform` : uniform (default)
                                                              !! * `sa_mode_normal` : normal (Gaussian)
                                                              !! * `sa_mode_cauchy` : cauchy
                                                              !! * `sa_mode_triangular` : triangular
                                                              !! * `sa_mode_bipareto` : bipareto (two-sided pareto)

      ! distribution parameters (per-variable, used depending on distribution_mode):
      real(wp), dimension(:), allocatable :: dist_std_dev   !! standard deviation for normal/truncated_normal
      real(wp), dimension(:), allocatable :: dist_scale     !! scale parameter for cauchy/bipareto distributions
      real(wp), dimension(:), allocatable :: dist_shape     !! shape parameter for bipareto distribution

      ! serial function evaluation:
      procedure(sa_func),pointer :: fcn => null()  !! the user's function

      ! parallel function evaluation (optional):
      logical :: parallel_mode = .false. !! if true, the user wants to evaluate multiple points in parallel (e.g. on a GPU or with OpenMP). if false, the user will evaluate one point at a time. if true, the user must provide all of n_inputs_to_send, fcn_parallel_input and fcn_parallel_output.
      procedure(sa_func_parallel_inputs),pointer :: n_inputs_to_send => null()  !! if the user wants to evaluate multiple points in parallel, this function should return the number of input points that will be sent for evaluation at a time.
      procedure(sa_func_parallel_inputs_func),pointer :: fcn_parallel_input => null()  !! this function receives the input points for evaluation
      procedure(sa_func_parallel_output_func),pointer :: fcn_parallel_output => null()  !! this function receives the output points from the parallel evaluation. The `x` here will be one of the ones sent via `sa_func_parallel_inputs_func`, The algorithm will only process one at the time.

   contains

      private

      procedure,public :: initialize => initialize_sa
      procedure,public :: optimize   => sa
      procedure,public :: destroy    => destroy_sa

      procedure :: func
      procedure :: perturb_and_evaluate
      procedure,public :: perturb_variable  !! this is public so we can use it in the tests

   end type simulated_annealing_type
   !*******************************************************

   !*******************************************************
   abstract interface
      subroutine sa_func(me, x, f, istat)
         !! interface to function to be maximized/minimized
         import :: wp, simulated_annealing_type
         implicit none
         class(simulated_annealing_type),intent(inout) :: me
         real(wp),dimension(:),intent(in) :: x
         real(wp), intent(out)            :: f
         integer, intent(out)             :: istat !! status flag:
                                                   !!
                                                   !! * 0 : valid function evaluation
                                                   !! * -1 : invalid function evaluation.
                                                   !!   try a different input vector.
                                                   !! * -2 : stop the optimization process
      end subroutine sa_func

      subroutine sa_func_parallel_inputs(me, n_inputs)
         !! interface for parallel function evaluations.
         !! This is used when the user wants to evaluate
         !! multiple points in parallel (e.g. on a GPU or with OpenMP)
         !! This function should return the number of input points to be received
         import :: simulated_annealing_type
         implicit none
         class(simulated_annealing_type),intent(inout) :: me
         integer,intent(out) :: n_inputs !! number of inputs that can be received
      end subroutine sa_func_parallel_inputs

      subroutine sa_func_parallel_inputs_func(me, x)
         !! interface for parallel function evaluations.
         !! This function receives the input points for evaluation
         import :: wp, simulated_annealing_type
         implicit none
         class(simulated_annealing_type),intent(inout) :: me
         real(wp),dimension(:,:),intent(in) :: x  !! input points (n_inputs x n)
      end subroutine sa_func_parallel_inputs_func

      subroutine sa_func_parallel_output_func(me, x, f, istat)
         !! interface for parallel function evaluations.
         !! This function receives the output points from the
         !! parallel evaluation. Returns both the `x` and `f` values
         !! for a completed evaluation (typically from a queue).
         import :: wp, simulated_annealing_type
         implicit none
         class(simulated_annealing_type),intent(inout) :: me
         real(wp),dimension(:),intent(out) :: x
         real(wp), intent(out)            :: f
         integer, intent(out)             :: istat !! status flag:
                                                   !!
                                                   !! * 0 : valid function evaluation
                                                   !! * -1 : invalid function evaluation.
                                                   !!   try a different input vector.
                                                   !! * -2 : stop the optimization process
      end subroutine sa_func_parallel_output_func

      function dist_func(me, lower, upper) result(r)
         !! interface for distribution functions used for perturbations
         import :: wp, simulated_annealing_type
         implicit none
         class(simulated_annealing_type),intent(inout) :: me
         real(wp),intent(in) :: lower  !! lower bound
         real(wp),intent(in) :: upper  !! upper bound
         real(wp) :: r  !! random value within [lower, upper]
      end function dist_func
   end interface
   !*******************************************************

   public :: print_vector

contains
!********************************************************************************

!********************************************************************************
!>
!  Destructor.

   subroutine destroy_sa(me)

      class(simulated_annealing_type),intent(out) :: me

   end subroutine destroy_sa
!********************************************************************************

!********************************************************************************
!>
!  Initialize the class.
!
!  See the definition of [[simulated_annealing_type]] for the options.
!  Any options not set here will use the default values in the type.

   subroutine initialize_sa(me,n,lb,ub,c,&
                            maximize,eps,ns,nt,neps,maxevl,&
                            iprint,iseed1,iseed2,step_mode,vms,iunit,&
                            use_initial_guess,n_resets,&
                            cooling_schedule,cooling_param,&
                            optimal_f_specified,optimal_f,optimal_f_tol,&
                            distribution_mode,dist_std_dev,&
                            dist_scale,dist_shape,&
                            fcn,n_inputs_to_send,fcn_parallel_input,fcn_parallel_output)

      class(simulated_annealing_type),intent(inout) :: me
      integer, intent(in)                           :: n
      real(wp), dimension(n), intent(in)            :: lb
      real(wp), dimension(n), intent(in)            :: ub
      real(wp), dimension(n), intent(in),optional   :: c
      logical, intent(in),optional                  :: maximize
      real(wp), intent(in),optional                 :: eps
      integer, intent(in),optional                  :: ns
      integer, intent(in),optional                  :: nt
      integer, intent(in),optional                  :: neps
      integer, intent(in),optional                  :: maxevl
      integer, intent(in),optional                  :: iprint
      integer, intent(in),optional                  :: iseed1
      integer, intent(in),optional                  :: iseed2
      integer,intent(in),optional                   :: step_mode
      real(wp), intent(in), optional                :: vms
      integer, intent(in), optional                 :: iunit
      logical, intent(in), optional                 :: use_initial_guess
      integer, intent(in), optional                 :: n_resets
      integer, intent(in), optional                 :: cooling_schedule   !! temperature reduction schedule (1-5)
      real(wp), intent(in), optional                :: cooling_param      !! parameter for cooling schedules 3 and 5
      logical, intent(in), optional  :: optimal_f_specified !! if the optional `f` value is known,
                                                            !! it can be specified by `optimal_f`.
                                                            !! [Default is False]
      real(wp), intent(in), optional  :: optimal_f          !! if `optimal_f_specified=True` the solver will stop
                                                            !! if this value is achieved.
      real(wp), intent(in), optional  :: optimal_f_tol      !! absolute tolerance for the `optimal_f` check
      integer, dimension(:), intent(in), optional   :: distribution_mode  !! distribution for perturbations (per variable):
                                                                          !!
                                                                          !! * `sa_mode_uniform` : uniform (default)
                                                                          !! * `sa_mode_normal` : normal (Gaussian)
                                                                          !! * `sa_mode_cauchy` : cauchy
                                                                          !! * `sa_mode_triangular` : triangular
                                                                          !! * `sa_mode_bipareto` : bipareto (two-sided pareto)
                                                                          !!
                                                                          !! Can be scalar (broadcast to all) or array(n)
      real(wp), dimension(:), intent(in), optional  :: dist_std_dev       !! std dev for normal/truncated_normal (per variable or scalar)
      real(wp), dimension(:), intent(in), optional  :: dist_scale         !! scale for cauchy/pareto (per variable or scalar)
      real(wp), dimension(:), intent(in), optional  :: dist_shape         !! shape for pareto (per variable or scalar)
      procedure(sa_func),optional                       :: fcn                  !! the user's function for scalar evaluation.
                                                                                !! if not present, the user must provide all of
                                                                                !! n_inputs_to_send, fcn_parallel_input and fcn_parallel_output
                                                                                !! for parallel evaluation.
      procedure(sa_func_parallel_inputs),optional       :: n_inputs_to_send     !! if the user wants to evaluate multiple points in parallel,
                                                                                !! this function should return the number of input points
                                                                                !! that will be sent for evaluation at a time.
      procedure(sa_func_parallel_inputs_func),optional  :: fcn_parallel_input   !! this function receives the input points for evaluation
      procedure(sa_func_parallel_output_func),optional  :: fcn_parallel_output  !! this function receives the output points from the parallel
                                                                                !! evaluation. The `x` here will be one of the ones sent via
                                                                                !! `sa_func_parallel_inputs_func`,
                                                                                !! The algorithm will only process one at the time.

      call me%destroy()

      if (present(fcn)) then
         ! serial function evaluation:
         me%fcn => fcn
         me%parallel_mode = .false.
      else if (present(n_inputs_to_send) .and. &
         present(fcn_parallel_input) .and. &
         present(fcn_parallel_output)) then
         ! parallel function evaluation:
         me%n_inputs_to_send    => n_inputs_to_send
         me%fcn_parallel_input  => fcn_parallel_input
         me%fcn_parallel_output => fcn_parallel_output
         me%parallel_mode = .true.
      else
         error stop 'Error: either fcn (serial mode) or all of n_inputs_to_send, '//&
                    'fcn_parallel_input and fcn_parallel_output (parallel mode) must be provided.'
      end if

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
         me%c = 2.0_wp
      end if

      if (present(maximize)            ) me%maximize = maximize
      if (present(eps)                 ) me%eps = abs(eps)
      if (present(ns)                  ) me%ns = ns
      if (present(nt)                  ) me%nt = nt
      if (present(neps)                ) me%neps = neps
      if (present(maxevl)              ) me%maxevl = maxevl
      if (present(iprint)              ) me%iprint = iprint
      if (present(iseed1)              ) me%iseed1 = iseed1
      if (present(iseed2)              ) me%iseed2 = iseed2
      if (present(step_mode)           ) me%step_mode = step_mode
      if (present(vms)                 ) me%vms = vms
      if (present(iunit)               ) me%iunit = iunit
      if (present(use_initial_guess)   ) me%use_initial_guess = use_initial_guess
      if (present(n_resets)            ) me%n_resets = n_resets
      if (present(cooling_schedule)    ) me%cooling_schedule = cooling_schedule
      if (present(cooling_param)       ) me%cooling_param = cooling_param
      if (present(optimal_f_specified) ) me%optimal_f_specified = optimal_f_specified
      if (present(optimal_f)           ) me%optimal_f = optimal_f
      if (present(optimal_f_tol)       ) me%optimal_f_tol = abs(optimal_f_tol)

      ! allocate and set distribution parameters (per-variable):
      allocate(me%distribution_mode(n))
      allocate(me%dist_std_dev(n))
      allocate(me%dist_scale(n))
      allocate(me%dist_shape(n))

      ! set default values (uniform distribution with default parameters):
      me%distribution_mode = sa_mode_uniform
      me%dist_std_dev = 1.0_wp
      me%dist_scale = 1.0_wp
      me%dist_shape = 1.0_wp

      ! override with user-supplied values (broadcast if scalar):
      if (present(distribution_mode)) then
         if (size(distribution_mode) == 1) then
            me%distribution_mode = distribution_mode(1)  ! broadcast scalar
         else if (size(distribution_mode) == n) then
            me%distribution_mode = distribution_mode
         else
            error stop 'Error: distribution_mode must be scalar or size n.'
         end if
      end if
      if (present(dist_std_dev)) then
         if (size(dist_std_dev) == 1) then
            me%dist_std_dev = abs(dist_std_dev(1))
         else if (size(dist_std_dev) == n) then
            me%dist_std_dev = abs(dist_std_dev)
         else
            error stop 'Error: dist_std_dev must be scalar or size n.'
         end if
      end if
      if (present(dist_scale)) then
         if (size(dist_scale) == 1) then
            me%dist_scale = abs(dist_scale(1))
         else if (size(dist_scale) == n) then
            me%dist_scale = abs(dist_scale)
         else
            error stop 'Error: dist_scale must be scalar or size n.'
         end if
      end if
      if (present(dist_shape)) then
         if (size(dist_shape) == 1) then
            me%dist_shape = abs(dist_shape(1))
         else if (size(dist_shape) == n) then
            me%dist_shape = abs(dist_shape)
         else
            error stop 'Error: dist_shape must be scalar or size n.'
         end if
      end if

      ! validate distribution modes:
      if (any(me%distribution_mode < 1 .or. me%distribution_mode > n_distribution_modes)) then
         error stop 'Error: invalid distribution_mode.'
      end if

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

      class(simulated_annealing_type),intent(inout) :: me
      real(wp), dimension(me%n), intent(inout)  :: x       !! on input: the starting values for the variables of the
                                                           !! function to be optimized. [Will be replaced by final point]
      real(wp), intent(in)                      :: rt      !! the temperature reduction factor.  the value suggested by
                                                           !! corana et al. is .85. see goffe et al. for more advice.
      real(wp), intent(inout)                   :: t       !! on input, the initial temperature. see goffe et al. for advice.
                                                           !! on output, the final temperature. Note that if `t=0`, then
                                                           !! all downhill steps will be rejected
      real(wp), dimension(me%n), intent(inout)  :: vm      !! the step length vector. on input it should encompass the region of
                                                           !! interest given the starting value x.  for point x(i), the next
                                                           !! trial point is selected is from x(i) - vm(i)  to  x(i) + vm(i).
                                                           !! since vm is adjusted so that about half of all points are accepted,
                                                           !! the input value is not very important (i.e. is the value is off,
                                                           !! sa adjusts vm to the correct value).
      real(wp), dimension(me%n), intent(out)  :: xopt      !! the variables that optimize the function.
      real(wp), intent(out)                   :: fopt      !! the optimal value of the function.
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

      real(wp),dimension(:),allocatable :: xp           !! perturbed `x` vector (size `n`)
      real(wp),dimension(:),allocatable :: fstar        !! array of optimals found so far (size `neps`)
      real(wp),dimension(:),allocatable :: vm_original  !! original input value of `vm` (size `n`)
      real(wp)  :: f             !! function value
      real(wp)  :: fp            !! function value
      real(wp)  :: p             !! for metropolis criteria
      real(wp)  :: pp            !! for metropolis criteria
      real(wp)  :: ratio         !! ratio of `nacp` to `ns`
      integer   :: nup           !! number of uphill steps
      integer   :: ndown         !! number of accepted downhill steps
      integer   :: nrej          !! number of rejected downhill steps
      integer   :: nnew          !! new optimal for this temperature
      integer   :: i, j, m, k    !! counter
      integer   :: totmov        !! total moves
      integer   :: nacp          !! number of accepts steps in a complete `ns` cycle
      logical   :: quit          !! if the termination criteria was achieved
      integer   :: unit          !! unit for printing
      real(wp)  :: t_original    !! original input value of `t`
      logical   :: first         !! indicates first function eval
      logical   :: abort         !! indicates the known solution has been found
      real(wp)  :: t0            !! initial temperature (for non-geometric cooling schedules)
      integer   :: temp_iter     !! temperature iteration counter (for non-geometric cooling schedules)

      ! initialize:
      allocate(xp(me%n))
      allocate(fstar(me%neps))
      allocate(vm_original(me%n))
      unit        = me%iunit
      nfcnev      = 0
      ier         = 99
      xopt        = x
      fopt        = huge(1.0_wp)
      vm          = abs(vm)
      t_original  = t
      vm_original = vm
      abort       = .false.

      ! if the initial temperature is not positive, notify the user and
      ! return to the calling routine.
      if (t < 0.0_wp) then
         if (me%iprint>0) write(unit,'(A)') 'Error: The initial temperature must be non-negative.'
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
      if (me%iprint >= 1) then
         write(unit, '(A)') '  '
         call print_vector(unit,xp,me%n,'initial x')
         write(unit, '(A,1X,G25.18)')  '  initial f: ', me%func(f)
      end if

      reset: do k = 1, me%n_resets

         if (me%iprint >= 1) then
            write(unit,'(A)')       ''
            write(unit,'(A)')       '=========================================='
            write(unit,'(A,1X,I5)') ' main loop ', k
            write(unit,'(A)')       '=========================================='
            write(unit,'(A)')       ''
         end if

         ! restart the main loop
         fstar = huge(1.0_wp)
         t     = t_original
         vm    = vm_original
         nacc  = 0
         temp_iter = 0  ! reset temperature iteration counter
         nacp  = 0

         ! the first function eval for a new main cycle
         first = k > 1

         main : do
            ! start the main loop. note that it terminates if (i) the algorithm
            ! succesfully optimizes the function or (ii) there are too many
            ! function evaluations (more than maxevl).
            nup = 0
            nrej = 0
            nnew = 0
            ndown = 0

            nt_loop : do m = 1, me%nt
               ns_loop : do j = 1, me%ns

                  ! evaluate the function with the trial point xp and return as fp.
                  call me%perturb_and_evaluate(x,vm,xp,fp,nfcnev,ier,first=first)
                  first = .false.
                  select case (ier)
                   case(1,4)
                     x    = xopt
                     fopt = me%func(fopt)
                     return
                  end select

                  fp = me%func(fp)
                  if (me%iprint >= 3) then
                     write(unit,'(A)') ' '
                     call print_vector(unit,x,me%n,'current x')
                     write(unit,'(A,G25.18)') ' current f: ', me%func(f)
                     call print_vector(unit,xp,me%n,'trial x')
                     write(unit,'(A,G25.18)') ' resulting f: ', me%func(fp)
                  end if

                  if (fp >= f) then
                     ! function value increased: accept the new point

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

                  else
                     ! if the point is lower, use the metropolis criteria to decide on
                     ! acceptance or rejection.
                     ! NOTE: if t=0, all downhill steps are rejected (monotonic)

                     if (t/=0.0_wp) then
                        p = exprep((fp - f)/t)
                        pp = uniform_random_number()
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
                           cycle  ! accepted
                        end if
                     end if

                     ! rejected:
                     nrej = nrej + 1
                     if (me%iprint >= 3) then
                        if (me%maximize) then
                           write(unit,'(A)') '  lower point rejected'
                        else
                           write(unit,'(A)') '  higher point rejected'
                        end if
                     end if

                  end if

               end do ns_loop ! j - ns loop

               select case (me%step_mode)
                case(1)
                  ! adjust vm so that approximately half of all evaluations are accepted.
                  ratio = real(nacp,wp) /real(me%ns,wp)
                  do i = 1, me%n
                     if (ratio > 0.6_wp) then
                        vm(i) = vm(i)*(1.0_wp + me%c(i)*(ratio - 0.6_wp)/0.4_wp)
                     else if (ratio < 0.4_wp) then
                        vm(i) = vm(i)/(1.0_wp + me%c(i)*((0.4_wp - ratio)/0.4_wp))
                     end if
                     if (vm(i) > (me%ub(i)-me%lb(i))) then
                        vm(i) = me%ub(i) - me%lb(i)
                     end if
                  end do
                  if (me%iprint >= 2) then
                     write(unit,'(/A)') '---------------------------------------------------'
                     write(unit,'(A)')  ' intermediate results after step length adjustment '
                     write(unit,'(A/)') '---------------------------------------------------'
                     call print_vector(unit, vm,   me%n, 'new step length (vm)')
                     call print_vector(unit, xopt, me%n, 'current optimal x')
                     call print_vector(unit, x,    me%n, 'current x')
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

            end do nt_loop ! m - nt loop

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
               call print_vector(unit,xopt, me%n, 'current optimal x')
               call print_vector(unit,vm,   me%n, 'step length (vm)')
               write(unit, '(A)') ' '
            end if

            ! check termination criteria.
            quit = .false.
            ! first, check `f` if the known value is available
            if (me%optimal_f_specified) then
               quit = abs(me%func(fopt) - me%optimal_f) <= abs(me%optimal_f_tol)
               if (quit) abort = .true.
            end if
            ! check the convergence tolerances:
            if (.not. quit) then
               fstar(1) = f
               quit = ((fopt-f)<=me%eps) .and. (all(abs(f-fstar)<=me%eps))
            end if

            ! terminate sa if appropriate.
            if (quit) then
               x = xopt
               ier = 0
               fopt = me%func(fopt)
               if (me%iprint >= 1) then
                  write(unit,'(/A)') '----------------------------------------------'
                  write(unit,'(A)')  '  sa achieved termination criteria. ier = 0.  '
                  write(unit,'(A/)') '----------------------------------------------'
               end if
               if (abort) then
                  exit reset ! the solution has been found, so can ignore
                  ! the reset loop if that is being used
               else
                  exit
               end if
            end if

            !  if termination criteria is not met, prepare for another loop.
            temp_iter = temp_iter + 1

            ! apply the selected cooling schedule
            select case (me%cooling_schedule)
             case (1)
               ! geometric (classical): T(k+1) = rt * T(k)
               t = abs(rt)*t
             case (2)
               ! fast annealing (Cauchy): T(k) = T0 / (1 + k)
               t = t_original / real(1 + temp_iter, wp)
             case (3)
               ! huang: T(k) = T0 / (1 + c*k)^n
               t = t_original / (1.0_wp + me%cooling_param * real(temp_iter, wp))**me%n
             case (4)
               ! boltzmann: T(k) = T0 / log(1 + k + e)
               t = t_original / log(1.0_wp + real(temp_iter, wp) + exp(1.0_wp))
             case (5)
               ! logarithmic: T(k) = T0 / (1 + c*log(1+k))
               t = t_original / (1.0_wp + me%cooling_param * log(1.0_wp + real(temp_iter, wp)))
             case default
               ! fallback to geometric
               t = abs(rt)*t
            end select

            do i = me%neps, 2, -1
               fstar(i) = fstar(i-1)
            end do
            f = fopt
            x = xopt

         end do main

      end do reset

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

      class(simulated_annealing_type),intent(inout) :: me
      real(wp),dimension(:),intent(in)  :: x       !! input optimization variable vector to perturb
      real(wp),dimension(:),intent(in)  :: vm      !! step length vector
      real(wp),dimension(:),intent(out) :: xp      !! the perturbed `x` value
      real(wp),intent(out)              :: fp      !! the value of the user function at `xp`
      integer,intent(inout)             :: nfcnev  !! total number of function evaluations
      integer,intent(inout)             :: ier     !! status output code
      logical,intent(in),optional       :: first   !! to use the input `x` the first time

      integer :: i          !! counter
      integer :: istat      !! user function status code
      logical :: first_try  !! local copy of `first`
      integer :: n_inputs   !! number of inputs to send to the user function for parallel evaluation
      logical :: reallocate !! whether to reallocate `xp_mat` for parallel evaluation
      real(wp),dimension(:,:),allocatable :: xp_mat !! array of `xp` vectors for parallel evaluation

      if (present(first)) then
         first_try = first
      else
         first_try = .false.
      end if

      do  ! we must return at least one that doesn't fail (or we exceed the max number of function evaluations)

         if (me%parallel_mode) then
            ! parallel mode: generate and send multiple x vectors
            call me%n_inputs_to_send(n_inputs)  !! get the inputs to send to the user function for parallel evaluation
            n_inputs = max(1, n_inputs)  !! there must be at least one
            ! do we need to reallocate xp_mat? (if n_inputs has changed):
            reallocate = .true.
            if (allocated(xp_mat)) then
               if (size(xp_mat,1) == me%n .and. size(xp_mat,2) == n_inputs) then
                  reallocate = .false.
               end if
            end if
            if (reallocate) then
               if (allocated(xp_mat)) deallocate(xp_mat)
               allocate(xp_mat(me%n, n_inputs))
            end if
            do i = 1, n_inputs
               xp_mat(:,i) = get_xp()  ! get a perturbed `x` value for each input
            end do
            call me%fcn_parallel_input(xp_mat)  !! send the perturbed `x` values to the user function for parallel evaluation
            call me%fcn_parallel_output(xp, fp, istat)  !! get the function values and status code back from the user function for parallel evaluation
         else
            ! serial mode: generate and send a single x vector
            xp = get_xp()  ! single funciton evaluation
            call me%fcn(xp, fp, istat)  ! evaluate the function with the trial point xp and return as fp.
         end if
         nfcnev = nfcnev + 1  ! function eval counter (note: for parallel runs,
                              ! only count one evaluation per returned values.
                              ! others can still be running)

         ! if too many function evaluations occur, terminate the algorithm.
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

   contains

      function get_xp() result(xp)
         !! get a perturbed `x` value
         real(wp),dimension(me%n) :: xp     !! the perturbed `x` value
         real(wp)                 :: lower  !! lower bound to use for random interval
         real(wp)                 :: upper  !! upper bound to use for random interval
         integer                  :: i      !! counter

         if (first_try) then
            if (me%use_initial_guess) then
               ! use the initial guess
               ! [note if this evauation fails, the subsequent ones
               !  are perturbed from this one]
               xp = x
               first_try = .false.
            else
               ! a random point in the bounds:
               ! [if it fails, a new random one is tried next time]
               do i = 1, me%n
                  xp(i) = me%perturb_variable(i, x(i), me%distribution_mode(i), &
                     me%lb(i), me%ub(i))
               end do
            end if
         else
            ! perturb all of them using per-variable distributions:
            do i = 1, me%n
               lower = max( me%lb(i), x(i) - vm(i) )
               upper = min( me%ub(i), x(i) + vm(i) )
               xp(i) = me%perturb_variable(i, x(i), me%distribution_mode(i), &
                  lower, upper)
            end do
         end if
      end function get_xp

   end subroutine perturb_and_evaluate
!********************************************************************************

!********************************************************************************
!>
!  if the function is to be minimized, switch the sign of the function.
!  note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.

   pure function func(me,f)

      class(simulated_annealing_type),intent(in) :: me
      real(wp),intent(in) :: f
      real(wp) :: func

      if (me%maximize) then
         func = f
      else
         func = -f
      end if

   end function func
!********************************************************************************

!********************************************************************************
!>
!  Perturb a single variable using its assigned distribution and parameters.

   function perturb_variable(me, ivar, x, mode, lower, upper) result(r)

      class(simulated_annealing_type),intent(inout) :: me
      integer,intent(in)  :: ivar   !! variable index
      real(wp),intent(in) :: x      !! variable value to perturb
      integer,intent(in)  :: mode   !! perturbation distribution mode (see `distribution_mode` for details)
      real(wp),intent(in) :: lower  !! lower bound
      real(wp),intent(in) :: upper  !! upper bound
      real(wp)            :: r      !! perturbed value

      integer :: i !! counter
      integer,parameter :: max_tries = 1000 !! max tries for rejection sampling

      ! select distribution based on the variable's distribution_mode:
      select case (mode)

       case(sa_mode_uniform)  ! uniform
         r = uniform(lower, upper)

       case(sa_mode_normal)  ! normal (truncated)
         ! center the distribution on the current value of the variable
         r = truncated_normal(x, me%dist_std_dev(ivar), lower, upper)

       case(sa_mode_cauchy)  ! cauchy
         ! center the distribution on the current value of the variable
         ! rejection sampling to ensure within bounds
         do i = 1, max_tries
            r = cauchy(x, me%dist_scale(ivar))
            if (r >= lower .and. r <= upper) return
         end do
         ! fallback to uniform if rejection sampling fails
         r = uniform(lower, upper)

       case(sa_mode_triangular)  ! triangular
         ! center the peak at the current value of the variable
         ! compute normalized position of x in [lower, upper]
         r = lower + (upper - lower) * triangular_dist((x - lower) / (upper - lower))

       case(sa_mode_bipareto)  ! bipareto
         ! center the distribution on the current value of the variable
         r = bipareto(x, me%dist_scale(ivar), &
                      me%dist_shape(ivar), lower, upper)

       case default
         error stop 'Error: invalid distribution_mode in perturb_variable'
      end select

   end function perturb_variable
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

      logical,dimension(2) :: flags
      type(ieee_flag_type),parameter,dimension(2) :: out_of_range = &
         [ieee_overflow,ieee_underflow]

      real(wp), intent(in) :: x
      real(wp) :: f

      call ieee_set_halting_mode(out_of_range,.false.)

      f = exp(x)

      call ieee_get_flag(out_of_range,flags)
      if (any(flags)) then
         call ieee_set_flag(out_of_range,.false.)
         if (flags(1)) then
            f = huge(1.0_wp)
         else
            f = 0.0_wp
         end if
      end if

   end function exprep
!********************************************************************************

!********************************************************************************
!>
!  Initialize the intrinsic random number generator.
!
!### Author
!  * Jacob Williams, 8/30/2019

   subroutine rand_init(seed1,seed2)

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
!  Get a new uniform random number from [0,1].
!
!### Author
!  * Jacob Williams, 8/30/2019

   function uniform_random_number() result(f)

      real(wp) :: f

      call random_number(f)

   end function uniform_random_number
!********************************************************************************

!********************************************************************************
!>
!  Uniform random number on the interval `[xl,xu]`.

   function uniform(xl,xu)

      real(wp),intent(in) :: xl !! lower bound
      real(wp),intent(in) :: xu !! upper bound
      real(wp) :: uniform

      uniform = xl + (xu-xl)*uniform_random_number()

   end function uniform
!********************************************************************************

!********************************************************************************
!>
!  Normal (Gaussian) random number with specified mean and standard deviation.
!  Uses the Box-Muller transform.

   function normal(mean, std_dev)

      real(wp),intent(in) :: mean    !! mean of the distribution
      real(wp),intent(in) :: std_dev !! standard deviation
      real(wp) :: normal

      real(wp) :: u1, u2

      u1 = uniform_random_number()
      u2 = uniform_random_number()

      normal = mean + std_dev * sqrt(-2.0_wp * log(u1)) * cos(twopi * u2)

   end function normal
!********************************************************************************

!********************************************************************************
!>
!  Cauchy random number with specified location and scale parameters.
!  The Cauchy distribution has heavier tails than the normal distribution,
!  which can be useful for occasional large jumps in simulated annealing.

   function cauchy(location, scale)

      real(wp),intent(in) :: location !! location parameter (median)
      real(wp),intent(in) :: scale    !! scale parameter
      real(wp) :: cauchy

      real(wp) :: u

      u = uniform_random_number()
      cauchy = location + scale * tan(pi * (u - 0.5_wp))

   end function cauchy
!********************************************************************************

!********************************************************************************
!>
!  Truncated normal distribution within bounds [xl, xu].
!  Uses rejection sampling to ensure the value stays within bounds.
!
!### Note
!  For efficiency, if bounds are more than a few standard deviations
!  from the mean, consider using uniform distribution instead.

   function truncated_normal(mean, std_dev, xl, xu)

      real(wp),intent(in) :: mean    !! mean of the underlying normal distribution
      real(wp),intent(in) :: std_dev !! standard deviation
      real(wp),intent(in) :: xl      !! lower bound
      real(wp),intent(in) :: xu      !! upper bound
      real(wp) :: truncated_normal

      integer,parameter :: max_tries = 1000
      integer :: i

      ! rejection sampling
      do i = 1, max_tries
         truncated_normal = normal(mean, std_dev)
         if (truncated_normal >= xl .and. truncated_normal <= xu) return
      end do

      ! fallback to uniform if rejection sampling fails
      truncated_normal = uniform(xl, xu)

   end function truncated_normal
!********************************************************************************

!********************************************************************************
!>
!  Triangular distribution on [0,1] with specified mode.
!  Uses inverse transform sampling - very efficient, no rejection needed.
!
!### Note
!  The mode parameter should be in [0,1]. Common choices:
!  - mode = 0.5: symmetric triangular
!  - mode < 0.5: left-skewed
!  - mode > 0.5: right-skewed

   function triangular_dist(mode)

      real(wp),intent(in) :: mode !! mode of the distribution (should be in [0,1])
      real(wp) :: triangular_dist

      real(wp) :: u, mode_clipped

      ! ensure mode is in [0,1]
      mode_clipped = max(0.0_wp, min(1.0_wp, mode))

      u = uniform_random_number()

      if (u < mode_clipped) then
         ! left side of triangle
         triangular_dist = sqrt(u * mode_clipped)
      else
         ! right side of triangle
         triangular_dist = 1.0_wp - sqrt((1.0_wp - u) * (1.0_wp - mode_clipped))
      end if

   end function triangular_dist
!********************************************************************************

!********************************************************************************
!>
!  Bi-polar (two-sided) Pareto distribution.
!  Generates a symmetric heavy-tailed distribution centered at a location.
!
!### Note
!  Unlike the standard Pareto which is one-sided (heavy tail on right only),
!  the bi-polar Pareto is symmetric and allows large jumps in both directions.
!  This is particularly useful for optimization where perturbations should be
!  unbiased with respect to direction.
!
!  The distribution:
!  - Randomly chooses a direction (positive or negative with equal probability)
!  - Generates a Pareto-distributed magnitude
!  - Applies the signed magnitude to the center location
!  - Uses rejection sampling to ensure the result stays within [xl, xu]

   function bipareto(center, scale, shape, xl, xu)

      real(wp),intent(in) :: center !! center location of the distribution
      real(wp),intent(in) :: scale  !! scale parameter for the Pareto distribution
      real(wp),intent(in) :: shape  !! shape parameter (smaller values = heavier tails)
      real(wp),intent(in) :: xl     !! lower bound
      real(wp),intent(in) :: xu     !! upper bound
      real(wp) :: bipareto

      integer,parameter :: max_tries = 1000
      integer :: i
      real(wp) :: magnitude, sign_val, u

      ! rejection sampling to ensure within bounds
      do i = 1, max_tries
         ! randomly choose direction: -1 or +1
         if (uniform_random_number() < 0.5_wp) then
            sign_val = -1.0_wp
         else
            sign_val = 1.0_wp
         end if

         ! generate Pareto-distributed magnitude (inlined)
         u = uniform_random_number()
         magnitude = (scale / u**(1.0_wp / shape)) - scale

         ! apply signed magnitude to center
         bipareto = center + sign_val * magnitude

         ! check if within bounds
         if (bipareto >= xl .and. bipareto <= xu) return
      end do

      ! fallback to uniform if rejection sampling fails
      bipareto = uniform(xl, xu)

   end function bipareto
!********************************************************************************

!********************************************************************************
!>
!  this subroutine prints the double precision vector named vector.
!  elements 1 thru ncols will be printed. name is a character variable
!  that describes vector. note that if name is given in the call to
!  print_vector, it must be enclosed in quotes. if there are more than 10
!  elements in vector, 10 elements will be printed on each line.

   subroutine print_vector(iunit, vector, ncols, name)

      integer, intent(in)                     :: iunit
      integer, intent(in)                     :: ncols
      real(wp), dimension(ncols), intent(in)  :: vector
      character(len=*), intent(in)            :: name

      integer :: i, lines, ll

      write(iunit,'(/,25X,A)') trim(name)

      if (ncols > 10) then
         lines = int(ncols/10.0_wp)
         do i = 1, lines
            ll = 10*(i - 1)
            write(iunit,'(10(G12.5,1X))') vector(1+ll:10+ll)
         end do
         write(iunit,'(10(G12.5,1X))') vector(11+ll:ncols)
      else
         write(iunit,'(10(G12.5,1X))') vector(1:ncols)
      end if

   end subroutine print_vector
!********************************************************************************

!********************************************************************************
end module simulated_annealing_module
!********************************************************************************
