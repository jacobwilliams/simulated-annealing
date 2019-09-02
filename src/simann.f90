!  This file is an example of the Corana et al. simulated annealing algorithm
!  for multimodal and robust optimization as implemented and modified by
!  Goffe, Ferrier and Rogers.  Counting the above line
!  ABSTRACT as 1, the routine itself (SA), with its supplementary
!  routines, is on lines 232-990. A multimodal example from Judge et al.
!  (FCN) is on lines 150-231. The rest of this file (lines 1-149) is a
!  driver routine with values appropriate for the Judge example.  Thus, this
!  example is ready to run.

!  To understand the algorithm, the documentation for SA on lines 236-
!  484 should be read along with the parts of the paper that describe
!  simulated annealing. Then the following lines will then aid the user
!  in becomming proficient with this implementation of simulated annealing.

!  Learning to use SA:
!      Use the sample function from Judge with the following suggestions
!  to get a feel for how SA works. When you've done this, you should be
!  ready to use it on most any function with a fair amount of expertise.
!    1. Run the program as is to make sure it runs okay. Take a look at
!       the intermediate output and see how it optimizes as temperature
!       (T) falls.  Notice how the optimal point is reached and how
!       falling T reduces VM.
!    2. Look through the documentation to SA so the following makes a
!       bit of sense. In line with the paper, it shouldn't be that hard
!       to figure out. The core of the algorithm is described on pp. 68-70
!       and on pp. 94-95. Also see Corana et al. pp. 264-9.
!    3. To see how it selects points and makes decisions about uphill and
!       downhill moves, set IPRINT = 3 (very detailed intermediate output)
!       and MAXEVL = 100 (only 100 function evaluations to limit output).
!    4. To see the importance of different temperatures, try starting
!       with a very low one (say T = 10E-5).  You'll see (i) it never
!       escapes from the local optima (in annealing terminology, it
!       quenches) & (ii) the step length (VM) will be quite small.  This
!       is a key part of the algorithm: as temperature (T) falls, step
!       length does too. In a minor point here, note how VM is quickly
!       reset from its initial value. Thus, the input VM is not very
!       important.  This is all the more reason to examine VM once the
!       algorithm is underway.
!    5. To see the effect of different parameters and their effect on
!       the speed of the algorithm, try RT = .95 & RT = .1.  Notice the
!       vastly different speed for optimization. Also try NT = 20.  Note
!       that this sample function is quite easy to optimize, so it will
!       tolerate big changes in these parameters.  RT and NT are the
!       parameters one should adjust to modify the runtime of the
!       algorithm and its robustness.
!    6. Try constraining the algorithm with either LB or UB.


MODULE simulated_anneal

! ABSTRACT:
!   Simulated annealing is a global optimization method that distinguishes
!   between different local optima.  Starting from an initial point, the
!   algorithm takes a step and the function is evaluated. When minimizing a
!   function, any downhill step is accepted and the process repeats from this
!   new point. An uphill step may be accepted. Thus, it can escape from local
!   optima. This uphill decision is made by the Metropolis criteria. As the
!   optimization process proceeds, the length of the steps decline and the
!   algorithm closes in on the global optimum. Since the algorithm makes very
!   few assumptions regarding the function to be optimized, it is quite
!   robust with respect to non-quadratic surfaces. The degree of robustness
!   can be adjusted by the user. In fact, simulated annealing can be used as
!   a local optimizer for difficult functions.

!   This implementation of simulated annealing was used in "Global Optimizatio
!   of Statistical Functions with Simulated Annealing," Goffe, Ferrier and
!   Rogers, Journal of Econometrics, vol. 60, no. 1/2, Jan./Feb. 1994, pp.
!   65-100. Briefly, we found it competitive, if not superior, to multiple
!   restarts of conventional optimization routines for difficult optimization
!   problems.

!   For more information on this routine, contact its author:
!   Bill Goffe, bgoffe@whale.st.usm.edu

! This version in Fortran 90 has been prepared by Alan Miller.
! It is compatible with Lahey's ELF90 compiler.
! N.B. The 3 last arguments have been removed from subroutine sa.   these
!      were work arrays and are now internal to the routine.
! e-mail:  amiller@bigpond.net.au
! URL   :  http://users.bigpond.net.au/amiller

! Latest revision of Fortran 90 version - 3 August 1997

IMPLICIT NONE

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

! The following variables were in COMMON /raset1/
REAL, SAVE    :: u(97), cc, cd, cm
INTEGER, SAVE :: i97, j97


CONTAINS

SUBROUTINE sa(n, x, max, rt, eps, ns, nt, neps, maxevl, lb, ub, c, iprint,  &
              iseed1, iseed2, t, vm, xopt, fopt, nacc, nfcnev, nobds, ier)

!  Version: 3.2
!  Date: 1/22/94.
!  Differences compared to Version 2.0:
!     1. If a trial is out of bounds, a point is randomly selected
!        from LB(i) to UB(i). Unlike in version 2.0, this trial is
!        evaluated and is counted in acceptances and rejections.
!        All corresponding documentation was changed as well.
!  Differences compared to Version 3.0:
!     1. If VM(i) > (UB(i) - LB(i)), VM is set to UB(i) - LB(i).
!        The idea is that if T is high relative to LB & UB, most
!        points will be accepted, causing VM to rise. But, in this
!        situation, VM has little meaning; particularly if VM is
!        larger than the acceptable region. Setting VM to this size
!        still allows all parts of the allowable region to be selected.
!  Differences compared to Version 3.1:
!     1. Test made to see if the initial temperature is positive.
!     2. WRITE statements prettied up.
!     3. References to paper updated.

!  Synopsis:
!  This routine implements the continuous simulated annealing global
!  optimization algorithm described in Corana et al.'s article "Minimizing
!  Multimodal Functions of Continuous Variables with the "Simulated Annealing"
!  Algorithm" in the September 1987 (vol. 13, no. 3, pp. 262-280) issue of
!  the ACM Transactions on Mathematical Software.

!  A very quick (perhaps too quick) overview of SA:
!     SA tries to find the global optimum of an N dimensional function.
!  It moves both up and downhill and as the optimization process
!  proceeds, it focuses on the most promising area.
!     To start, it randomly chooses a trial point within the step length
!  VM (a vector of length N) of the user selected starting point. The
!  function is evaluated at this trial point and its value is compared
!  to its value at the initial point.
!     In a maximization problem, all uphill moves are accepted and the
!  algorithm continues from that trial point. Downhill moves may be
!  accepted; the decision is made by the Metropolis criteria. It uses T
!  (temperature) and the size of the downhill move in a probabilistic
!  manner. The smaller T and the size of the downhill move are, the more
!  likely that move will be accepted. If the trial is accepted, the
!  algorithm moves on from that point. If it is rejected, another point
!  is chosen instead for a trial evaluation.
!     Each element of VM periodically adjusted so that half of all
!  function evaluations in that direction are accepted.
!     A fall in T is imposed upon the system with the RT variable by
!  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
!  downhill moves are less likely to be accepted and the percentage of
!  rejections rise. Given the scheme for the selection for VM, VM falls.
!  Thus, as T declines, VM falls and SA focuses upon the most promising
!  area for optimization.

!  The importance of the parameter T:
!     The parameter T is crucial in using SA successfully. It influences
!  VM, the step length over which the algorithm searches for optima. For
!  a small intial T, the step length may be too small; thus not enough
!  of the function might be evaluated to find the global optima. The user
!  should carefully examine VM in the intermediate output (set IPRINT =
!  1) to make sure that VM is appropriate. The relationship between the
!  initial temperature and the resulting step length is function
!  dependent.
!     To determine the starting temperature that is consistent with
!  optimizing a function, it is worthwhile to run a trial run first. Set
!  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
!  rises as well. Then select the T that produces a large enough VM.

!  For modifications to the algorithm and many details on its use,
!  (particularly for econometric applications) see Goffe, Ferrier
!  and Rogers, "Global Optimization of Statistical Functions with
!  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2,
!  Jan./Feb. 1994, pp. 65-100.
!  For more information, contact
!              Bill Goffe
!              Department of Economics and International Business
!              University of Southern Mississippi
!              Hattiesburg, MS  39506-5072
!              (601) 266-4484 (office)
!              (601) 266-4920 (fax)
!              bgoffe@whale.st.usm.edu (Internet)

!  As far as possible, the parameters here have the same name as in
!  the description of the algorithm on pp. 266-8 of Corana et al.

!  In this description, SP is single precision, DP is double precision,
!  INT is integer, L is logical and (N) denotes an array of length n.
!  Thus, DP(N) denotes a double precision array of length n.

!  Input Parameters:
!    Note: The suggested values generally come from Corana et al. To
!          drastically reduce runtime, see Goffe et al., pp. 90-1 for
!          suggestions on choosing the appropriate RT and NT.
!    N - Number of variables in the function to be optimized. (INT)
!    X - The starting values for the variables of the function to be
!        optimized. (DP(N))
!    MAX - Denotes whether the function should be maximized or minimized.
!          A true value denotes maximization while a false value denotes
!          minimization.  Intermediate output (see IPRINT) takes this into
!          account. (L)
!    RT - The temperature reduction factor.  The value suggested by
!         Corana et al. is .85. See Goffe et al. for more advice. (DP)
!    EPS - Error tolerance for termination. If the final function
!          values from the last neps temperatures differ from the
!          corresponding value at the current temperature by less than
!          EPS and the final function value at the current temperature
!          differs from the current optimal function value by less than
!          EPS, execution terminates and IER = 0 is returned. (EP)
!    NS - Number of cycles.  After NS*N function evaluations, each element of
!         VM is adjusted so that approximately half of all function evaluations
!         are accepted.  The suggested value is 20. (INT)
!    NT - Number of iterations before temperature reduction. After
!         NT*NS*N function evaluations, temperature (T) is changed
!         by the factor RT.  Value suggested by Corana et al. is
!         MAX(100, 5*N).  See Goffe et al. for further advice. (INT)
!    NEPS - Number of final function values used to decide upon termi-
!           nation.  See EPS.  Suggested value is 4. (INT)
!    MAXEVL - The maximum number of function evaluations.  If it is
!             exceeded, IER = 1. (INT)
!    LB - The lower bound for the allowable solution variables. (DP(N))
!    UB - The upper bound for the allowable solution variables. (DP(N))
!         If the algorithm chooses X(I) .LT. LB(I) or X(I) .GT. UB(I),
!         I = 1, N, a point is from inside is randomly selected. This
!         This focuses the algorithm on the region inside UB and LB.
!         Unless the user wishes to concentrate the search to a particular
!         region, UB and LB should be set to very large positive
!         and negative values, respectively.  Note that the starting
!         vector X should be inside this region.  Also note that LB and
!         UB are fixed in position, while VM is centered on the last
!         accepted trial set of variables that optimizes the function.
!    C - Vector that controls the step length adjustment.  The suggested
!        value for all elements is 2.0. (DP(N))
!    IPRINT - controls printing inside SA. (INT)
!             Values: 0 - Nothing printed.
!                     1 - Function value for the starting value and
!                         summary results before each temperature
!                         reduction. This includes the optimal
!                         function value found so far, the total
!                         number of moves (broken up into uphill,
!                         downhill, accepted and rejected), the
!                         number of out of bounds trials, the
!                         number of new optima found at this
!                         temperature, the current optimal X and
!                         the step length VM. Note that there are
!                         N*NS*NT function evalutations before each
!                         temperature reduction. Finally, notice is
!                         is also given upon achieveing the termination
!                         criteria.
!                     2 - Each new step length (VM), the current optimal
!                         X (XOPT) and the current trial X (X). This
!                         gives the user some idea about how far X
!                         strays from XOPT as well as how VM is adapting
!                         to the function.
!                     3 - Each function evaluation, its acceptance or
!                         rejection and new optima. For many problems,
!                         this option will likely require a small tree
!                         if hard copy is used. This option is best
!                         used to learn about the algorithm. A small
!                         value for MAXEVL is thus recommended when
!                         using IPRINT = 3.
!             Suggested value: 1
!             Note: For a given value of IPRINT, the lower valued
!                   options (other than 0) are utilized.
!    ISEED1 - The first seed for the random number generator RANMAR.
!             0 <= ISEED1 <= 31328. (INT)
!    ISEED2 - The second seed for the random number generator RANMAR.
!             0 <= ISEED2 <= 30081. Different values for ISEED1
!             and ISEED2 will lead to an entirely different sequence
!             of trial points and decisions on downhill moves (when
!             maximizing).  See Goffe et al. on how this can be used
!             to test the results of SA. (INT)

!  Input/Output Parameters:
!    T - On input, the initial temperature. See Goffe et al. for advice.
!        On output, the final temperature. (DP)
!    VM - The step length vector. On input it should encompass the region of
!         interest given the starting value X.  For point X(I), the next
!         trial point is selected is from X(I) - VM(I)  to  X(I) + VM(I).
!         Since VM is adjusted so that about half of all points are accepted,
!         the input value is not very important (i.e. is the value is off,
!         SA adjusts VM to the correct value). (DP(N))

!  Output Parameters:
!    XOPT - The variables that optimize the function. (DP(N))
!    FOPT - The optimal value of the function. (DP)
!    NACC - The number of accepted function evaluations. (INT)
!    NFCNEV - The total number of function evaluations. In a minor
!             point, note that the first evaluation is not used in the
!             core of the algorithm; it simply initializes the
!             algorithm. (INT).
!    NOBDS - The total number of trial function evaluations that
!            would have been out of bounds of LB and UB. Note that
!            a trial point is randomly selected between LB and UB. (INT)
!    IER - The error return number. (INT)
!          Values: 0 - Normal return; termination criteria achieved.
!                  1 - Number of function evaluations (NFCNEV) is
!                      greater than the maximum number (MAXEVL).
!                  2 - The starting value (X) is not inside the
!                      bounds (LB and UB).
!                  3 - The initial temperature is not positive.
!                  99 - Should not be seen; only used internally.

!  Work arrays that must be dimensioned in the calling routine:
!       RWK1 (DP(NEPS))  (FSTAR in SA)
!       RWK2 (DP(N))     (XP    "  " )
!       IWK  (INT(N))    (NACP  "  " )
!  N.B. In the Fortran 90 version, these are automatic arrays.

!  Required Functions (included):
!    EXPREP - Replaces the function EXP to avoid under- and overflows.
!             It may have to be modified for non IBM-type main-
!             frames. (DP)
!    RMARIN - Initializes the random number generator RANMAR.
!    RANMAR - The actual random number generator. Note that
!             RMARIN must run first (SA does this). It produces uniform
!             random numbers on [0,1]. These routines are from
!             Usenet's comp.lang.fortran. For a reference, see
!             "Toward a Universal Random Number Generator"
!             by George Marsaglia and Arif Zaman, Florida State
!             University Report: FSU-SCRI-87-50 (1987).
!             It was later modified by F. James and published in
!             "A Review of Pseudo-random Number Generators." For
!             further information, contact stuart@ads.com. These
!             routines are designed to be portable on any machine
!             with a 24-bit or more mantissa. I have found it produces
!             identical results on a IBM 3081 and a Cray Y-MP.

!  Required Subroutines (included):
!    PRTVEC - Prints vectors.
!    PRT1 ... PRT10 - Prints intermediate output.
!    FCN - Function to be optimized. The form is
!            SUBROUTINE FCN(N, X, F)
!            IMPLICIT NONE
!            INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
!            INTEGER, INTENT(IN)    :: N
!            REAL (dp), INTENT(IN)  :: X(N)
!            REAL (dp), INTENT(OUT) :: F
!            ...
!            function code with F = F(X)
!            ...
!            RETURN
!            END
!  Note: This is the same form used in the multivariable
!        minimization algorithms in the IMSL edition 10 library.

!  Machine Specific Features:
!    1. EXPREP may have to be modified if used on non-IBM type main-
!       frames. Watch for under- and overflows in EXPREP.
!    2. Some FORMAT statements use G25.18; this may be excessive for
!       some machines.
!    3. RMARIN and RANMAR are designed to be protable; they should not
!       cause any problems.

!  Type all external variables.
REAL (dp), INTENT(IN)     :: lb(:), ub(:), c(:), eps, rt
REAL (dp), INTENT(IN OUT) :: x(:), t, vm(:)
REAL (dp), INTENT(OUT)    :: xopt(:), fopt
INTEGER, INTENT(IN)       :: n, ns, nt, neps, maxevl, iprint, iseed1, iseed2
INTEGER, INTENT(OUT)      :: nacc, nfcnev, nobds, ier
LOGICAL, INTENT(IN)       :: max

!  Type all internal variables.
REAL (dp) :: f, fp, p, pp, ratio, xp(n), fstar(neps)
INTEGER   :: nup, ndown, nrej, nnew, lnobds, h, i, j, m, nacp(n)
LOGICAL   :: quit

INTERFACE
  SUBROUTINE fcn(n, theta, h)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)    :: n
    REAL (dp), INTENT(IN)  :: theta(:)
    REAL (dp), INTENT(OUT) :: h
  END SUBROUTINE fcn
END INTERFACE

!  Initialize the random number generator RANMAR.
CALL rmarin(iseed1, iseed2)

!  Set initial values.
nacc = 0
nobds = 0
nfcnev = 0
ier = 99

DO i = 1, n
  xopt(i) = x(i)
  nacp(i) = 0
END DO

fstar = 1.0D+20

!  If the initial temperature is not positive, notify the user and
!  return to the calling routine.
IF (t <= 0.0) THEN
  WRITE(*,'(/, "  THE INITIAL TEMPERATURE IS NOT POSITIVE. "/,  &
        &    "  reset the variable t. "/)')
  ier = 3
  RETURN
END IF

!  If the initial value is out of bounds, notify the user and return
!  to the calling routine.
DO i = 1, n
  IF ((x(i) > ub(i)) .OR. (x(i) < lb(i))) THEN
    CALL prt1()
    ier = 2
    RETURN
  END IF
END DO

!  Evaluate the function with input X and return value as F.
CALL fcn(n, x, f)

!  If the function is to be minimized, switch the sign of the function.
!  Note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.
IF(.NOT. max) f = -f
nfcnev = nfcnev + 1
fopt = f
fstar(1) = f
IF(iprint >= 1) CALL prt2(max, n, x, f)

!  Start the main loop. Note that it terminates if (i) the algorithm
!  succesfully optimizes the function or (ii) there are too many
!  function evaluations (more than MAXEVL).
100 nup = 0
nrej = 0
nnew = 0
ndown = 0
lnobds = 0

DO m = 1, nt
  DO j = 1, ns
    DO h = 1, n
      
!  Generate XP, the trial value of X. Note use of VM to choose XP.
      DO i = 1, n
        IF (i == h) THEN
          xp(i) = x(i) + (ranmar()*2. - 1.) * vm(i)
        ELSE
          xp(i) = x(i)
        END IF
        
!  If XP is out of bounds, select a point in bounds for the trial.
        IF((xp(i) < lb(i)) .OR. (xp(i) > ub(i))) THEN
          xp(i) = lb(i) + (ub(i) - lb(i))*ranmar()
          lnobds = lnobds + 1
          nobds = nobds + 1
          IF(iprint >= 3) CALL prt3(max, n, xp, x, f)
        END IF
      END DO
      
!  Evaluate the function with the trial point XP and return as FP.
      CALL fcn(n, xp, fp)
      IF(.NOT. max) fp = -fp
      nfcnev = nfcnev + 1
      IF(iprint >= 3) CALL prt4(max, n, xp, x, fp, f)
      
!  If too many function evaluations occur, terminate the algorithm.
      IF(nfcnev >= maxevl) THEN
        CALL prt5()
        IF (.NOT. max) fopt = -fopt
        ier = 1
        RETURN
      END IF
      
!  Accept the new point if the function value increases.
      IF(fp >= f) THEN
        IF(iprint >= 3) THEN
          WRITE(*,'("  POINT ACCEPTED")')
        END IF
        x(1:n) = xp(1:n)
        f = fp
        nacc = nacc + 1
        nacp(h) = nacp(h) + 1
        nup = nup + 1
        
!  If greater than any other point, record as new optimum.
        IF (fp > fopt) THEN
          IF(iprint >= 3) THEN
            WRITE(*,'("  NEW OPTIMUM")')
          END IF
          xopt(1:n) = xp(1:n)
          fopt = fp
          nnew = nnew + 1
        END IF
        
!  If the point is lower, use the Metropolis criteria to decide on
!  acceptance or rejection.
      ELSE
        p = exprep((fp - f)/t)
        pp = ranmar()
        IF (pp < p) THEN
          IF(iprint >= 3) CALL prt6(max)
          x(1:n) = xp(1:n)
          f = fp
          nacc = nacc + 1
          nacp(h) = nacp(h) + 1
          ndown = ndown + 1
        ELSE
          nrej = nrej + 1
          IF(iprint >= 3) CALL prt7(max)
        END IF
      END IF
      
    END DO
  END DO
  
!  Adjust VM so that approximately half of all evaluations are accepted.
  DO i = 1, n
    ratio = DBLE(nacp(i)) /DBLE(ns)
    IF (ratio > .6) THEN
      vm(i) = vm(i)*(1. + c(i)*(ratio - .6)/.4)
    ELSE IF (ratio < .4) THEN
      vm(i) = vm(i)/(1. + c(i)*((.4 - ratio)/.4))
    END IF
    IF (vm(i) > (ub(i)-lb(i))) THEN
      vm(i) = ub(i) - lb(i)
    END IF
  END DO
  
  IF(iprint >= 2) THEN
    CALL prt8(n, vm, xopt, x)
  END IF
  
  nacp(1:n) = 0
  
END DO

IF(iprint >= 1) THEN
  CALL prt9(max,n,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew)
END IF

!  Check termination criteria.
quit = .false.
fstar(1) = f
IF ((fopt - fstar(1)) <= eps) quit = .true.
DO i = 1, neps
  IF (ABS(f - fstar(i)) > eps) quit = .false.
END DO

!  Terminate SA if appropriate.
IF (quit) THEN
  x(1:n) = xopt(1:n)
  ier = 0
  IF (.NOT. max) fopt = -fopt
  IF(iprint >= 1) CALL prt10()
  RETURN
END IF

!  If termination criteria is not met, prepare for another loop.
t = rt*t
DO i = neps, 2, -1
  fstar(i) = fstar(i-1)
END DO
f = fopt
x(1:n) = xopt(1:n)

!  Loop again.
GO TO 100

END SUBROUTINE sa


FUNCTION exprep(rdum) RESULT(fn_val)
!  This function replaces exp to avoid under- and overflows and is
!  designed for IBM 370 type machines. It may be necessary to modify
!  it for other machines. Note that the maximum and minimum values of
!  EXPREP are such that they has no effect on the algorithm.

REAL (dp), INTENT(IN) :: rdum
REAL (dp)             :: fn_val

IF (rdum > 174._dp) THEN
  fn_val = 3.69D+75
ELSE IF (rdum < -180.+dp) THEN
  fn_val = 0.0_dp
ELSE
  fn_val = EXP(rdum)
END IF

RETURN
END FUNCTION exprep


SUBROUTINE rmarin(ij, kl)
!  This subroutine and the next function generate random numbers. See
!  the comments for SA for more information. The only changes from the
!  orginal code is that (1) the test to make sure that RMARIN runs first
!  was taken out since SA assures that this is done (this test didn't
!  compile under IBM's VS Fortran) and (2) typing ivec as integer was
!  taken out since ivec isn't used. With these exceptions, all following
!  lines are original.

! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081

INTEGER, INTENT(IN) :: ij, kl

INTEGER :: i, j, k, l, ii, jj, m
REAL    :: s, t

IF( ij < 0  .OR.  ij > 31328  .OR. kl < 0  .OR.  kl > 30081 ) THEN
  WRITE(*, '(A)') ' The first random number seed must have a value ',  &
               'between 0 AND 31328'
  WRITE(*, '(A)') ' The second seed must have a value between 0 and 30081'
  STOP
END IF
i = MOD(ij/177, 177) + 2
j = MOD(ij    , 177) + 2
k = MOD(kl/169, 178) + 1
l = MOD(kl,     169)
DO ii = 1, 97
  s = 0.0
  t = 0.5
  DO jj = 1, 24
    m = MOD(MOD(i*j, 179)*k, 179)
    i = j
    j = k
    k = m
    l = MOD(53*l+1, 169)
    IF (MOD(l*m, 64) >= 32) THEN
      s = s + t
    END IF
    t = 0.5 * t
  END DO
  u(ii) = s
END DO
cc = 362436.0 / 16777216.0
cd = 7654321.0 / 16777216.0
cm = 16777213.0 /16777216.0
i97 = 97
j97 = 33
RETURN
END SUBROUTINE rmarin


FUNCTION ranmar() RESULT(fn_val)
REAL :: fn_val

! Local variable
REAL :: uni

uni = u(i97) - u(j97)
IF( uni < 0.0 ) uni = uni + 1.0
u(i97) = uni
i97 = i97 - 1
IF(i97 == 0) i97 = 97
j97 = j97 - 1
IF(j97 == 0) j97 = 97
cc = cc - cd
IF( cc < 0.0 ) cc = cc + cm
uni = uni - cc
IF( uni < 0.0 ) uni = uni + 1.0
fn_val = uni

RETURN
END FUNCTION ranmar


SUBROUTINE prt1()
!  This subroutine prints intermediate output, as does PRT2 through
!  PRT10. Note that if SA is minimizing the function, the sign of the
!  function value and the directions (up/down) are reversed in all
!  output to correspond with the actual function optimization. This
!  correction is because SA was written to maximize functions and
!  it minimizes by maximizing the negative a function.

WRITE(*, '(/, "  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS "/,  &
      &     "  (lb AND ub). execution terminated without any"/,  &
      &     "  optimization. respecify x, ub OR lb so that  "/,  &
      &     "  lb(i) < x(i) < ub(i), i = 1, n. "/)')

RETURN
END SUBROUTINE prt1


SUBROUTINE prt2(max, n, x, f)

REAL (dp), INTENT(IN) ::  x(:), f
INTEGER, INTENT(IN)   ::  n
LOGICAL, INTENT(IN)   ::  max

WRITE(*, '("  ")')
CALL prtvec(x,n,'INITIAL X')
IF (max) THEN
  WRITE(*, '("  INITIAL F: ",/, G25.18)') f
ELSE
  WRITE(*, '("  INITIAL F: ",/, G25.18)') -f
END IF

RETURN
END SUBROUTINE prt2


SUBROUTINE prt3(max, n, xp, x, f)

REAL (dp), INTENT(IN) ::  xp(:), x(:), f
INTEGER, INTENT(IN)   ::  n
LOGICAL, INTENT(IN)   ::  max

WRITE(*, '("  ")')
CALL prtvec(x, n, 'CURRENT X')
IF (max) THEN
  WRITE(*, '("  CURRENT F: ", G25.18)') f
ELSE
  WRITE(*, '("  CURRENT F: ", G25.18)') -f
END IF
CALL prtvec(xp, n, 'TRIAL X')
WRITE(*, '("  POINT REJECTED SINCE OUT OF BOUNDS")')

RETURN
END SUBROUTINE prt3


SUBROUTINE prt4(max, n, xp, x, fp, f)

REAL (dp), INTENT(IN) ::  xp(:), x(:), fp, f
INTEGER, INTENT(IN)   ::  n
LOGICAL, INTENT(IN)   ::  max

WRITE(*,'("  ")')
CALL prtvec(x,n,'CURRENT X')
IF (max) THEN
  WRITE(*,'("  CURRENT F: ",G25.18)') f
  CALL prtvec(xp,n,'TRIAL X')
  WRITE(*,'("  RESULTING F: ",G25.18)') fp
ELSE
  WRITE(*,'("  CURRENT F: ",G25.18)') -f
  CALL prtvec(xp,n,'TRIAL X')
  WRITE(*,'("  RESULTING F: ",G25.18)') -fp
END IF

RETURN
END SUBROUTINE prt4


SUBROUTINE prt5()

WRITE(*, '(/, "  TOO MANY FUNCTION EVALUATIONS; CONSIDER "/,  &
      &     "  increasing maxevl OR eps, OR decreasing "/,  &
      &     "  nt OR rt. these results are likely TO be "/, "  poor.",/)')

RETURN
END SUBROUTINE prt5


SUBROUTINE prt6(max)
LOGICAL, INTENT(IN) ::  max

IF (max) THEN
  WRITE(*,'("  THOUGH LOWER, POINT ACCEPTED")')
ELSE
  WRITE(*,'("  THOUGH HIGHER, POINT ACCEPTED")')
END IF

RETURN
END SUBROUTINE prt6


SUBROUTINE prt7(max)
LOGICAL, INTENT(IN) :: max

IF (max) THEN
  WRITE(*,'("  LOWER POINT REJECTED")')
ELSE
  WRITE(*,'("  HIGHER POINT REJECTED")')
END IF

RETURN
END SUBROUTINE prt7


SUBROUTINE prt8(n, vm, xopt, x)

REAL (dp), INTENT(IN) :: vm(:), xopt(:), x(:)
INTEGER, INTENT(IN)   :: n

WRITE(*,'(/, " intermediate results after step length adjustment", /)')
CALL prtvec(vm, n, 'NEW STEP LENGTH (VM)')
CALL prtvec(xopt, n, 'CURRENT OPTIMAL X')
CALL prtvec(x, n, 'CURRENT X')
WRITE(*,'(" ")')

RETURN
END SUBROUTINE prt8


SUBROUTINE prt9(max, n, t, xopt, vm, fopt, nup, ndown, nrej, lnobds, nnew)

REAL (dp), INTENT(IN) :: xopt(:), vm(:), t, fopt
INTEGER, INTENT(IN)   :: n, nup, ndown, nrej, lnobds, nnew
LOGICAL, INTENT(IN)   :: max

! Local variable
INTEGER :: totmov

totmov = nup + ndown + nrej

WRITE(*,'(/," intermediate results before next temperature reduction",/)')
WRITE(*,'("  CURRENT TEMPERATURE:            ",G12.5)') t
IF (max) THEN
  WRITE(*, '("  MAX FUNCTION VALUE SO FAR:  ",G25.18)') fopt
  WRITE(*, '("  TOTAL MOVES:                ",I8)') totmov
  WRITE(*, '("     UPHILL:                  ",I8)') nup
  WRITE(*, '("     ACCEPTED DOWNHILL:       ",I8)') ndown
  WRITE(*, '("     REJECTED DOWNHILL:       ",I8)') nrej
  WRITE(*, '("  OUT OF BOUNDS TRIALS:       ",I8)') lnobds
  WRITE(*, '("  NEW MAXIMA THIS TEMPERATURE:",I8)') nnew
ELSE
  WRITE(*, '("  MIN FUNCTION VALUE SO FAR:  ",G25.18)') -fopt
  WRITE(*, '("  TOTAL MOVES:                ",I8)') totmov
  WRITE(*, '("     DOWNHILL:                ",I8)')  nup
  WRITE(*, '("     ACCEPTED UPHILL:         ",I8)')  ndown
  WRITE(*, '("     REJECTED UPHILL:         ",I8)')  nrej
  WRITE(*, '("  TRIALS OUT OF BOUNDS:       ",I8)')  lnobds
  WRITE(*, '("  NEW MINIMA THIS TEMPERATURE:",I8)')  nnew
END IF
CALL prtvec(xopt, n, 'CURRENT OPTIMAL X')
CALL prtvec(vm, n, 'STEP LENGTH (VM)')
WRITE(*, '(" ")')

RETURN
END SUBROUTINE prt9


SUBROUTINE prt10()

WRITE(*, '(/, "  SA ACHIEVED TERMINATION CRITERIA. IER = 0. ",/)')

RETURN
END SUBROUTINE prt10


SUBROUTINE prtvec(vector, ncols, name)
!  This subroutine prints the double precision vector named VECTOR.
!  Elements 1 thru NCOLS will be printed. NAME is a character variable
!  that describes VECTOR. Note that if NAME is given in the call to
!  PRTVEC, it must be enclosed in quotes. If there are more than 10
!  elements in VECTOR, 10 elements will be printed on each line.

INTEGER, INTENT(IN)           :: ncols
REAL (dp), INTENT(IN)         :: vector(ncols)
CHARACTER (LEN=*), INTENT(IN) :: name

INTEGER :: i, lines, ll

WRITE(*,1001) NAME

IF (ncols > 10) THEN
  lines = INT(ncols/10.)
  
  DO i = 1, lines
    ll = 10*(i - 1)
    WRITE(*,1000) vector(1+ll:10+ll)
  END DO
  
  WRITE(*,1000) vector(11+ll:ncols)
ELSE
  WRITE(*,1000) vector(1:ncols)
END IF

1000 FORMAT( 10(g12.5, ' '))
1001 FORMAT(/, 25(' '), a)

RETURN
END SUBROUTINE prtvec

END MODULE simulated_anneal



PROGRAM simann

USE simulated_anneal
IMPLICIT NONE
INTEGER, PARAMETER :: n = 2, neps = 4

REAL (dp)   :: lb(n), ub(n), x(n), xopt(n), c(n), vm(n), t, eps, rt, fopt

INTEGER     :: ns, nt, nfcnev, ier, iseed1, iseed2, i, maxevl, iprint,  &
               nacc, nobds

LOGICAL     :: max

!  Set underflows to zero on IBM mainframes.
!     CALL XUFLOW(0)

!  Set input parameters.
max = .false.
eps = 1.0D-6
rt = .5
iseed1 = 1
iseed2 = 2
ns = 20
nt = 5
maxevl = 100000
iprint = 1
DO i = 1, n
  lb(i) = -1.0D25
  ub(i) =  1.0D25
  c(i) = 2.0
END DO

!  Note start at local, but not global, optima of the Judge function.
x(1) =  2.354471
x(2) = -0.319186

!  Set input values of the input/output parameters.
t = 5.0
vm(1:n) = 1.0

WRITE(*,1000) n, max, t, rt, eps, ns, nt, neps, maxevl, iprint, iseed1, iseed2

CALL prtvec(x, n, 'STARTING VALUES')
CALL prtvec(vm, n, 'INITIAL STEP LENGTH')
CALL prtvec(lb, n, 'LOWER BOUND')
CALL prtvec(ub, n, 'UPPER BOUND')
CALL prtvec(c, n, 'C VECTOR')
WRITE(*, '(/, "  ****   END OF DRIVER ROUTINE OUTPUT   ****"/,  &
      &     "  ****   before CALL TO sa.             ****")')

CALL sa(n, x, max, rt, eps, ns, nt, neps, maxevl, lb, ub, c, iprint, iseed1,  &
        iseed2, t, vm, xopt, fopt, nacc, nfcnev, nobds, ier)

WRITE(*, '(/, "  ****   RESULTS AFTER SA   ****   ")')
CALL prtvec(xopt, n, 'SOLUTION')
CALL prtvec(vm, n, 'FINAL STEP LENGTH')
WRITE(*,1001) fopt, nfcnev, nacc, nobds, t, ier

1000 FORMAT(/,' SIMULATED ANNEALING EXAMPLE',/,/,  &
             ' NUMBER OF PARAMETERS: ',i3,'   MAXIMIZATION: ',l5, /, &
             ' INITIAL TEMP: ', g8.2, '   RT: ',g8.2, '   EPS: ',g8.2, /, &
             ' NS: ',i3, '   NT: ',i2, '   NEPS: ',i2, /,  &
             ' MAXEVL: ',i10, '   IPRINT: ',i1, '   ISEED1: ',i4,  &
             '   ISEED2: ',i4)
1001 FORMAT(/,' OPTIMAL FUNCTION VALUE: ',g20.13  &
            /,' NUMBER OF FUNCTION EVALUATIONS:     ',i10,  &
            /,' NUMBER OF ACCEPTED EVALUATIONS:     ',i10,  &
            /,' NUMBER OF OUT OF BOUND EVALUATIONS: ',i10,  &
            /,' FINAL TEMP: ', g20.13,'  IER: ', i3)

STOP
END PROGRAM simann


SUBROUTINE fcn(n, theta, h)
!  This subroutine is from the example in Judge et al., The Theory and
!  Practice of Econometrics, 2nd ed., pp. 956-7. There are two optima:
!  F(.864,1.23) = 16.0817 (the global minumum) and F(2.35,-.319) = 20.9805.

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: theta(:)
REAL (dp), INTENT(OUT) :: h

! Local variables
INTEGER   :: i
REAL (dp) :: y(20), x2(20), x3(20)

y(1) = 4.284
y(2) = 4.149
y(3) = 3.877
y(4) = 0.533
y(5) = 2.211
y(6) = 2.389
y(7) = 2.145
y(8) = 3.231
y(9) = 1.998
y(10) = 1.379
y(11) = 2.106
y(12) = 1.428
y(13) = 1.011
y(14) = 2.179
y(15) = 2.858
y(16) = 1.388
y(17) = 1.651
y(18) = 1.593
y(19) = 1.046
y(20) = 2.152

x2(1) =  .286
x2(2) =  .973
x2(3) =  .384
x2(4) =  .276
x2(5) =  .973
x2(6) =  .543
x2(7) =  .957
x2(8) =  .948
x2(9) =  .543
x2(10) =  .797
x2(11) =  .936
x2(12) =  .889
x2(13) =  .006
x2(14) =  .828
x2(15) =  .399
x2(16) =  .617
x2(17) =  .939
x2(18) =  .784
x2(19) =  .072
x2(20) =  .889

x3(1) = .645
x3(2) = .585
x3(3) = .310
x3(4) = .058
x3(5) = .455
x3(6) = .779
x3(7) = .259
x3(8) = .202
x3(9) = .028
x3(10) = .099
x3(11) = .142
x3(12) = .296
x3(13) = .175
x3(14) = .180
x3(15) = .842
x3(16) = .039
x3(17) = .103
x3(18) = .620
x3(19) = .158
x3(20) = .704

h = 0.0_dp
DO i = 1, 20
  h = (theta(1) + theta(n)*x2(i) + (theta(n)**2)*x3(i) - y(i))**2 + h
END DO

RETURN
END SUBROUTINE fcn
