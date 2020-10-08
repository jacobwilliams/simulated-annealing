program simann

use iso_fortran_env, only: output_unit
use simulated_annealing_module, only: simulated_annealing_type, dp, prtvec

implicit none

integer, parameter :: n = 2, neps = 4

real (dp)   :: lb(n), ub(n), x(n), xopt(n), c(n), vm(n), t, eps, rt, fopt, vms
integer     :: ns, nt, nfcnev, ier, iseed1, iseed2, i, maxevl, iprint,  &
               nacc, nobds, step_mode, iunit
logical     :: max

type(simulated_annealing_type) :: sa

!  set input parameters.
max = .false.
eps = 1.0e-6_dp
rt = .5_dp
iseed1 = 1
iseed2 = 2
ns = 20
nt = 5
maxevl = 100000
iprint = 1
do i = 1, n
  lb(i) = -1.0e25_dp
  ub(i) =  1.0e25_dp
  c(i) = 2.0_dp
end do

nobds = 0 ! JW : no longer used

!  note start at local, but not global, optima of the judge function.
x(1) =  2.354471_dp
x(2) = -0.319186_dp

!  set input values of the input/output parameters.
t = 5.0_dp
vm(1:n) = 1.0_dp

! others:
step_mode = 1
vms = 0.1_dp
iunit = output_unit

write(output_unit,&
        '(A,//A,I3,/,A,L5,/A,G9.2,/A,G9.2,/A,G9.2,/A,I4,/A,I4,/A,I4,/A,I10,/A,I4,/A,I4,/A,I4)') &
        ' simulated annealing example',&
        '   number of parameters: ', n, &
        '   maximization:         ', max, &
        '   initial temp:         ', t, &
        '   rt:                   ', rt, &
        '   eps:                  ', eps, &
        '   ns:                   ', ns, &
        '   nt:                   ', nt, &
        '   neps:                 ', neps, &
        '   maxevl:               ', maxevl, &
        '   iprint:               ', iprint, &
        '   iseed1:               ', iseed1, &
        '   iseed2:               ', iseed2

call prtvec(output_unit,x, n, 'starting values')
call prtvec(output_unit,vm, n, 'initial step length')
call prtvec(output_unit,lb, n, 'lower bound')
call prtvec(output_unit,ub, n, 'upper bound')
call prtvec(output_unit,c, n, 'c vector')
write(output_unit, '(A)') '  ****   end of driver routine output   ****'
write(output_unit, '(A)') '  ****   before call to sa.             ****'

call sa%initialize(fcn,n,lb,ub,c,&
                   max,eps,ns,nt,neps,maxevl,&
                   iprint,iseed1,iseed2,step_mode,vms,iunit)

call sa%optimize(x, rt, t, vm, xopt, fopt, nacc, nfcnev, ier)

write(output_unit, '(A)') '  ****   results after sa   ****   '
call prtvec(output_unit,xopt, n, 'solution')
call prtvec(output_unit,vm, n, 'final step length')
write(output_unit,'(/A,G20.13,/A,I10,/A,I10,/A,I10,/A,G20.13,/A,I3/)') &
      ' optimal function value: ', fopt, &
      ' number of function evaluations:     ', nfcnev,&
      ' number of accepted evaluations:     ', nacc, &
      ' number of out of bound evaluations: ', nobds, &
      ' final temp: ', t, &
      '  ier: ', ier

contains

    subroutine fcn(me, theta, h, istat)
    !  this subroutine is from the example in judge et al., the theory and
    !  practice of econometrics, 2nd ed., pp. 956-7. there are two optima:
    !  f(.864,1.23) = 16.0817 (the global minumum) and f(2.35,-.319) = 20.9805.

    implicit none

    class(simulated_annealing_type),intent(inout) :: me
    real (dp), intent(in)  :: theta(:)
    real (dp), intent(out) :: h
    integer,intent(out) :: istat

    integer   :: i
    real (dp) :: y(20), x2(20), x3(20)

    y(1) = 4.284_dp
    y(2) = 4.149_dp
    y(3) = 3.877_dp
    y(4) = 0.533_dp
    y(5) = 2.211_dp
    y(6) = 2.389_dp
    y(7) = 2.145_dp
    y(8) = 3.231_dp
    y(9) = 1.998_dp
    y(10) = 1.379_dp
    y(11) = 2.106_dp
    y(12) = 1.428_dp
    y(13) = 1.011_dp
    y(14) = 2.179_dp
    y(15) = 2.858_dp
    y(16) = 1.388_dp
    y(17) = 1.651_dp
    y(18) = 1.593_dp
    y(19) = 1.046_dp
    y(20) = 2.152_dp

    x2(1) =  .286_dp
    x2(2) =  .973_dp
    x2(3) =  .384_dp
    x2(4) =  .276_dp
    x2(5) =  .973_dp
    x2(6) =  .543_dp
    x2(7) =  .957_dp
    x2(8) =  .948_dp
    x2(9) =  .543_dp
    x2(10) =  .797_dp
    x2(11) =  .936_dp
    x2(12) =  .889_dp
    x2(13) =  .006_dp
    x2(14) =  .828_dp
    x2(15) =  .399_dp
    x2(16) =  .617_dp
    x2(17) =  .939_dp
    x2(18) =  .784_dp
    x2(19) =  .072_dp
    x2(20) =  .889_dp

    x3(1) = .645_dp
    x3(2) = .585_dp
    x3(3) = .310_dp
    x3(4) = .058_dp
    x3(5) = .455_dp
    x3(6) = .779_dp
    x3(7) = .259_dp
    x3(8) = .202_dp
    x3(9) = .028_dp
    x3(10) = .099_dp
    x3(11) = .142_dp
    x3(12) = .296_dp
    x3(13) = .175_dp
    x3(14) = .180_dp
    x3(15) = .842_dp
    x3(16) = .039_dp
    x3(17) = .103_dp
    x3(18) = .620_dp
    x3(19) = .158_dp
    x3(20) = .704_dp

    h = 0.0_dp
    do i = 1, 20
        h = (theta(1) + theta(n)*x2(i) + (theta(n)**2)*x3(i) - y(i))**2 + h
    end do

    istat = 0

    end subroutine fcn

end program simann

