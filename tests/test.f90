program simann

use iso_fortran_env, only: output_unit
use simulated_annealing_module, only: simulated_annealing_type, dp, print_vector

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

call print_vector(output_unit,x, n, 'starting values')
call print_vector(output_unit,vm, n, 'initial step length')
call print_vector(output_unit,lb, n, 'lower bound')
call print_vector(output_unit,ub, n, 'upper bound')
call print_vector(output_unit,c, n, 'c vector')
write(output_unit, '(A)') '  ****   end of driver routine output   ****'
write(output_unit, '(A)') '  ****   before call to sa.             ****'

call sa%initialize(fcn,n,lb,ub,c,&
                   max,eps,ns,nt,neps,maxevl,&
                   iprint,iseed1,iseed2,step_mode,vms,iunit)

call sa%optimize(x, rt, t, vm, xopt, fopt, nacc, nfcnev, ier)

write(output_unit, '(A)') '  ****   results after sa   ****   '
call print_vector(output_unit,xopt, n, 'solution')
call print_vector(output_unit,vm, n, 'final step length')
write(output_unit,'(/A,G20.13,/A,I10,/A,I10,/A,I10,/A,G20.13,/A,I3/)') &
      ' optimal function value: ', fopt, &
      ' number of function evaluations:     ', nfcnev,&
      ' number of accepted evaluations:     ', nacc, &
      ' number of out of bound evaluations: ', nobds, &
      ' final temp: ', t, &
      '  ier: ', ier

contains

    subroutine fcn(me, theta, h, istat)
    !!  this subroutine is from the example in judge et al., the theory and
    !!  practice of econometrics, 2nd ed., pp. 956-7. there are two optima:
    !!  f(.864,1.23) = 16.0817 (the global minumum) and f(2.35,-.319) = 20.9805.

    implicit none

    class(simulated_annealing_type),intent(inout) :: me
    real (dp), dimension(:), intent(in)  :: theta
    real (dp), intent(out) :: h
    integer,intent(out) :: istat

    integer,parameter :: n = 20 !! size of arrays

    integer :: i !! counter
    real(dp),dimension(n),parameter :: y = [4.284_dp, &
                                            4.149_dp, &
                                            3.877_dp, &
                                            0.533_dp, &
                                            2.211_dp, &
                                            2.389_dp, &
                                            2.145_dp, &
                                            3.231_dp, &
                                            1.998_dp, &
                                            1.379_dp, &
                                            2.106_dp, &
                                            1.428_dp, &
                                            1.011_dp, &
                                            2.179_dp, &
                                            2.858_dp, &
                                            1.388_dp, &
                                            1.651_dp, &
                                            1.593_dp, &
                                            1.046_dp, &
                                            2.152_dp ]
    real(dp),dimension(n),parameter :: x2 = [.286_dp, &
                                             .973_dp, &
                                             .384_dp, &
                                             .276_dp, &
                                             .973_dp, &
                                             .543_dp, &
                                             .957_dp, &
                                             .948_dp, &
                                             .543_dp, &
                                             .797_dp, &
                                             .936_dp, &
                                             .889_dp, &
                                             .006_dp, &
                                             .828_dp, &
                                             .399_dp, &
                                             .617_dp, &
                                             .939_dp, &
                                             .784_dp, &
                                             .072_dp, &
                                             .889_dp ]

    real(dp),dimension(n),parameter :: x3 = [.645_dp, &
                                             .585_dp, &
                                             .310_dp, &
                                             .058_dp, &
                                             .455_dp, &
                                             .779_dp, &
                                             .259_dp, &
                                             .202_dp, &
                                             .028_dp, &
                                             .099_dp, &
                                             .142_dp, &
                                             .296_dp, &
                                             .175_dp, &
                                             .180_dp, &
                                             .842_dp, &
                                             .039_dp, &
                                             .103_dp, &
                                             .620_dp, &
                                             .158_dp, &
                                             .704_dp ]

    h = 0.0_dp
    do i = 1, n
        h = (theta(1) + theta(2)*x2(i) + (theta(2)**2)*x3(i) - y(i))**2 + h
    end do

    istat = 0

    end subroutine fcn

end program simann

