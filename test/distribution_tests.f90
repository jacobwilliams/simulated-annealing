program distribution_tests

!! generate a data file with the results from perturbing variables
!! using all the different distributions. This is used to test the distribution
!! functions in the simulated annealing module.

use iso_fortran_env, only: output_unit
use simulated_annealing_module, only: simulated_annealing_type, dp => simann_wp

implicit none

type(simulated_annealing_type) :: sa
integer, parameter :: n_samples = 10000  !! number of samples to generate per distribution
integer, parameter :: n_dists = 10      !! number of distribution types
real(dp), dimension(n_samples, n_dists) :: samples
real(dp) :: lower, upper
integer :: i, j, k, iunit
character(len=100), dimension(n_dists) :: dist_names
character(len=*), parameter :: csv_file = 'distribution_samples.csv'

! Distribution names for header
dist_names(1)  = 'uniform'
dist_names(2)  = 'normal'
dist_names(3)  = 'cauchy'
dist_names(4)  = 'exponential'
dist_names(5)  = 'pareto'
dist_names(6)  = 'truncated_normal'
dist_names(7)  = 'beta'
dist_names(8)  = 'triangular'
dist_names(9)  = 'kumaraswamy'
dist_names(10) = 'bipareto'

! Set bounds for perturbation
lower = -5.0_dp
upper = 5.0_dp

write(output_unit, '(A)') 'Generating distribution samples...'
write(output_unit, '(A,I0,A)') 'Number of samples per distribution: ', n_samples
write(output_unit, '(A,F0.2,A,F0.2,A)') 'Bounds: [', lower, ', ', upper, ']'
write(output_unit, '(A)') ''

! Generate samples for each distribution
do j = 1, n_dists
    write(output_unit, '(A,I0,A,A)') 'Generating samples for distribution ', j, ': ', trim(dist_names(j))

    ! Initialize with appropriate distribution and parameters for each type
    select case(j)
    case(1)  ! uniform
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[1])
    case(2)  ! normal
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[2], dist_mean=[0.0_dp], dist_std_dev=[2.0_dp])
    case(3)  ! cauchy
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[3], dist_location=[0.0_dp], dist_scale=[1.0_dp])
    case(4)  ! exponential
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[4], dist_rate=[1.0_dp])
    case(5)  ! pareto
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[5], dist_scale=[1.0_dp], dist_shape=[2.0_dp])
    case(6)  ! truncated_normal
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[6], dist_mean=[0.0_dp], dist_std_dev=[2.0_dp])
    case(7)  ! beta
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[7], dist_alpha=[2.0_dp], dist_beta=[5.0_dp])
    case(8)  ! triangular
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[8], dist_mode=[0.3_dp])
    case(9)  ! kumaraswamy
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[9], dist_a=[2.0_dp], dist_b=[5.0_dp])
    case(10)  ! bipareto
        call sa%initialize(fcn=dummy_fcn, n=1, lb=[lower], ub=[upper], &
                          distribution_mode=[10], dist_scale=[1.0_dp], dist_shape=[2.0_dp])
    end select

    ! Generate samples
    do i = 1, n_samples
        samples(i, j) = sa%perturb_variable(1, lower, upper)
    end do

    ! Clean up for next iteration
    call sa%destroy()
end do

! Write to CSV file
write(output_unit, '(A)') ''
write(output_unit, '(A)') 'Writing samples to ' // csv_file

open(newunit=iunit, file=csv_file, status='replace', action='write')

! Write header
write(iunit, '(A)', advance='no') 'sample'
do k = 1, n_dists
    write(iunit, '(A,A)', advance='no') ',', trim(dist_names(k))
end do
write(iunit, '(A)') ''

! Write data
do i = 1, n_samples
    write(iunit, '(I0)', advance='no') i
    do k = 1, n_dists
        write(iunit, '(A,ES15.8)', advance='no') ',', samples(i, k)
    end do
    write(iunit, '(A)') ''
end do

close(iunit)

write(output_unit, '(A)') 'Done!'
write(output_unit, '(A)') ''
write(output_unit, '(A)') 'Statistics summary:'
write(output_unit, '(A)') '-------------------'
do k = 1, n_dists
        write(output_unit, '(A,A)') trim(dist_names(k)), ':'
        write(output_unit, '(A,F10.4)') '  Mean:   ', sum(samples(:, k)) / n_samples
        write(output_unit, '(A,F10.4)') '  Min:    ', minval(samples(:, k))
        write(output_unit, '(A,F10.4)') '  Max:    ', maxval(samples(:, k))
        write(output_unit, '(A,G10.4)') '  StdDev: ', sqrt(sum((samples(:, k) - sum(samples(:, k))/n_samples)**2) / n_samples)
    end do

contains

    subroutine dummy_fcn(me, x, f, istat)
        !! Dummy function (not actually used for this test)
        use simulated_annealing_module, only: simulated_annealing_type, dp => simann_wp
        implicit none
        class(simulated_annealing_type), intent(inout) :: me
        real(dp), dimension(:), intent(in) :: x
        real(dp), intent(out) :: f
        integer, intent(out) :: istat

        f = sum(x**2)  ! simple quadratic
        istat = 0
    end subroutine dummy_fcn

end program distribution_tests