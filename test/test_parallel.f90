!********************************************************************************
!>
!  Test program for parallel simulated annealing.
!
!  This demonstrates how to use the parallel evaluation interface to
!  manage a pool of workers. Results are stored in a FIFO queue and
!  returned in the order they were completed.

program test_parallel

use iso_fortran_env, only: output_unit
use simulated_annealing_module, only: simulated_annealing_type, dp => simann_wp, print_vector

implicit none

integer, parameter :: n = 2  ! number of variables
integer, parameter :: max_pool_size = 10  ! maximum number of workers

! worker pool data structures:
type :: worker_type
    real(dp), dimension(n) :: x      ! input point
    real(dp) :: f                     ! function value
    logical :: busy                   ! is this worker busy?
end type worker_type

! result queue for completed evaluations:
type :: result_type
    real(dp), dimension(n) :: x      ! input point
    real(dp) :: f                     ! function value
end type result_type

type(worker_type), dimension(max_pool_size) :: workers
type(result_type), dimension(max_pool_size*10) :: result_queue  ! queue for completed results
integer :: n_busy_workers
integer :: queue_head, queue_tail, queue_size  ! queue management

! SA parameters:
real(dp)   :: lb(n), ub(n), x(n), xopt(n), c(n), vm(n), t, eps, rt, fopt
integer    :: ns, nt, nfcnev, ier, iseed1, iseed2, i, maxevl, iprint, nacc, n_resets

type(simulated_annealing_type) :: sa

! initialize worker pool and result queue:
workers%busy = .false.
n_busy_workers = 0
queue_head = 1
queue_tail = 1
queue_size = 0

! set SA parameters for Rosenbrock function:
! minimum at (1,1) with f=0
eps = 1.0e-6_dp
rt = 0.85_dp
iseed1 = 1
iseed2 = 2
ns = 20
nt = 5
maxevl = 10000
iprint = 1
n_resets = 1

! bounds:
lb(1:n) = -5.0_dp
ub(1:n) = 5.0_dp
c(1:n) = 2.0_dp

! starting point:
x(1) = 0.0_dp
x(2) = 0.0_dp

! initial temperature and step length:
t = 5.0_dp
vm(1:n) = 1.0_dp

write(output_unit, '(/A/)') ' =========================================='
write(output_unit, '(A)')   ' Parallel Simulated Annealing Example'
write(output_unit, '(A)')   ' Test function: Rosenbrock'
write(output_unit, '(A)')   ' Global minimum: f(1,1) = 0'
write(output_unit, '(A/)')  ' =========================================='
write(output_unit, '(A,I3)') '   Number of parameters:  ', n
write(output_unit, '(A,I3)') '   Max pool size:         ', max_pool_size
write(output_unit, '(A,G9.2)') '   Initial temperature:   ', t
write(output_unit, '(A,G9.2)') '   Temperature reduction: ', rt
write(output_unit, '(A,G9.2)') '   Convergence tolerance: ', eps
write(output_unit, '(A,I4)')   '   NS (cycle length):     ', ns
write(output_unit, '(A,I4)')   '   NT (temp reduction):   ', nt
write(output_unit, '(A,I10)')  '   Max evaluations:       ', maxevl
write(output_unit, '(A/)')     ' =========================================='

call print_vector(output_unit, x, n, 'Starting point')
call print_vector(output_unit, lb, n, 'Lower bounds')
call print_vector(output_unit, ub, n, 'Upper bounds')

! initialize SA with parallel function evaluation:
call sa%initialize(n=n, &
                   lb=lb, &
                   ub=ub, &
                   c=c, &
                   maximize=.false., &
                   eps=eps, &
                   ns=ns, &
                   nt=nt, &
                   maxevl=maxevl, &
                   iprint=iprint, &
                   iseed1=iseed1, &
                   iseed2=iseed2, &
                   n_resets=n_resets,&
                   use_initial_guess=.true., &
                   n_inputs_to_send=get_n_available_workers, &
                   fcn_parallel_input=send_work_to_pool, &
                   fcn_parallel_output=get_result_from_pool)

! run optimization:
call sa%optimize(x, rt, t, vm, xopt, fopt, nacc, nfcnev, ier)

! print results:
write(output_unit, '(/A/)')  ' =========================================='
write(output_unit, '(A)')    '   Final Results'
write(output_unit, '(A/)')   ' =========================================='
write(output_unit, '(A,I4)') '   Exit code:             ', ier
write(output_unit, '(A,I10)') '   Function evaluations:  ', nfcnev
write(output_unit, '(A,I10)') '   Accepted moves:        ', nacc
write(output_unit, '(A,G25.18)') '   Final f value:         ', fopt
call print_vector(output_unit, xopt, n, 'Optimal x')
call print_vector(output_unit, vm, n, 'Final step length')

if (ier == 0) then
    write(output_unit, '(/A)') ' Optimization successful!'
    write(output_unit, '(A,G12.5)') ' Error from global minimum: ', &
        sqrt((xopt(1)-1.0_dp)**2 + (xopt(2)-1.0_dp)**2)
else
    write(output_unit, '(/A,I4)') ' Optimization failed with error code: ', ier
end if

contains

    !******************************************************************************
    subroutine get_n_available_workers(me, n_inputs)
        !! Return the number of available workers in the pool
        !! (simulates checking how many workers are idle)
        class(simulated_annealing_type), intent(inout) :: me
        integer, intent(out) :: n_inputs

        integer :: i, n_available

        ! simulate some workers finishing their work:
        ! in a real implementation, this would check actual worker status
        call simulate_worker_completion()

        ! count available workers:
        n_available = 0
        do i = 1, max_pool_size
            if (.not. workers(i)%busy) then
                n_available = n_available + 1
            end if
        end do

        ! return at least 1, but typically we want to keep the pool busy
        n_inputs = max(1, min(3, n_available))  ! send up to 3 at a time

    end subroutine get_n_available_workers
    !******************************************************************************

    !******************************************************************************
    subroutine send_work_to_pool(me, x_array)
        !! Send work to the worker pool
        !! In a real implementation, this would dispatch to actual workers
        class(simulated_annealing_type), intent(inout) :: me
        real(dp), dimension(:,:), intent(in) :: x_array  ! (n, n_inputs)

        integer :: i, j

        ! find available workers and assign them work:
        j = 0
        do i = 1, max_pool_size
            if (.not. workers(i)%busy) then
                j = j + 1
                if (j > size(x_array, 2)) exit

                ! assign work to this worker:
                workers(i)%x(:) = x_array(:, j)
                workers(i)%busy = .true.
                n_busy_workers = n_busy_workers + 1

                ! simulate starting the evaluation:
                ! in a real implementation, this would send the work to a GPU,
                ! another process, or thread pool
            end if
        end do

    end subroutine send_work_to_pool
    !******************************************************************************

    !******************************************************************************
    subroutine get_result_from_pool(me, x, f, istat)
        !! Get a completed result from the result queue
        !! Returns results in FIFO order (first completed, first returned)
        class(simulated_annealing_type), intent(inout) :: me
        real(dp), dimension(:), intent(out) :: x  ! returns the x that was evaluated
        real(dp), intent(out) :: f
        integer, intent(out) :: istat

        ! wait for at least one result to be in the queue:
        do while (queue_size == 0)
            call simulate_worker_completion()
        end do

        ! dequeue the oldest result (FIFO):
        x = result_queue(queue_head)%x
        f = result_queue(queue_head)%f
        istat = 0  ! success

        ! advance queue head (with wraparound):
        queue_head = queue_head + 1
        if (queue_head > size(result_queue)) queue_head = 1
        queue_size = queue_size - 1

    end subroutine get_result_from_pool
    !******************************************************************************

    !******************************************************************************
    subroutine simulate_worker_completion()
        !! Simulate workers completing their evaluations and add results to queue
        !! In a real implementation, this would check actual worker status
        integer :: i
        real(dp) :: rand

        ! randomly complete some busy workers:
        do i = 1, max_pool_size
            if (workers(i)%busy) then
                call random_number(rand)
                if (rand > 0.5_dp) then  ! 50% chance of completion
                    ! evaluate the function:
                    workers(i)%f = rosenbrock(workers(i)%x)

                    ! add result to queue:
                    result_queue(queue_tail)%x = workers(i)%x
                    result_queue(queue_tail)%f = workers(i)%f
                    queue_tail = queue_tail + 1
                    if (queue_tail > size(result_queue)) queue_tail = 1
                    queue_size = queue_size + 1

                    ! mark worker as available:
                    workers(i)%busy = .false.
                    n_busy_workers = n_busy_workers - 1
                end if
            end if
        end do

    end subroutine simulate_worker_completion
    !******************************************************************************

    !******************************************************************************
    function rosenbrock(x) result(f)
        !! Rosenbrock function: f(x) = (1-x1)^2 + 100*(x2-x1^2)^2
        !! Global minimum at (1,1) with f=0
        real(dp), dimension(:), intent(in) :: x
        real(dp) :: f

        f = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2

    end function rosenbrock
    !******************************************************************************

end program test_parallel
