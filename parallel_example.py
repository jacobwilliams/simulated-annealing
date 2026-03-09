#!/usr/bin/env python
"""
Example using the parallel interface with multiprocessing.
Demonstrates how to evaluate multiple function evaluations in parallel.

This example shows how to implement custom parallelization by:
1. Creating parallel callback functions
2. Setting up a worker pool (using multiprocessing here)
3. Passing the callbacks to sa_fortran.initialize()

Note: the shared library must be built first using:
`fpm install --prefix ./simulated_annealing_fortran/lib --profile release`
"""

import numpy as np
import multiprocessing as mp
from multiprocessing import Queue, Process
import time
import sys

# Import the class and callback types
from simulated_annealing_fortran import (sa_fortran,
                                         CALLBACK_N_INPUTS,
                                         CALLBACK_PARALLEL_INPUT,
                                         CALLBACK_PARALLEL_OUTPUT)

def rastrigin(x):
    """
    Rastrigin function - multimodal test function.
    Global minimum: f(0, 0, ..., 0) = 0
    """
    A = 10.0
    n = len(x)
    return A * n + sum(xi**2 - A * np.cos(2 * np.pi * xi) for xi in x)

################
# from scipy:
################
# import numpy as np
# from scipy.optimize import dual_annealing
# func = lambda x: np.sum(x*x - 10*np.cos(2*np.pi*x)) + 10*np.size(x)
# lw = [-5.12] * 10
# up = [5.12] * 10
# ret = dual_annealing(func, bounds=list(zip(lw, up)))
# ret.x
# array([-4.26437714e-09, -3.91699361e-09, -1.86149218e-09, -3.97165720e-09,
#        -6.29151648e-09, -6.53145322e-09, -3.93616815e-09, -6.55623025e-09,
#        -6.05775280e-09, -5.00668935e-09]) # random
# ret.fun
# 0.000000



def worker_process(input_queue, output_queue):
    """
    Worker process that evaluates functions in parallel.
    Reads (job_id, x) from input_queue, computes f, writes (job_id, x, f) to output_queue.
    """
    while True:
        item = input_queue.get()
        if item is None:  # Poison pill
            break

        job_id, x = item
        try:
            f = rastrigin(x)
            output_queue.put((job_id, x, f, 0))  # 0 = success
        except Exception as e:
            print(f"Worker error on job {job_id}: {e}")
            output_queue.put((job_id, x, 0.0, -1))  # -1 = error


class ParallelEvaluator:
    """Manages parallel function evaluation using multiprocessing."""

    def __init__(self, n_workers=4):
        self.n_workers = n_workers
        self.input_queue = Queue()  # No maxsize limit to avoid blocking in callbacks
        self.output_queue = Queue()
        self.workers = []
        self.job_counter = 0
        self.pending_jobs = {}

        # Start worker processes
        for i in range(n_workers):
            p = Process(target=worker_process, args=(self.input_queue, self.output_queue))
            p.start()
            self.workers.append(p)

    def submit_batch(self, x_batch):
        """
        Submit a batch of inputs for evaluation.
        x_batch: (n_inputs, n) array
        """
        n_inputs = x_batch.shape[0]
        job_ids = []

        for i in range(n_inputs):
            job_id = self.job_counter
            self.job_counter += 1
            x = x_batch[i, :].copy()
            self.input_queue.put((job_id, x))
            self.pending_jobs[job_id] = True
            job_ids.append(job_id)

        return job_ids

    def get_result(self, timeout=10.0):
        """Get one completed result (blocking)."""
        try:
            job_id, x, f, istat = self.output_queue.get(timeout=timeout)
            if job_id in self.pending_jobs:
                del self.pending_jobs[job_id]
            return x, f, istat
        except:
            return None, 0.0, -2  # timeout

    def shutdown(self):
        """Shutdown worker processes."""
        # Send poison pills
        for _ in range(self.n_workers):
            self.input_queue.put(None)

        # Wait for workers to finish (short timeout, then terminate)
        for p in self.workers:
            p.join(timeout=0.5)
            if p.is_alive():
                p.terminate()
                p.join()

        # Cancel queue threads to avoid blocking on cleanup
        self.input_queue.cancel_join_thread()
        self.output_queue.cancel_join_thread()

        # Close queues
        self.input_queue.close()
        self.output_queue.close()


def main():
    n = 1  # Problem dimension
    n_workers = 4  # Number of parallel workers

    print(f"Solving {n}-dimensional Rastrigin function with {n_workers} parallel workers")
    print(f"Global minimum at x = [0, 0, ..., 0] with f = 0")
    print()

    # Create the parallel evaluator (worker pool)
    evaluator = ParallelEvaluator(n_workers=n_workers)

    # Define the parallel callbacks
    @CALLBACK_N_INPUTS
    def n_inputs_callback(iproblem, n_inputs_ptr):
        """Tell Fortran how many inputs we can process at once."""
        n_inputs_ptr[0] = n_workers

    @CALLBACK_PARALLEL_INPUT
    def parallel_input_callback(iproblem, x_ptr, n_val, n_inputs_val):
        """Receive batch of inputs and submit to worker pool."""
        total_size = n_val * n_inputs_val
        x_flat = np.ctypeslib.as_array(x_ptr, shape=(total_size,))

        # Reshape: Fortran's x(n, n_inputs) in column-major order
        x_batch = np.zeros((n_inputs_val, n_val))
        for i in range(n_inputs_val):
            x_batch[i, :] = x_flat[i * n_val:(i + 1) * n_val]

        # Submit to worker pool
        evaluator.submit_batch(x_batch)

    @CALLBACK_PARALLEL_OUTPUT
    def parallel_output_callback(iproblem, x_ptr, n_val, f_ptr, istat_ptr):
        """Get one result from the completed queue."""
        x_result, f, istat = evaluator.get_result()

        if istat >= 0 and x_result is not None:
            x_array = np.ctypeslib.as_array(x_ptr, shape=(n_val,))
            x_array[:] = x_result
            f_ptr[0] = f
            istat_ptr[0] = istat
        else:
            f_ptr[0] = 0.0
            istat_ptr[0] = istat

    # Create optimizer and initialize in parallel mode
    optimizer = sa_fortran()
    optimizer.initialize(
        n=n,
        lb=[-5.12] * n,
        ub=[5.12] * n,
        fcn=None,  # Not used in parallel mode
        c=[2.0] * n,  # Step adjustment factor (not used in constant mode)
        maximize=False,
        eps=1e-2,  # Convergence tolerance
        ns=20,  # Cycles before step adjustment
        nt=100,  # Iterations before temperature reduction
        neps=4,  # Number of final function values for convergence check
        maxevl=10000,  # Function evaluations budget
        iprint=2,
        iseed1=1234,
        iseed2=5678,
        step_mode=1,  # 1=adaptive step size (default)
        cooling_schedule=2,  # 2=fast annealing (terminates quicker)
        n_resets=3,  # Restart search from best point found
        optimal_f_specified=False,  # Don't use optimal stopping criterion
        # Parallel mode callbacks
        n_inputs_to_send=n_inputs_callback,
        fcn_parallel_input=parallel_input_callback,
        fcn_parallel_output=parallel_output_callback,
    )

    print("Optimizer initialized in parallel mode.")
    print()

    # Solve
    start_time = time.time()
    sys.stdout.flush()  # Flush Python output before Fortran writes

    result = optimizer.solve(
        x0=[1.0] * n,
        rt=0.85,  # Standard cooling rate
        t0=1.0,  # Initial temperature
        vm=[2.0] * n,  # Fixed step size for each variable
    )

    sys.stdout.flush()  # Flush after Fortran is done
    elapsed = time.time() - start_time

    # Results
    print()
    print("=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    print(f"Exit code (ier): {result['ier']}")
    print(f"Function evaluations: {result['nfcnev']}")
    print(f"Accepted moves: {result['nacc']}")
    print(f"Final temperature: {result['t']:.6e}")
    print(f"Optimal function value: {result['f']:.10e}")
    print(f"Optimal x: {result['x']}")
    print(f"Wall time: {elapsed:.2f} seconds")
    print(f"Time per evaluation: {elapsed/result['nfcnev']*1000:.3f} ms")
    print()

    error = np.linalg.norm(result['x'])
    print(f"Distance from true optimum [0,0,...,0]: {error:.6e}")

    # Cleanup
    optimizer.destroy()
    evaluator.shutdown()
    print("\nOptimizer destroyed, workers shutdown.")


if __name__ == "__main__":
    # Required for multiprocessing on macOS
    mp.set_start_method('spawn', force=True)
    main()
