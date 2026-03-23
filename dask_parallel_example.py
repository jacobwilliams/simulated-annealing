#!/usr/bin/env python3
"""
Dask-based parallel example for sa_fortran library.

Demonstrates how to use Dask distributed computing with the parallel
interface of the sa_fortran library. This example uses a work queue
pattern where workers on Dask cluster nodes evaluate objective functions.

Features:
---------
- Works with both local clusters (development) and SLURM clusters (HPC)
- Uses distributed work and result queues for bidirectional communication
- Persistent workers that continuously process function evaluations
- Parallel batch evaluation of objective functions
- Real-time monitoring of function evaluations across cluster nodes

Architecture:
-------------
Main Process:
    - Initializes sa_fortran optimizer with parallel callbacks
    - Manages work_queue (sends x arrays to evaluate)
    - Manages result_queue (receives f values)
    - Callbacks bridge sa_fortran ↔ Dask queues

Worker Processes (on cluster nodes):
    - Pull x arrays from work_queue
    - Evaluate objective function f(x)
    - Send results back via result_queue
    - Continue until shutdown signal

Note: The shared library must be built first using:
`fpm install --prefix ./sa_fortran/lib --profile release`

Usage:
------
    # Local mode (default)
    python dask_parallel_example.py

    # Local mode with more workers
    python dask_parallel_example.py --mode local --workers 8

    # SLURM mode (requires dask-jobqueue and SLURM access)
    python dask_parallel_example.py --mode slurm --workers 12
"""

import argparse
import time
import numpy as np
import csv
import sys
import logging
from dask.distributed import Client, Queue

# Suppress expected Dask timeout warnings
logging.getLogger('distributed.core').setLevel(logging.CRITICAL)
logging.getLogger('distributed.queues').setLevel(logging.CRITICAL)

# Import the sa_fortran class and callback types
from sa_fortran import (sa_fortran,
                        CALLBACK_N_INPUTS,
                        CALLBACK_PARALLEL_INPUT,
                        CALLBACK_PARALLEL_OUTPUT)

# Configuration
PLOT_ITERATIONS = True  # Log function evaluations to CSV for later analysis

# Lists to store intermediate results
all_evaluations = []
best_values = []

def report_callback(x, f, istat):
    """
    Callback function to capture intermediate results.

    Parameters:
    -----------
    x : array
        Current point being evaluated
    f : float
        Function value at x
    istat : int
        Status: 1=function evaluation, 2=new optimal found
    """
    if istat == 1:
        # Regular function evaluation
        all_evaluations.append({'x': x.copy(), 'f': f})
    elif istat == 2:
        # New optimal value found
        best_values.append({'x': x.copy(), 'f': f})
        print(f"New best: f = {f:.6f}, x = {x}")


def rastrigin(x):
    """
    Rastrigin function - multimodal test function.
    Global minimum: f(0, 0, ..., 0) = 0

    This is a challenging test function with many local minima.

    Standard form: f(x) = A*n + Σ(xi² - A*cos(2π*xi))
    where A = 10.0
    """
    A = 10.0
    n = len(x)

    res = A * n + sum(xi**2 - A * np.cos(2 * np.pi * xi) for xi in x)

    if PLOT_ITERATIONS:
        # Log all function evaluations to CSV for later analysis
        with open('function_evaluations_dask.csv', 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(list(x) + [res])

    return res


def create_local_cluster(n_workers: int = 4):
    """
    Create a local Dask cluster for development and testing.

    Args:
        n_workers: Number of worker processes to start

    Returns:
        Client: Connected Dask client
    """
    print(f"Creating local cluster with {n_workers} workers...")
    client = Client(
        n_workers=n_workers,
        threads_per_worker=1,
        processes=True,
        dashboard_address=':8787'
    )
    print(f"Local cluster ready with {len(client.scheduler_info()['workers'])} workers")
    print(f"Dashboard: {client.dashboard_link}")

    ###################
    # leave this here for testing:  [open a dask dashboard in the browser]
    # import webbrowser
    # webbrowser.open(client.dashboard_link, new=0, autoraise=True)
    ###################

    return client


def create_slurm_cluster(n_workers: int = 12):
    """
    Create a SLURM cluster for HPC computing.

    Args:
        n_workers: Total number of Dask workers to request

    Returns:
        tuple: (cluster, client) - SLURMCluster and connected Client
    """
    from dask_jobqueue import SLURMCluster

    print(f"Configuring SLURM cluster for {n_workers} workers...")

    # Adjust these parameters for YOUR cluster!
    cluster = SLURMCluster(
        cores=4,
        memory="8GB",
        processes=4,
        walltime="00:30:00",
        queue="regular",  # CHANGE THIS to your partition name!
        name="dask-sa",
        job_extra_directives=[
            "--ntasks=1",
            "--cpus-per-task=4",
        ],
        scheduler_options={
            "dashboard_address": ":8787",
        },
    )

    # Calculate number of SLURM jobs needed
    workers_per_job = cluster.worker_spec['options']['processes']
    n_jobs = max(1, (n_workers + workers_per_job - 1) // workers_per_job)

    print(f"Submitting {n_jobs} SLURM jobs ({workers_per_job} workers/job)...")
    cluster.scale(jobs=n_jobs)

    print(f"\nConnecting to cluster...")
    client = Client(cluster)

    print(f"Waiting for SLURM jobs to start...")
    min_workers = min(4, n_workers)
    client.wait_for_workers(n_workers=min_workers, timeout=120)

    actual_workers = len(client.scheduler_info()['workers'])
    print(f"SLURM cluster ready with {actual_workers} workers")
    print(f"Dashboard: {client.dashboard_link}")

    return cluster, client


def worker_function_evaluator(worker_id: int, work_queue_name: str,
                              result_queue_name: str, sentinel: str = "STOP"):
    """
    Long-running worker that evaluates objective functions from work queue.

    This worker runs on a Dask cluster node and continuously pulls
    function evaluation requests from the distributed work queue.

    Work Item Format:
        {
            'job_id': int,
            'x': numpy array
        }

    Result Format:
        {
            'job_id': int,
            'x': numpy array,
            'f': float,
            'istat': int  (0 = success)
        }

    Args:
        worker_id: Unique identifier for this worker
        work_queue_name: Name of distributed Queue to pull work from
        result_queue_name: Name of distributed Queue to send results to
        sentinel: Value that signals worker to stop

    Returns:
        str: Summary message with tasks processed
    """
    from dask.distributed import get_client
    import socket

    # Connect to the Dask scheduler
    client = get_client()

    # Access distributed queues
    work_queue = Queue(work_queue_name, client=client)
    result_queue = Queue(result_queue_name, client=client)

    hostname = socket.gethostname()
    processed = 0

    # Signal ready
    result_queue.put({
        'type': 'ready',
        'worker_id': worker_id,
        'hostname': hostname
    })

    # Main processing loop
    while True:
        try:
            # Pull work from queue
            work_item = work_queue.get(timeout=2.0)

            # Check for shutdown signal
            if work_item == sentinel:
                result_queue.put({
                    'type': 'shutdown',
                    'worker_id': worker_id,
                    'hostname': hostname,
                    'processed': processed
                })
                break

            # Extract job information
            job_id = work_item['job_id']
            x = work_item['x']

            # Evaluate objective function
            try:
                f = rastrigin(x)
                istat = 0  # Success
            except Exception as e:
                print(f"Worker {worker_id} error evaluating job {job_id}: {e}")
                f = 0.0
                istat = -1  # Error

            # Send result back
            result_queue.put({
                'type': 'result',
                'job_id': job_id,
                'x': x,
                'f': f,
                'istat': istat,
                'worker_id': worker_id,
                'hostname': hostname
            })

            processed += 1

        except TimeoutError:
            # No work available, loop back
            continue
        except Exception as e:
            print(f"Worker {worker_id} unexpected error: {e}")
            continue

    return f"Worker {worker_id} on {hostname} processed {processed} evaluations"


class DaskParallelEvaluator:
    """
    Manages parallel function evaluation using Dask distributed queues.

    This class bridges the sa_fortran parallel callbacks with Dask's
    distributed computing infrastructure.
    """

    def __init__(self, client, n_workers):
        """
        Initialize the evaluator.

        Args:
            client: Dask distributed Client
            n_workers: Number of workers to use
        """
        self.client = client
        self.n_workers = n_workers
        self.job_counter = 0
        self.pending_jobs = {}
        self.pending_count = 0  # Track how many jobs are currently being processed

        # Create distributed queues
        self.work_queue = Queue('sa_work', client=self.client)
        self.result_queue = Queue('sa_results', client=self.client)

        # Submit long-running workers
        print(f"\nStarting {n_workers} persistent workers on cluster...")
        self.worker_futures = []
        for i in range(n_workers):
            future = self.client.submit(
                worker_function_evaluator,
                i,
                'sa_work',
                'sa_results',
                'STOP'
            )
            self.worker_futures.append(future)

        # Wait for workers to be ready
        print("Waiting for workers to initialize...")
        time.sleep(2)
        ready_count = 0
        for _ in range(n_workers * 2):  # Check multiple times
            try:
                msg = self.result_queue.get(timeout=1.0)
                if msg['type'] == 'ready':
                    ready_count += 1
                    print(f"  Worker {msg['worker_id']} ready on {msg['hostname']}")
            except:
                break

        print(f"✓ {ready_count} workers ready\n")

    def get_num_idle_workers(self):
        """
        Get the number of workers currently waiting for work.

        Returns:
            int: Number of idle workers (total workers - busy workers)
        """
        busy_workers = min(self.pending_count, self.n_workers)
        idle_workers = max(0, self.n_workers - busy_workers)
        return idle_workers

    def submit_batch(self, x_batch):
        """
        Submit a batch of inputs for evaluation.

        Args:
            x_batch: (n_inputs, n) array of input vectors

        Returns:
            list: Job IDs for submitted work
        """
        n_inputs = x_batch.shape[0]
        job_ids = []

        for i in range(n_inputs):
            job_id = self.job_counter
            self.job_counter += 1
            x = x_batch[i, :].copy()

            # Submit to work queue
            self.work_queue.put({
                'job_id': job_id,
                'x': x
            })

            self.pending_jobs[job_id] = True
            self.pending_count += 1
            job_ids.append(job_id)

        return job_ids

    def get_result(self, timeout=10.0):
        """
        Get one completed result (blocking).

        Args:
            timeout: Maximum time to wait in seconds

        Returns:
            tuple: (x, f, istat) or (None, 0.0, -2) on timeout
        """
        # Keep trying to get a result message, skipping non-result messages
        start_time = time.time()
        while True:
            remaining_time = timeout - (time.time() - start_time)
            if remaining_time <= 0:
                print("WARNING: get_result timeout!")
                return None, 0.0, -2

            try:
                msg = self.result_queue.get(timeout=min(remaining_time, 1.0))

                if msg['type'] == 'result':
                    job_id = msg['job_id']
                    if job_id in self.pending_jobs:
                        del self.pending_jobs[job_id]
                        self.pending_count -= 1  # Worker is now idle

                    return msg['x'], msg['f'], msg['istat']
                else:
                    # Skip non-result messages (ready, shutdown, etc.)
                    # Keep looping to find a result message
                    continue

            except Exception as e:
                # Timeout on this attempt, keep trying if we have time
                continue

    def shutdown(self):
        """Shutdown worker processes."""
        print("\nShutting down workers...")

        # Send stop signals
        for _ in range(self.n_workers):
            self.work_queue.put("STOP")

        # Wait for acknowledgments
        time.sleep(2)
        try:
            for _ in range(self.n_workers):
                msg = self.result_queue.get(timeout=1.0)
                if msg['type'] == 'shutdown':
                    print(f"  Worker {msg['worker_id']} [{msg['hostname']}] "
                          f"processed {msg['processed']} evaluations")
        except:
            pass

        print("✓ Workers shutdown complete")


def main(mode: str = "local", n_workers: int = None, problem_dim: int = 10):
    """
    Main function demonstrating Dask parallel optimization with sa_fortran.

    Args:
        mode: "local" or "slurm" cluster mode
        n_workers: Number of workers (default: 4 for local, 12 for SLURM)
        problem_dim: Dimension of the optimization problem
    """

    print("=" * 70)
    print(f"DASK PARALLEL SIMULATED ANNEALING - {mode.upper()} MODE")
    print("=" * 70)
    print(f"Problem: {problem_dim}-dimensional Rastrigin function")
    print(f"Global minimum at x = [0, 0, ..., 0] with f = 0")
    print()

    # Set default worker counts
    if n_workers is None:
        n_workers = 4 if mode == "local" else 12

    # Clear previous log file
    if PLOT_ITERATIONS:
        with open('function_evaluations_dask.csv', 'w') as f:
            pass

    # Create cluster based on mode
    cluster = None
    if mode == "local":
        client = create_local_cluster(n_workers=n_workers)
    elif mode == "slurm":
        cluster, client = create_slurm_cluster(n_workers=n_workers)
    else:
        raise ValueError(f"Invalid mode: {mode}. Use 'local' or 'slurm'")

    try:
        # Create the parallel evaluator
        evaluator = DaskParallelEvaluator(client, n_workers)

        # Define the parallel callbacks for sa_fortran
        @CALLBACK_N_INPUTS
        def n_inputs_callback(iproblem, n_inputs_ptr):
            """Tell Fortran how many inputs we can process at once (idle workers only)."""
            n_idle = evaluator.get_num_idle_workers()
            # DEBUG: Monitor idle worker count
            # if evaluator.job_counter < 20 or evaluator.job_counter % 500 == 0:
            #     print(f"DEBUG: n_inputs returning {n_idle} idle workers (pending={evaluator.pending_count})")
            n_inputs_ptr[0] = n_idle

        @CALLBACK_PARALLEL_INPUT
        def parallel_input_callback(iproblem, x_ptr, n_val, n_inputs_val):
            """Receive batch of inputs and submit to Dask worker pool."""
            total_size = n_val * n_inputs_val
            x_flat = np.ctypeslib.as_array(x_ptr, shape=(total_size,))

            # Reshape: Fortran's x(n, n_inputs) in column-major order
            x_batch = np.zeros((n_inputs_val, n_val))
            for i in range(n_inputs_val):
                x_batch[i, :] = x_flat[i * n_val:(i + 1) * n_val]

            # Submit to Dask worker pool
            evaluator.submit_batch(x_batch)

        @CALLBACK_PARALLEL_OUTPUT
        def parallel_output_callback(iproblem, x_ptr, n_val, f_ptr, istat_ptr):
            """Get one result from Dask workers - BLOCKS until result available."""
            # Block with long timeout to ensure we get a result
            # This is critical: we must wait for workers to complete before returning
            # to n_inputs_callback, ensuring workers become idle
            x_result, f, istat = evaluator.get_result(timeout=300.0)

            if istat >= 0 and x_result is not None:
                x_array = np.ctypeslib.as_array(x_ptr, shape=(n_val,))
                x_array[:] = x_result
                f_ptr[0] = f
                istat_ptr[0] = istat
            else:
                # Timeout or error - should not happen with long timeout
                print(f"WARNING: parallel_output_callback timeout or error, istat={istat}")
                f_ptr[0] = 0.0
                istat_ptr[0] = istat

        # Create and initialize the optimizer
        print("Initializing optimizer in parallel mode...")
        optimizer = sa_fortran()
        optimizer.initialize(
            n=problem_dim,
            lb=[-5.12] * problem_dim,
            ub=[5.12] * problem_dim,
            fcn=None,  # Not used in parallel mode
            c=[2.0] * problem_dim,
            maximize=False,
            eps=1e-2,
            ns=20,
            nt=100,
            neps=4,
            maxevl=10000,
            iprint=0,
            iseed1=1234,
            iseed2=5678,
            step_mode=1,  # Adaptive step size
            cooling_schedule=2,  # Fast annealing
            n_resets=5,
            optimal_f_specified=True,
            optimal_f=0.0,
            ireport=3,  # Report all evaluations and new optima
            report=report_callback,
            # Parallel callbacks
            n_inputs_to_send=n_inputs_callback,
            fcn_parallel_input=parallel_input_callback,
            fcn_parallel_output=parallel_output_callback,
        )

        print("✓ Optimizer initialized\n")
        print("=" * 70)
        print("STARTING OPTIMIZATION")
        print("=" * 70)

        # Solve
        start_time = time.time()
        sys.stdout.flush()

        result = optimizer.solve(
            x0=[-1.0] * problem_dim,
            rt=0.85,
            # t0=5.0,
            t0=10.0,
            # t0=0.0,  # monotonic test
            vm=[5.0] * problem_dim,
        )

        sys.stdout.flush()
        elapsed = time.time() - start_time

        # Print results
        print()
        print("=" * 70)
        print("FINAL RESULTS")
        print("=" * 70)
        print(f"Mode: {mode.upper()}")
        print(f"Workers: {n_workers}")
        print(f"Problem dimension: {problem_dim}")
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
        print("\n✓ Optimizer destroyed")

    finally:
        # Always cleanup cluster
        print("\nClosing cluster...")
        client.close()
        if cluster is not None:
            cluster.close()
            print("✓ SLURM cluster shut down")
        else:
            print("✓ Local cluster shut down")

    print('\nGenerating reports...')
    # Save results to CSV files (like the Fortran test does)
    import pandas as pd

    # Save all evaluations
    if all_evaluations:
        results = {f'x{i+1}': [e['x'][i] for e in all_evaluations] for i in range(problem_dim)}
        results['f'] = [e['f'] for e in all_evaluations]
        df_all = pd.DataFrame(results)
        df_all.to_csv('python_report_all.csv', index=False)
        print(f"All evaluations saved to python_report_all.csv ({len(all_evaluations)} evaluations)")
    else:
        print("No evaluations to save to python_report_all.csv")

    # Save best values
    if best_values:
        results = {f'x{i+1}': [e['x'][i] for e in best_values] for i in range(problem_dim)}
        results['f'] = [e['f'] for e in best_values]
        df_best = pd.DataFrame(results)
        df_best.to_csv('python_report_best.csv', index=False)
        print(f"Best values saved to python_report_best.csv ({len(best_values)} improvements)")
    else:
        print("No best values to save to python_report_best.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Dask parallel example for sa_fortran library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run on local cluster with default settings (4 workers, 10D problem)
  python dask_parallel_example.py

  # Run on local cluster with 8 workers
  python dask_parallel_example.py --mode local --workers 8

  # Run on local cluster with higher dimensional problem
  python dask_parallel_example.py --dimension 20

  # Run on SLURM cluster (requires dask-jobqueue)
  python dask_parallel_example.py --mode slurm --workers 12

  # SLURM with large problem
  python dask_parallel_example.py --mode slurm --workers 16 --dimension 50
        """
    )

    parser.add_argument(
        '--mode',
        type=str,
        default='local',
        choices=['local', 'slurm'],
        help='Cluster mode: "local" or "slurm" (default: local)'
    )

    parser.add_argument(
        '--workers',
        type=int,
        default=None,
        help='Number of workers (default: 4 for local, 12 for SLURM)'
    )

    parser.add_argument(
        '--dimension',
        type=int,
        default=10,
        help='Problem dimension (default: 10)'
    )

    args = parser.parse_args()

    # Run main with parsed arguments
    main(mode=args.mode, n_workers=args.workers, problem_dim=args.dimension)
