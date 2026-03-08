#!/usr/bin/env python3
"""
Example of calling the Fortran simulated annealing library from Python using ctypes.
Solves the Rosenbrock function minimization problem.
"""

import ctypes
import numpy as np
import os

# Load the shared library
lib_path = "./install/lib/libsimulated-annealing.dylib"
if not os.path.exists(lib_path):
    raise FileNotFoundError(f"Library not found: {lib_path}")

lib = ctypes.CDLL(lib_path)

# Define ctypes for the callback function signature
# void fcn(size_t ipointer, double* x, int n, double* f, int* istat)
CALLBACK_FUNC = ctypes.CFUNCTYPE(
    None,  # return type (void)
    ctypes.c_size_t,  # ipointer
    ctypes.POINTER(ctypes.c_double),  # x array
    ctypes.c_int,  # n
    ctypes.POINTER(ctypes.c_double),  # f
    ctypes.POINTER(ctypes.c_int)  # istat
)

# Define the Fortran library functions
lib.initialize_simulated_annealing.argtypes = [
    ctypes.POINTER(ctypes.c_size_t),  # ipointer (out)
    ctypes.c_int,  # n
    ctypes.POINTER(ctypes.c_double),  # lb
    ctypes.POINTER(ctypes.c_double),  # ub
    ctypes.POINTER(ctypes.c_double),  # c
    ctypes.c_bool,  # maximize
    ctypes.c_double,  # eps
    ctypes.c_int,  # ns
    ctypes.c_int,  # nt
    ctypes.c_int,  # neps
    ctypes.c_int,  # maxevl
    ctypes.c_int,  # iprint
    ctypes.c_int,  # iseed1
    ctypes.c_int,  # iseed2
    ctypes.c_int,  # step_mode
    ctypes.c_double,  # vms
    ctypes.c_int,  # iunit
    ctypes.c_bool,  # use_initial_guess
    ctypes.c_int,  # n_resets
    ctypes.c_int,  # cooling_schedule
    ctypes.c_double,  # cooling_param
    ctypes.c_bool,  # optimal_f_specified
    ctypes.c_double,  # optimal_f
    ctypes.c_double,  # optimal_f_tol
    ctypes.POINTER(ctypes.c_int),  # distribution_mode
    ctypes.POINTER(ctypes.c_double),  # dist_std_dev
    ctypes.POINTER(ctypes.c_double),  # dist_scale
    ctypes.POINTER(ctypes.c_double),  # dist_shape
    ctypes.c_void_p,  # fcn (function pointer)
    ctypes.c_void_p,  # n_inputs_to_send (NULL for serial)
    ctypes.c_void_p,  # fcn_parallel_input (NULL for serial)
    ctypes.c_void_p,  # fcn_parallel_output (NULL for serial)
]
lib.initialize_simulated_annealing.restype = None

lib.solve_simulated_annealing.argtypes = [
    ctypes.c_size_t,  # ipointer (pass by value)
    ctypes.c_int,  # n
    ctypes.POINTER(ctypes.c_double),  # x (inout)
    ctypes.c_double,  # rt
    ctypes.POINTER(ctypes.c_double),  # t (inout)
    ctypes.POINTER(ctypes.c_double),  # vm (inout)
    ctypes.POINTER(ctypes.c_double),  # xopt (out)
    ctypes.POINTER(ctypes.c_double),  # fopt (out)
    ctypes.POINTER(ctypes.c_int),  # nacc (out)
    ctypes.POINTER(ctypes.c_int),  # nfcnev (out)
    ctypes.POINTER(ctypes.c_int),  # ier (out)
]
lib.solve_simulated_annealing.restype = None

lib.destroy_simulated_annealing.argtypes = [ctypes.c_size_t]  # ipointer (pass by value)
lib.destroy_simulated_annealing.restype = None


def rosenbrock(x):
    """
    Rosenbrock function: f(x) = sum_{i=1}^{n-1} [100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2]
    Global minimum: f(1, 1, ..., 1) = 0
    """
    return sum(100.0 * (x[i+1] - x[i]**2)**2 + (1.0 - x[i])**2
               for i in range(len(x) - 1))


# Create the callback function
@CALLBACK_FUNC
def objective_function(ipointer, x_ptr, n, f_ptr, istat_ptr):
    """
    Callback function called by Fortran.

    Args:
        ipointer: integer pointer (for context, not used here)
        x_ptr: pointer to input array
        n: size of x
        f_ptr: pointer to output function value
        istat_ptr: pointer to status flag
    """
    try:
        # Convert C array to numpy array
        x = np.ctypeslib.as_array(x_ptr, shape=(n,))

        # Evaluate the function
        f = rosenbrock(x)

        # Set the output
        f_ptr[0] = f
        istat_ptr[0] = 0  # success

    except Exception as e:
        print(f"Error in objective function: {e}")
        istat_ptr[0] = -1  # error


def main():
    # Problem setup
    n = 5  # number of variables

    # Bounds
    lb = np.array([-5.0] * n, dtype=np.float64)
    ub = np.array([5.0] * n, dtype=np.float64)

    # Initial guess
    x0 = np.array([0.0] * n, dtype=np.float64)

    # Step adjustment parameter
    c = np.array([2.0] * n, dtype=np.float64)

    # Distribution mode (1 = uniform for all variables)
    distribution_mode = np.array([1] * n, dtype=np.int32)
    dist_std_dev = np.array([1.0] * n, dtype=np.float64)
    dist_scale = np.array([1.0] * n, dtype=np.float64)
    dist_shape = np.array([1.0] * n, dtype=np.float64)

    # Initialize the SA instance
    ipointer = ctypes.c_size_t()

    print(f"Solving {n}-dimensional Rosenbrock function")
    print(f"Global minimum at x = [1, 1, ..., 1] with f = 0")
    print(f"Initial guess: {x0}")
    print()

    lib.initialize_simulated_annealing(
        ctypes.byref(ipointer),  # ipointer (out)
        ctypes.c_int(n),  # n
        lb.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),  # lb
        ub.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),  # ub
        c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),  # c
        ctypes.c_bool(False),  # maximize (False = minimize)
        ctypes.c_double(1.0e-6),  # eps
        ctypes.c_int(20),  # ns
        ctypes.c_int(max(100, 5*n)),  # nt
        ctypes.c_int(4),  # neps
        ctypes.c_int(100000),  # maxevl
        ctypes.c_int(1),  # iprint (1 = summary output)
        ctypes.c_int(1234),  # iseed1
        ctypes.c_int(5678),  # iseed2
        ctypes.c_int(1),  # step_mode
        ctypes.c_double(0.1),  # vms
        ctypes.c_int(6),  # iunit (stdout)
        ctypes.c_bool(True),  # use_initial_guess
        ctypes.c_int(1),  # n_resets
        ctypes.c_int(1),  # cooling_schedule (1 = geometric)
        ctypes.c_double(1.0),  # cooling_param
        ctypes.c_bool(True),  # optimal_f_specified
        ctypes.c_double(0.0),  # optimal_f (known optimum)
        ctypes.c_double(1.0e-4),  # optimal_f_tol
        distribution_mode.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        dist_std_dev.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        dist_scale.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        dist_shape.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.cast(objective_function, ctypes.c_void_p),  # fcn
        None,  # n_inputs_to_send (NULL for serial mode)
        None,  # fcn_parallel_input (NULL)
        None,  # fcn_parallel_output (NULL)
    )

    print("SA instance initialized.")
    print()

    # Solve the optimization problem
    x = x0.copy()
    rt = 0.85  # temperature reduction factor
    t = ctypes.c_double(1.0)  # initial temperature
    vm = np.ones(n, dtype=np.float64)  # initial step lengths
    xopt = np.zeros(n, dtype=np.float64)
    fopt = ctypes.c_double(0.0)
    nacc = ctypes.c_int(0)
    nfcnev = ctypes.c_int(0)
    ier = ctypes.c_int(0)

    lib.solve_simulated_annealing(
        ipointer,
        ctypes.c_int(n),
        x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_double(rt),
        ctypes.byref(t),
        vm.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        xopt.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(fopt),
        ctypes.byref(nacc),
        ctypes.byref(nfcnev),
        ctypes.byref(ier),
    )

    # Print results
    print()
    print("=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    print(f"Exit code (ier): {ier.value}")
    print(f"  0 = success, 1 = max evals exceeded, 4 = user stop")
    print(f"Function evaluations: {nfcnev.value}")
    print(f"Accepted moves: {nacc.value}")
    print(f"Final temperature: {t.value:.6e}")
    print(f"Optimal function value: {fopt.value:.10e}")
    print(f"Optimal x: {xopt}")
    print()

    # Verify the solution
    error = np.linalg.norm(xopt - 1.0)
    print(f"Distance from true optimum [1,1,...,1]: {error:.6e}")

    # Clean up
    lib.destroy_simulated_annealing(ipointer)
    print("\nSA instance destroyed.")


if __name__ == "__main__":
    main()
