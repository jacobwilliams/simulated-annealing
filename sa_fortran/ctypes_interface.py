"""
    Main code for the Python interface to the Fortran simulated annealing library.
    This module defines the sa_fortran class, which provides methods to
    initialize the optimizer, solve optimization problems, and clean up resources.
    It uses ctypes to call the Fortran library functions and to define callback
    function signatures for both serial and parallel modes.
"""

import ctypes
import os
import numpy as np
import platform
from pathlib import Path


if platform.system() == "Windows":
    SHARED_LIB_EXTENSION = ".dll"
elif platform.system() == "Darwin":
    SHARED_LIB_EXTENSION = ".dylib"
else:
    SHARED_LIB_EXTENSION = ".so"


# Define ctypes for the callback function signature
# void fcn(size_t iproblem, double* x, int n, double* f, int* istat)
CALLBACK_FUNC = ctypes.CFUNCTYPE(
    None,  # return type (void)
    ctypes.c_size_t,  # iproblem
    ctypes.POINTER(ctypes.c_double),  # x array
    ctypes.c_int,  # n
    ctypes.POINTER(ctypes.c_double),  # f
    ctypes.POINTER(ctypes.c_int)  # istat
        )

# Define callback types for parallel mode
CALLBACK_N_INPUTS = ctypes.CFUNCTYPE(
    None,
    ctypes.c_size_t,  # iproblem
    ctypes.POINTER(ctypes.c_int)  # n_inputs
)

CALLBACK_PARALLEL_INPUT = ctypes.CFUNCTYPE(
    None,
    ctypes.c_size_t,  # iproblem
    ctypes.POINTER(ctypes.c_double),  # x array
    ctypes.c_int,  # n
    ctypes.c_int  # n_inputs
)

CALLBACK_PARALLEL_OUTPUT = ctypes.CFUNCTYPE(
    None,
    ctypes.c_size_t,  # iproblem
    ctypes.POINTER(ctypes.c_double),  # x
    ctypes.c_int,  # n
    ctypes.POINTER(ctypes.c_double),  # f
    ctypes.POINTER(ctypes.c_int)  # istat
)


class sa_fortran():
    """
    This module provides an interface to the Fortran simulated annealing library.
    It allows users to call the Fortran library from Python using ctypes.
    The Fortran library is built as a shared library and can be used in Python by loading it with ctypes.
    The example provided demonstrates how to use the library to minimize the Rosenbrock function.
    """

    def __init__(self):

        # Load the shared library
        lib_path = str(Path(__file__).parent / 'lib' / 'lib' /
                       f'libsimulated-annealing{SHARED_LIB_EXTENSION}')
        if not os.path.exists(lib_path):
            raise FileNotFoundError(f"Library not found: {lib_path}")

        self.lib = ctypes.CDLL(lib_path)

        # Define the Fortran library functions
        self.lib.initialize_simulated_annealing.argtypes = [
            ctypes.POINTER(ctypes.c_size_t),  # iproblem (out)
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
        self.lib.initialize_simulated_annealing.restype = None

        self.lib.solve_simulated_annealing.argtypes = [
            ctypes.c_size_t,  # iproblem (pass by value)
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
        self.lib.solve_simulated_annealing.restype = None

        self.lib.destroy_simulated_annealing.argtypes = [ctypes.c_size_t]  # iproblem (pass by value)
        self.lib.destroy_simulated_annealing.restype = None

        # Instance variables
        self.iproblem = None
        self.callback_ref = None  # Keep reference to prevent garbage collection
        self.n = None

    def initialize(self, n, lb, ub, fcn=None, c=None, maximize=False, eps=1e-6,
                   ns=20, nt=None, neps=4, maxevl=100000, iprint=1,
                   iseed1=1234, iseed2=5678, step_mode=1, vms=0.1, iunit=6,
                   use_initial_guess=True, n_resets=1, cooling_schedule=1,
                   cooling_param=1.0, optimal_f_specified=False, optimal_f=0.0,
                   optimal_f_tol=1e-4, distribution_mode=None, dist_std_dev=None,
                   dist_scale=None, dist_shape=None,
                   n_inputs_to_send=None, fcn_parallel_input=None, fcn_parallel_output=None):
        """
        Initialize the simulated annealing optimizer.

        Parameters
        ----------
        n : int
            Number of variables
        lb : array-like
            Lower bounds (size n)
        ub : array-like
            Upper bounds (size n)
        fcn : callable, optional
            Objective function with signature: fcn(x) -> f
            where x is array of size n, returns scalar f.
            Required for serial mode, should be None for parallel mode.
        c : array-like, optional
            Step adjustment parameter (default: 2.0 for all variables)
        maximize : bool, optional
            True to maximize, False to minimize (default: False)
        eps : float, optional
            Convergence tolerance (default: 1e-6)
        ns : int, optional
            Number of cycles before step adjustment (default: 20)
        nt : int, optional
            Number of iterations before temperature reduction (default: max(100, 5*n))
        neps : int, optional
            Number of final values for termination check (default: 4)
        maxevl : int, optional
            Maximum number of function evaluations (default: 100000)
        iprint : int, optional
            Print level: 0=none, 1=summary, 2=detailed, 3=verbose (default: 1)
        iseed1, iseed2 : int, optional
            Random number seeds (default: 1234, 5678)
        step_mode : int, optional
            Step adjustment mode: 1=adaptive, 2=constant, 3=factor (default: 1)
        vms : float, optional
            Step adjustment factor for step_mode=3 (default: 0.1)
        iunit : int, optional
            Output unit number (default: 6 for stdout)
        use_initial_guess : bool, optional
            Use initial guess in solve() (default: True)
        n_resets : int, optional
            Number of main loop resets (default: 1)
        cooling_schedule : int, optional
            Temperature schedule: 1=geometric, 2=fast, 3=huang, 4=boltzmann, 5=logarithmic (default: 1)
        cooling_param : float, optional
            Parameter for cooling schedules 3 and 5 (default: 1.0)
        optimal_f_specified : bool, optional
            Whether optimal value is known (default: False)
        optimal_f : float, optional
            Known optimal value (default: 0.0)
        optimal_f_tol : float, optional
            Tolerance for optimal_f check (default: 1e-4)
        distribution_mode : array-like, optional
            Distribution for perturbations per variable: 1=uniform, 2=normal, 3=cauchy, 4=triangular, 5=bipareto (default: 1 for all)
        dist_std_dev : array-like, optional
            Standard deviation for normal distribution (default: 1.0 for all)
        dist_scale : array-like, optional
            Scale for cauchy/bipareto (default: 1.0 for all)
        dist_shape : array-like, optional
            Shape for bipareto (default: 1.0 for all)
        n_inputs_to_send : callable, optional
            Callback for parallel mode: n_inputs_to_send(iproblem) -> n_inputs
            Returns number of inputs ready to evaluate. Required for parallel mode.
        fcn_parallel_input : callable, optional
            Callback for parallel mode: fcn_parallel_input(iproblem, x, n, n_inputs)
            Receives batch of inputs to evaluate. Required for parallel mode.
        fcn_parallel_output : callable, optional
            Callback for parallel mode: fcn_parallel_output(iproblem, x, n, f, istat)
            Returns one completed result. Required for parallel mode.
        """
        self.n = n

        # Convert arrays to numpy and ensure correct dtype
        lb = np.asarray(lb, dtype=np.float64)
        ub = np.asarray(ub, dtype=np.float64)

        if c is None:
            c = np.full(n, 2.0, dtype=np.float64)
        else:
            c = np.asarray(c, dtype=np.float64)

        if nt is None:
            nt = max(100, 5 * n)

        # Distribution parameters
        if distribution_mode is None:
            distribution_mode = np.ones(n, dtype=np.int32)
        else:
            distribution_mode = np.asarray(distribution_mode, dtype=np.int32)

        if dist_std_dev is None:
            dist_std_dev = np.ones(n, dtype=np.float64)
        else:
            dist_std_dev = np.asarray(dist_std_dev, dtype=np.float64)

        if dist_scale is None:
            dist_scale = np.ones(n, dtype=np.float64)
        else:
            dist_scale = np.asarray(dist_scale, dtype=np.float64)

        if dist_shape is None:
            dist_shape = np.ones(n, dtype=np.float64)
        else:
            dist_shape = np.asarray(dist_shape, dtype=np.float64)

        # Determine mode: serial or parallel
        is_parallel = (n_inputs_to_send is not None and
                      fcn_parallel_input is not None and
                      fcn_parallel_output is not None)

        if is_parallel:
            # Parallel mode: use provided callbacks directly
            if fcn is not None:
                print("Warning: fcn is ignored in parallel mode")

            # Keep references to prevent garbage collection
            self.callback_ref = (n_inputs_to_send, fcn_parallel_input, fcn_parallel_output)

            fcn_ptr = ctypes.c_void_p(0)
            n_inputs_ptr = ctypes.cast(n_inputs_to_send, ctypes.c_void_p)
            parallel_input_ptr = ctypes.cast(fcn_parallel_input, ctypes.c_void_p)
            parallel_output_ptr = ctypes.cast(fcn_parallel_output, ctypes.c_void_p)

        else:
            # Serial mode: create callback wrapper for fcn
            if fcn is None:
                raise ValueError("fcn is required for serial mode")

            @CALLBACK_FUNC
            def callback_wrapper(iproblem, x_ptr, n_val, f_ptr, istat_ptr):
                try:
                    x = np.ctypeslib.as_array(x_ptr, shape=(n_val,))
                    f = fcn(x)
                    f_ptr[0] = float(f)
                    istat_ptr[0] = 0
                except Exception as e:
                    print(f"Error in objective function: {e}")
                    istat_ptr[0] = -1

            # Keep reference to prevent garbage collection
            self.callback_ref = callback_wrapper

            fcn_ptr = ctypes.cast(callback_wrapper, ctypes.c_void_p)
            n_inputs_ptr = ctypes.c_void_p(0)
            parallel_input_ptr = ctypes.c_void_p(0)
            parallel_output_ptr = ctypes.c_void_p(0)

        # Initialize iproblem
        self.iproblem = ctypes.c_size_t()

        # Call Fortran initialize
        self.lib.initialize_simulated_annealing(
            ctypes.byref(self.iproblem),
            ctypes.c_int(n),
            lb.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ub.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_bool(maximize),
            ctypes.c_double(eps),
            ctypes.c_int(ns),
            ctypes.c_int(nt),
            ctypes.c_int(neps),
            ctypes.c_int(maxevl),
            ctypes.c_int(iprint),
            ctypes.c_int(iseed1),
            ctypes.c_int(iseed2),
            ctypes.c_int(step_mode),
            ctypes.c_double(vms),
            ctypes.c_int(iunit),
            ctypes.c_bool(use_initial_guess),
            ctypes.c_int(n_resets),
            ctypes.c_int(cooling_schedule),
            ctypes.c_double(cooling_param),
            ctypes.c_bool(optimal_f_specified),
            ctypes.c_double(optimal_f),
            ctypes.c_double(optimal_f_tol),
            distribution_mode.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            dist_std_dev.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            dist_scale.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            dist_shape.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            fcn_ptr,
            n_inputs_ptr,
            parallel_input_ptr,
            parallel_output_ptr,
        )

    def destroy(self):
        """
        Destroy the simulated annealing instance and free memory.
        """
        if self.iproblem is not None:
            self.lib.destroy_simulated_annealing(self.iproblem)
            self.iproblem = None
            self.callback_ref = None
            self.n = None

    def solve(self, x0, rt=0.85, t0=1.0, vm=None):
        """
        Solve the optimization problem.

        Parameters
        ----------
        x0 : array-like
            Initial guess (size n)
        rt : float, optional
            Temperature reduction factor (default: 0.85)
        t0 : float, optional
            Initial temperature (default: 1.0)
        vm : array-like, optional
            Initial step lengths (default: ones(n))

        Returns
        -------
        result : dict
            Dictionary containing:
            - 'x': array, optimal solution
            - 'f': float, optimal function value
            - 'xopt': array, optimal solution (same as 'x')
            - 'fopt': float, optimal function value (same as 'f')
            - 'nacc': int, number of accepted moves
            - 'nfcnev': int, number of function evaluations
            - 'ier': int, exit code (0=success, 1=maxevl exceeded, 4=user stop)
            - 't': float, final temperature
            - 'vm': array, final step lengths
        """
        if self.iproblem is None:
            raise RuntimeError("Must call initialize() before solve()")

        # Convert inputs to numpy arrays
        x = np.asarray(x0, dtype=np.float64).copy()

        if vm is None:
            vm = np.ones(self.n, dtype=np.float64)
        else:
            vm = np.asarray(vm, dtype=np.float64).copy()

        # Prepare output variables
        t = ctypes.c_double(t0)
        xopt = np.zeros(self.n, dtype=np.float64)
        fopt = ctypes.c_double(0.0)
        nacc = ctypes.c_int(0)
        nfcnev = ctypes.c_int(0)
        ier = ctypes.c_int(0)

        # Call Fortran solve
        self.lib.solve_simulated_annealing(
            self.iproblem,
            ctypes.c_int(self.n),
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

        # Return results as dictionary
        return {
            'x': xopt,
            'f': fopt.value,
            'xopt': xopt,
            'fopt': fopt.value,
            'nacc': nacc.value,
            'nfcnev': nfcnev.value,
            'ier': ier.value,
            't': t.value,
            'vm': vm,
        }

    def __del__(self):
        """Destructor to ensure cleanup."""
        self.destroy()

