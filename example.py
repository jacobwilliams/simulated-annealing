#!/usr/bin/env python
"""
Example using the Pythonic sa_fortran class interface.
Solves the Rosenbrock function minimization problem.

Note: the shared library must be built first using:
`fpm install --prefix ./sa_fortran/lib --profile release`
"""

import numpy as np
from sa_fortran import sa_fortran


def rosenbrock(x):
    """
    Rosenbrock function: f(x) = sum_{i=1}^{n-1} [100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2]
    Global minimum: f(1, 1, ..., 1) = 0
    """
    return sum(100.0 * (x[i+1] - x[i]**2)**2 + (1.0 - x[i])**2
               for i in range(len(x) - 1))


def main():
    # Problem setup
    n = 5
    lb = [-5.0] * n
    ub = [5.0] * n
    x0 = [0.0] * n

    print(f"Solving {n}-dimensional Rosenbrock function")
    print(f"Global minimum at x = [1, 1, ..., 1] with f = 0")
    print(f"Initial guess: {x0}")
    print()

    # Create optimizer instance
    optimizer = sa_fortran()

    # Initialize with problem setup
    optimizer.initialize(
        n=n,
        lb=lb,
        ub=ub,
        fcn=rosenbrock,
        maximize=False,  # minimize
        eps=1e-6,
        ns=20,
        nt=max(100, 5*n),
        neps=4,
        maxevl=100000,
        iprint=1,  # summary output
        iseed1=1234,
        iseed2=5678,
        cooling_schedule=1,  # geometric
        optimal_f_specified=True,
        optimal_f=0.0,
        optimal_f_tol=1e-4,
    )

    print("Optimizer initialized.")
    print()

    # Solve the optimization problem
    result = optimizer.solve(
        x0=x0,
        rt=0.85,  # temperature reduction factor
        t0=1.0,   # initial temperature
    )

    # Print results
    print()
    print("=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    print(f"Exit code (ier): {result['ier']}")
    print(f"  0 = success, 1 = max evals exceeded, 4 = user stop")
    print(f"Function evaluations: {result['nfcnev']}")
    print(f"Accepted moves: {result['nacc']}")
    print(f"Final temperature: {result['t']:.6e}")
    print(f"Optimal function value: {result['f']:.10e}")
    print(f"Optimal x: {result['x']}")
    print()

    # Verify the solution
    error = np.linalg.norm(result['x'] - 1.0)
    print(f"Distance from true optimum [1,1,...,1]: {error:.6e}")

    # Clean up
    optimizer.destroy()
    print("\nOptimizer destroyed.")


if __name__ == "__main__":
    main()
