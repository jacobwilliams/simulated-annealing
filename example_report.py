"""
Example demonstrating the report callback feature.
This shows how to capture intermediate results during optimization.
"""

import numpy as np
from sa_fortran import sa_fortran

# Lists to store intermediate results
all_evaluations = []
best_values = []


def rosenbrock(x):
    """Rosenbrock function to minimize."""
    return sum(100.0 * (x[1:] - x[:-1]**2)**2 + (1 - x[:-1])**2)


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
        print(f"New best: f = {f:.6f}")


# Problem setup
n = 2
lb = np.array([-5.0, -5.0])
ub = np.array([5.0, 5.0])
x0 = np.array([2.0, -1.0])

# Create and initialize optimizer with report callback
sa = sa_fortran()
sa.initialize(
    n=n,
    lb=lb,
    ub=ub,
    fcn=rosenbrock,
    maximize=False,
    maxevl=5000,
    iprint=0,  # Disable print output to keep output clean
    ireport=3,  # Report all evaluations (1) and new optima (2)
    report=report_callback
)

# Solve the problem
result = sa.solve(x0=x0, rt=0.85, t0=5.0)

# Print summary
print(f"\nOptimization complete!")
print(f"Final solution: x = {result['x']}")
print(f"Final value: f = {result['f']:.6f}")
print(f"Function evaluations: {result['nfcnev']}")
print(f"Total evaluations captured: {len(all_evaluations)}")
print(f"New optima found: {len(best_values)}")

# Save results to CSV files (like the Fortran test does)
import pandas as pd

# Save all evaluations
if all_evaluations:
    df_all = pd.DataFrame({
        'x1': [e['x'][0] for e in all_evaluations],
        'x2': [e['x'][1] for e in all_evaluations],
        'f': [e['f'] for e in all_evaluations]
    })
    df_all.to_csv('python_report_all.csv', index=False)
    print(f"\nAll evaluations saved to python_report_all.csv")

# Save best values
if best_values:
    df_best = pd.DataFrame({
        'x1': [e['x'][0] for e in best_values],
        'x2': [e['x'][1] for e in best_values],
        'f': [e['f'] for e in best_values]
    })
    df_best.to_csv('python_report_best.csv', index=False)
    print(f"Best values saved to python_report_best.csv")

# Cleanup
sa.destroy()
