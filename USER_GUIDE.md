# Simulated Annealing User Guide

A comprehensive guide to using the modern Fortran simulated annealing optimization library.

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Basic Usage](#basic-usage)
5. [Configuration Parameters](#configuration-parameters)
6. [Advanced Features](#advanced-features)
7. [Parallel Evaluation](#parallel-evaluation)
8. [Best Practices](#best-practices)
9. [Examples](#examples)
10. [Troubleshooting](#troubleshooting)

---

## Introduction

### What is Simulated Annealing?

Simulated annealing is a probabilistic optimization technique inspired by the annealing process in metallurgy. It's particularly effective for finding global optima of complex, multimodal functions where traditional gradient-based methods might get stuck in local minima.

**Key characteristics:**
- **Global optimization**: Can escape local minima to find global optimum
- **Derivative-free**: No gradient information needed
- **Robust**: Works well with noisy, discontinuous, or non-smooth functions
- **Flexible**: Minimal assumptions about the objective function

**How it works:**
1. Start from an initial point and temperature
2. Randomly perturb the current solution
3. Accept improvements (downhill moves for minimization)
4. Probabilistically accept worse solutions (uphill moves) based on temperature
5. Gradually reduce temperature, focusing search on promising regions
6. Terminate when convergence criteria are met

### When to Use Simulated Annealing

**Good for:**
- Multimodal functions with many local optima
- Black-box optimization (no derivatives available)
- Discontinuous or non-smooth objective functions
- Problems where finding a near-optimal solution is acceptable
- Constrained optimization within bounded regions

**Not ideal for:**
- High-dimensional problems (>50-100 variables) without parallelization
- Problems requiring guaranteed global optimum
- Real-time optimization (can be computationally expensive)
- Problems where gradient information is readily available and function is convex

---

## Installation

### Prerequisites

- **Fortran compiler**: gfortran 8.0+ or ifort
- **FPM** (Fortran Package Manager): [Installation instructions](https://fpm.fortran-lang.org/install/index.html)
- **Python 3.7+** (optional, for Python interface): numpy, ctypes

### Building the Library

#### Using FPM (Recommended)

```bash
# Clone the repository
git clone https://github.com/jacobwilliams/simulated-annealing.git
cd simulated-annealing

# Build the library
fpm build --profile release

# Run tests to verify installation
fpm test --profile release
```

#### For Python Interface

Build and install the shared library:

```bash
fpm install --prefix ./sa_fortran/lib --profile release
```

This creates the shared library needed by the Python interface in `sa_fortran/lib/`.

#### As FPM Dependency

Add to your `fpm.toml`:

```toml
[dependencies]
simulated-annealing = { git="https://github.com/jacobwilliams/simulated-annealing.git" }
```

---

## Quick Start

### Python Quick Start

```python
import numpy as np
from sa_fortran import sa_fortran

# Define objective function
def sphere(x):
    return sum(xi**2 for xi in x)

# Setup
n = 5
optimizer = sa_fortran()
optimizer.initialize(
    n=n,
    lb=[-10.0] * n,
    ub=[10.0] * n,
    fcn=sphere,
    maximize=False,
    maxevl=10000
)

# Optimize
result = optimizer.solve(
    x0=[5.0] * n,
    rt=0.85,
    t0=1.0
)

print(f"Optimal f: {result['f']}")
print(f"Optimal x: {result['x']}")

optimizer.destroy()
```

### Fortran Quick Start

```fortran
program simple_example
   use simulated_annealing_module
   implicit none

   type(simulated_annealing_type) :: solver
   real(wp), dimension(2) :: x, xopt, lb, ub, vm
   real(wp) :: fopt, rt, t
   integer :: nacc, nfcnev, ier

   ! Setup problem
   lb = [-5.0_wp, -5.0_wp]
   ub = [5.0_wp, 5.0_wp]
   x = [0.0_wp, 0.0_wp]
   vm = [1.0_wp, 1.0_wp]

   ! Initialize solver
   call solver%initialize(n=2, lb=lb, ub=ub, fcn=rosenbrock, maximize=.false.)

   ! Solve
   call solver%optimize(x, rt=0.85_wp, t=1.0_wp, vm=vm, &
                        xopt=xopt, fopt=fopt, nacc=nacc, nfcnev=nfcnev, ier=ier)

   print *, 'Optimal f:', fopt
   print *, 'Optimal x:', xopt

contains

   subroutine rosenbrock(me, x, f, istat)
      class(simulated_annealing_type), intent(inout) :: me
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: f
      integer, intent(out) :: istat

      f = 100.0_wp * (x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2
      istat = 0  ! success
   end subroutine rosenbrock

end program simple_example
```

---

## Basic Usage

### Defining Your Objective Function

#### Python

```python
def my_function(x):
    """
    Args:
        x: numpy array of decision variables

    Returns:
        float: objective function value
    """
    # Your calculation here
    return result
```

The function should:
- Accept a numpy array (or list) of variables
- Return a single float value
- Raise an exception or return inf/nan for invalid inputs

#### Fortran

```fortran
subroutine my_function(me, x, f, istat)
   class(simulated_annealing_type), intent(inout) :: me
   real(wp), dimension(:), intent(in) :: x
   real(wp), intent(out) :: f
   integer, intent(out) :: istat

   ! Calculate objective function
   f = ! your calculation

   ! Set status
   istat = 0   ! success
   ! istat = -1  ! invalid point, try another
   ! istat = -2  ! stop optimization
end subroutine my_function
```

### Setting Up the Optimizer

#### Python Example

```python
optimizer = sa_fortran()

optimizer.initialize(
    n=5,                    # number of variables
    lb=[-10.0] * 5,        # lower bounds
    ub=[10.0] * 5,         # upper bounds
    fcn=my_function,       # objective function
    maximize=False,         # minimize (True for maximize)
    eps=1e-6,              # convergence tolerance
    ns=20,                 # cycles per temperature
    nt=100,                # iterations before cooling
    maxevl=100000,         # max function evaluations
    iprint=1,              # output level (0-3)
    cooling_schedule=1,    # cooling method
)
```

#### Running the Optimization

```python
result = optimizer.solve(
    x0=[1.0] * 5,   # initial guess
    rt=0.85,        # cooling rate (0 < rt < 1)
    t0=1.0,         # initial temperature
    vm=[2.0] * 5    # initial step sizes (optional)
)
```

#### Understanding Results

```python
result = {
    'x': optimized_variables,    # optimal solution
    'f': optimal_value,          # optimal objective value
    't': final_temperature,       # final temperature
    'vm': final_step_sizes,       # final step sizes
    'ier': exit_code,            # 0=success, 1=maxevl, 4=stop
    'nfcnev': n_evaluations,     # total function calls
    'nacc': n_accepted           # accepted moves
}
```

---

## Configuration Parameters

### Essential Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n` | int | required | Number of optimization variables |
| `lb` | array[n] | required | Lower bounds for each variable |
| `ub` | array[n] | required | Upper bounds for each variable |
| `fcn` | function | required | Objective function to optimize |
| `maximize` | bool | False | True for maximization, False for minimization |

### Convergence Control

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `eps` | float | 1e-9 | Convergence tolerance |
| `neps` | int | 4 | Number of temperatures for convergence check |
| `maxevl` | int | 10000 | Maximum function evaluations |
| `optimal_f_specified` | bool | False | Whether optimal value is known |
| `optimal_f` | float | 0.0 | Known optimal value (if specified) |
| `optimal_f_tol` | float | 0.0 | Tolerance for optimal_f check |

### Temperature Schedule

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cooling_schedule` | int | 1 | Cooling method (1-5) |
| `cooling_param` | float | 1.0 | Parameter for schedules 3, 5 |
| `cooling_exponent` | float | 1.0 | Exponent for schedule 3 |

**Cooling schedules:**
1. **Geometric** (default): T(k+1) = rt × T(k)
2. **Fast annealing** (Cauchy): T(k) = T₀ / (1 + k)
3. **Huang**: T(k) = T₀ / (1 + param × k)^exponent
4. **Boltzmann**: T(k) = T₀ / log(1 + k + e)
5. **Logarithmic**: T(k) = T₀ / (1 + param × log(1 + k))

### Algorithm Tuning

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ns` | int | 20 | Cycles before VM adjustment |
| `nt` | int | 100 | Iterations before temperature reduction |
| `c` | array[n] | [2.0] | Step size adjustment factor per variable |
| `step_mode` | int | 1 | VM adjustment method (1-3) |
| `vms` | float | 0.1 | Factor for step_mode=3 |
| `n_resets` | int | 2 | Number of restarts with reset conditions |
| `use_initial_guess` | bool | True | Use x0 or random start |

**Step modes:**
1. **Adaptive** (default): Adjust to ~50% acceptance
2. **Constant**: Keep VM fixed
3. **Factor**: Multiply VM by constant factor

### Random Number Control

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `iseed1` | int | 1234 | First random seed |
| `iseed2` | int | 5678 | Second random seed |

Different seeds produce different optimization paths.

### Output Control

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `iprint` | int | 1 | Console output level (0-3) |
| `ireport` | int | 0 | Callback reporting frequency (0-5) |
| `report` | function | None | Custom reporting callback |

**Output levels (iprint):**
- 0: Silent
- 1: Summary per temperature (recommended)
- 2: Step length adjustments
- 3: Every function evaluation (verbose)

**Report levels (ireport):**
- 0: No callbacks
- 1: Each valid evaluation
- 2: Each new optimum
- 3: Valid evaluations + new optima
- 4: All evaluations (including invalid)
- 5: All evaluations + new optima

---

## Advanced Features

### Custom Perturbation Distributions

Control how variables are perturbed at each step:

```python
optimizer.initialize(
    n=3,
    lb=[-5.0, -5.0, -5.0],
    ub=[5.0, 5.0, 5.0],
    fcn=my_function,
    distribution_mode=[
        1,  # uniform (default)
        2,  # normal/Gaussian
        3   # Cauchy (heavy-tailed)
    ],
    dist_std_dev=[1.0, 0.5, 0.5],  # for normal
    dist_scale=[1.0, 1.0, 2.0],     # for Cauchy
)
```

**Available distributions:**
- `0`: Constant (no perturbation)
- `1`: Uniform (default)
- `2`: Normal/Gaussian
- `3`: Cauchy (heavy-tailed, good for occasional large jumps)
- `4`: Triangular
- `5`: Bipareto (two-sided Pareto)

### Intermediate Result Reporting

Monitor optimization progress in real-time:

```python
def report_callback(x, f, istat):
    """
    Called during optimization.

    Args:
        x: current point
        f: function value
        istat: 1=function eval, 2=new optimum, 3=invalid eval
    """
    if istat == 2:
        print(f"New best: f={f:.6f}")

optimizer.initialize(
    ...,
    ireport=3,  # report all evals and new optima
    report=report_callback
)
```

### Early Stopping

Stop when known solution is found:

```python
optimizer.initialize(
    ...,
    optimal_f_specified=True,
    optimal_f=0.0,
    optimal_f_tol=1e-4  # stop when |f - optimal_f| < tol
)
```

### Variable-Specific Configuration

Different settings per variable:

```python
optimizer.initialize(
    n=3,
    lb=[-10.0, -1.0, -100.0],
    ub=[10.0, 1.0, 100.0],
    c=[2.0, 1.5, 3.0],  # different adjustment rates
    distribution_mode=[1, 2, 3],  # different distributions
    dist_std_dev=[1.0, 0.2, 2.0],  # different scales
    fcn=my_function
)
```

---

## Parallel Evaluation

For expensive objective functions, evaluate multiple points simultaneously using parallel workers (GPU, OpenMP, MPI, Dask, etc.).

### Architecture

```
Main Process:                  Worker Pool:
  SA Algorithm  ──────┬────────▶ Worker 1: evaluate f(x₁)
                      ├────────▶ Worker 2: evaluate f(x₂)
                      ├────────▶ Worker 3: evaluate f(x₃)
                      └────────▶ Worker 4: evaluate f(x₄)
      ▲
      └─────────────────────────── Results: [f₁, f₂, f₃, f₄]
```

### Python Parallel Example with Dask

See `dask_parallel_example.py` for a complete implementation:

```python
from dask.distributed import Client, Queue
from sa_fortran import CALLBACK_N_INPUTS, CALLBACK_PARALLEL_INPUT, CALLBACK_PARALLEL_OUTPUT

# Define callbacks
@CALLBACK_N_INPUTS
def n_inputs_callback(iproblem, n_inputs_ptr):
    """Tell SA how many points we can evaluate in parallel."""
    n_inputs_ptr[0] = n_idle_workers

@CALLBACK_PARALLEL_INPUT
def parallel_input_callback(iproblem, x_ptr, n_val, n_inputs_val):
    """Receive batch of points and submit to workers."""
    # Send to worker pool for evaluation
    for i in range(n_inputs_val):
        x = extract_x_from_array(x_ptr, i, n_val)
        work_queue.put({'x': x})

@CALLBACK_PARALLEL_OUTPUT
def parallel_output_callback(iproblem, x_ptr, n_val, f_ptr, istat_ptr):
    """Get one completed result from workers."""
    result = result_queue.get()  # blocking
    x_ptr[:] = result['x']
    f_ptr[0] = result['f']
    istat_ptr[0] = result['istat']

# Initialize optimizer with parallel callbacks
optimizer.initialize(
    n=n,
    lb=lb,
    ub=ub,
    fcn=None,  # not used in parallel mode
    n_inputs_to_send=n_inputs_callback,
    fcn_parallel_input=parallel_input_callback,
    fcn_parallel_output=parallel_output_callback,
)
```

### Running Parallel Example

```bash
# Local cluster with 8 workers
python dask_parallel_example.py --workers 8 --dimension 10

# SLURM cluster (HPC)
python dask_parallel_example.py --mode slurm --workers 48 --dimension 20
```

---

## Best Practices

### Choosing Initial Temperature

The temperature controls exploration vs. exploitation:

**Method 1: Trial runs**
```python
# Run with rt > 1 to increase temperature
result = optimizer.solve(x0, rt=1.5, t0=1.0)
# Observe VM (step sizes) in output
# Choose T that gives appropriate VM
```

**Method 2: Acceptance rate test**
- Start with T where ~50-80% of moves are accepted
- For typical functions: T ≈ Δf_typical (typical function change)

**Guidelines:**
- Higher T: More exploration, slower convergence
- Lower T: Less exploration, faster convergence (risk of local minima)
- Typical range: 0.1 to 10.0 for normalized problems

### Tuning Cooling Rate (rt)

Slower cooling = better quality, more time:

| rt | Speed | Quality | Use Case |
|----|-------|---------|----------|
| 0.95-0.99 | Slow | Best | Research, critical applications |
| 0.85-0.9 | Medium | Good | Production (recommended) |
| 0.7-0.8 | Fast | Fair | Quick exploration, testing |

### Setting NS and NT

Controls when adjustments happen:

```python
ns = 20         # cycles before VM adjustment (good default)
nt = max(100, 5*n)  # iterations before cooling (Corana et al.)
```

For n-dimensional problems:
- Small n (< 10): ns=20, nt=100
- Medium n (10-50): ns=20, nt=5*n
- Large n (> 50): Consider parallel evaluation

### Bounds Selection

```python
# Strategy 1: Physical bounds (if known)
lb = [0.0, -180.0, 1e-6]  # physical limits
ub = [100.0, 180.0, 1.0]

# Strategy 2: Exploratory bounds
range_estimate = 10.0
lb = [x0_i - range_estimate for x0_i in x0]
ub = [x0_i + range_estimate for x0_i in x0]

# Strategy 3: Very wide bounds (algorithm will focus search)
lb = [-1e6] * n
ub = [1e6] * n
```

### Function Evaluation Tips

```python
def robust_function(x):
    """Well-designed objective function."""

    # 1. Validate inputs
    if np.any(np.isnan(x)) or np.any(np.isinf(x)):
        return np.inf  # or raise ValueError

    # 2. Try-catch for numerical issues
    try:
        result = expensive_calculation(x)
    except (ValueError, RuntimeError):
        return np.inf  # penalty for invalid region

    # 3. Check output validity
    if np.isnan(result) or np.isinf(result):
        return 1e30  # large penalty

    return result
```

### Multi-start Strategy

For difficult problems, use multiple restarts:

```python
optimizer.initialize(
    ...,
    n_resets=5,  # run 5 times with different starting conditions
    iseed1=1234,  # change seeds for different runs
)
```

Or manually:
```python
best_f = np.inf
best_x = None

for run in range(10):
    optimizer = sa_fortran()
    optimizer.initialize(..., iseed1=1000+run, iseed2=2000+run)
    result = optimizer.solve(...)

    if result['f'] < best_f:
        best_f = result['f']
        best_x = result['x']

    optimizer.destroy()
```

### Convergence Monitoring

```python
evaluations = []

def monitor(x, f, istat):
    if istat == 1:  # each evaluation
        evaluations.append(f)

optimizer.initialize(..., ireport=1, report=monitor)
result = optimizer.solve(...)

# Plot convergence
import matplotlib.pyplot as plt
plt.plot(evaluations)
plt.xlabel('Evaluation')
plt.ylabel('Objective Value')
plt.yscale('log')
plt.show()
```

---

## Examples

### Example 1: Rosenbrock Function (5D)

Classic test problem with narrow valley:

```python
import numpy as np
from sa_fortran import sa_fortran

def rosenbrock(x):
    return sum(100.0 * (x[i+1] - x[i]**2)**2 + (1.0 - x[i])**2
               for i in range(len(x) - 1))

n = 5
optimizer = sa_fortran()
optimizer.initialize(
    n=n,
    lb=[-5.0] * n,
    ub=[5.0] * n,
    fcn=rosenbrock,
    maximize=False,
    ns=20,
    nt=max(100, 5*n),
    maxevl=100000,
    cooling_schedule=1,  # geometric
    optimal_f_specified=True,
    optimal_f=0.0,
    optimal_f_tol=1e-4,
)

result = optimizer.solve(
    x0=[0.0] * n,
    rt=0.85,
    t0=1.0,
)

print(f"Found: f={result['f']:.6e}, x={result['x']}")
print(f"Error from [1,1,...,1]: {np.linalg.norm(result['x'] - 1.0):.6e}")

optimizer.destroy()
```

### Example 2: Rastrigin Function (Multimodal)

Highly multimodal with many local minima:

```python
def rastrigin(x):
    A = 10.0
    n = len(x)
    return A * n + sum(xi**2 - A * np.cos(2 * np.pi * xi) for xi in x)

n = 10
optimizer = sa_fortran()
optimizer.initialize(
    n=n,
    lb=[-5.12] * n,
    ub=[5.12] * n,
    fcn=rastrigin,
    maximize=False,
    ns=20,
    nt=max(100, 5*n),
    maxevl=200000,
    cooling_schedule=2,  # fast annealing (good for multimodal)
    optimal_f_specified=True,
    optimal_f=0.0,
)

result = optimizer.solve(
    x0=[-2.0] * n,
    rt=0.9,  # slower cooling for difficult function
    t0=5.0,  # higher initial temperature
)

print(f"Global minimum: f={result['f']:.6e}")
optimizer.destroy()
```

### Example 3: Constrained Optimization with Penalty

Minimize x² + y² subject to x + y ≥ 1:

```python
def constrained_objective(x):
    # Objective
    obj = x[0]**2 + x[1]**2

    # Constraint: x + y >= 1
    constraint_violation = max(0, 1 - (x[0] + x[1]))

    # Penalty method
    penalty = 1000.0 * constraint_violation**2

    return obj + penalty

optimizer = sa_fortran()
optimizer.initialize(
    n=2,
    lb=[-10.0, -10.0],
    ub=[10.0, 10.0],
    fcn=constrained_objective,
    maximize=False,
    maxevl=50000,
)

result = optimizer.solve(
    x0=[0.0, 0.0],
    rt=0.85,
    t0=1.0,
)

print(f"Optimal: f={result['f']:.6f}, x={result['x']}")
print(f"Constraint x+y={sum(result['x']):.6f} (should be ≥ 1)")

optimizer.destroy()
```

### Example 4: Parameter Fitting

Fit model to noisy data:

```python
import numpy as np

# Generate noisy data
def true_model(x, a, b, c):
    return a * np.exp(-b * x) + c

x_data = np.linspace(0, 5, 50)
y_true = true_model(x_data, 2.0, 0.5, 1.0)
y_data = y_true + np.random.normal(0, 0.1, len(x_data))

# Optimization
def residual(params):
    a, b, c = params
    if b < 0:  # enforce physical constraint
        return 1e10
    y_model = a * np.exp(-b * x_data) + c
    return np.sum((y_data - y_model)**2)

optimizer = sa_fortran()
optimizer.initialize(
    n=3,
    lb=[0.0, 0.0, 0.0],
    ub=[10.0, 5.0, 5.0],
    fcn=residual,
    maximize=False,
    maxevl=50000,
)

result = optimizer.solve(
    x0=[1.0, 1.0, 1.0],
    rt=0.85,
    t0=0.5,
)

a_fit, b_fit, c_fit = result['x']
print(f"Fitted: a={a_fit:.3f}, b={b_fit:.3f}, c={c_fit:.3f}")
print(f"True:   a=2.0, b=0.5, c=1.0")
print(f"SSE: {result['f']:.6f}")

optimizer.destroy()
```

### Example 5: Parallel Evaluation (Simplified)

```python
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from queue import Queue

def expensive_function(x):
    """Simulated expensive calculation."""
    import time
    time.sleep(0.01)  # simulate computation
    return sum(xi**2 for xi in x)

# Worker pool
executor = ProcessPoolExecutor(max_workers=4)
work_queue = Queue()
result_queue = Queue()

def submit_work(x_batch):
    """Submit batch to workers."""
    for x in x_batch:
        future = executor.submit(expensive_function, x)
        work_queue.put(future)

def get_result():
    """Get completed result."""
    future = work_queue.get()
    f = future.result()
    return f

# Use with parallel callbacks...
# (See dask_parallel_example.py for complete implementation)
```

---

## Troubleshooting

### Problem: Algorithm Not Finding Global Minimum

**Solutions:**
1. Increase initial temperature: `t0=5.0` or higher
2. Slow down cooling: `rt=0.9` instead of 0.85
3. Increase iterations: larger `maxevl`, `nt`
4. Use multi-start: `n_resets=5` or multiple independent runs
5. Try different cooling schedule: `cooling_schedule=2` (fast annealing)

### Problem: Optimization Too Slow

**Solutions:**
1. Faster cooling: `rt=0.8` instead of 0.85
2. Reduce iterations: smaller `nt`, `maxevl`
3. Use parallel evaluation (see Parallel Section)
4. Tighten bounds if possible
5. Fast cooling schedule: `cooling_schedule=2`

### Problem: Poor Convergence

**Solutions:**
1. Check function for bugs (test with known inputs)
2. Scale variables to similar ranges
3. Adjust `ns` and `nt`: `nt = max(100, 5*n)`
4. Verify bounds contain optimum
5. Increase `neps` for stricter convergence

### Problem: Step Sizes (VM) Not Adapting

**Solutions:**
1. Adjust `c` vector: larger values = faster adaptation
2. Change `step_mode`: try mode 1 (adaptive)
3. Check `ns` value: too large = slow adaptation
4. Monitor with `iprint=2` to see VM evolution

### Problem: Getting Invalid Function Evaluations

**Solutions:**
1. Return `istat=-1` in Fortran or raise exception in Python
2. Tighten bounds to valid region
3. Add penalty for constraint violations
4. Check for numerical overflow/underflow

### Problem: Results Not Reproducible

**Cause:** Different random seeds

**Solution:**
```python
# Set fixed seeds for reproducibility
optimizer.initialize(..., iseed1=1234, iseed2=5678)

# For testing multiple runs:
for run in range(10):
    optimizer.initialize(..., iseed1=1000+run, iseed2=2000+run)
```

### Problem: High Memory Usage

**Cause:** Large arrays or parallel workers

**Solutions:**
1. Reduce problem dimension if possible
2. Limit parallel workers: `--workers 4` instead of 16
3. Monitor with `iprint=0` to reduce output
4. Use smaller `nt` and `neps`

### Problem: "vm_min threshold" Error

**Cause:** Step sizes collapsed (no progress possible)

**Solutions:**
1. Increase `vm_min` tolerance
2. Widen bounds
3. Check if stuck in flat region
4. Use different initial point
5. Increase temperature

---

## Additional Resources

### API Documentation

Full API documentation: https://jacobwilliams.github.io/simulated-annealing/

### Key Papers

1. **Corana et al. (1987)**: [Original algorithm](https://dl.acm.org/doi/10.1145/29380.29864)
2. **Goffe et al. (1994)**: [Practical guidance](https://www.sciencedirect.com/science/article/abs/pii/0304407694900388)
3. **Kirkpatrick et al. (1983)**: [Foundational paper](https://science.sciencemag.org/content/220/4598/671)

### Example Files

- `example.py`: Basic Python usage
- `parallel_example.py`: OpenMP-style parallel (ctypes)
- `dask_parallel_example.py`: Distributed parallel with Dask
- `test/test.f90`: Fortran examples

### Getting Help

- GitHub Issues: https://github.com/jacobwilliams/simulated-annealing/issues
- Documentation: https://jacobwilliams.github.io/simulated-annealing/

---

## Appendix: Exit Codes

| Code | Meaning | Action |
|------|---------|--------|
| 0 | Success - convergence achieved | Results valid |
| 1 | Max evaluations exceeded | Increase `maxevl` or loosen `eps` |
| 3 | Negative initial temperature | Fix `t0` (must be ≥ 0) |
| 4 | User stop in function | Check `istat=-2` logic |
| 5 | Step sizes collapsed | See troubleshooting above |
| 99 | Internal error | Should not occur |

## Appendix: Quick Reference

```python
# Minimal working example
from sa_fortran import sa_fortran

optimizer = sa_fortran()
optimizer.initialize(n=5, lb=[-10]*5, ub=[10]*5, fcn=my_func, maximize=False)
result = optimizer.solve(x0=[0]*5, rt=0.85, t0=1.0)
print(result['x'], result['f'])
optimizer.destroy()
```

**Common parameter sets:**

```python
# Fast exploration
optimizer.initialize(..., rt=0.8, t0=5.0, ns=10, maxevl=10000)

# Thorough optimization (recommended)
optimizer.initialize(..., rt=0.85, t0=1.0, ns=20, nt=max(100,5*n), maxevl=100000)

# High-quality (slow)
optimizer.initialize(..., rt=0.95, t0=1.0, ns=20, nt=max(100,5*n), maxevl=500000)
```

---

**Version**: 1.0
**Last Updated**: 2026-06-22
**License**: See repository LICENSE file
