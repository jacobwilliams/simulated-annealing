# Cooling Schedule References

This document provides academic references and implementation details for the cooling schedules available in the simulated annealing library.

## Overview

The cooling schedule (temperature reduction schedule) is critical for simulated annealing performance. It controls how quickly the algorithm transitions from exploration (high temperature) to exploitation (low temperature). The library implements five well-established cooling schedules from the literature.

---

## 1. Geometric (Classical) Schedule

**Formula:** `T(k+1) = rt × T(k)` where `0 < rt < 1`

**Parameters:**
- `rt`: Temperature reduction factor (typically 0.8-0.95)

**Description:**
The classical geometric cooling schedule, where temperature decreases exponentially. This is the most commonly used schedule and the default in most SA implementations.

**References:**

1. **Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983)**
   "Optimization by Simulated Annealing"
   *Science*, 220(4598), 671-680.
   [DOI: 10.1126/science.220.4598.671](https://doi.org/10.1126/science.220.4598.671)
   *Original paper introducing simulated annealing with geometric cooling*

2. **Corana, A., Marchesi, M., Martini, C., & Ridella, S. (1987)**
   "Minimizing Multimodal Functions of Continuous Variables with the Simulated Annealing Algorithm"
   *ACM Transactions on Mathematical Software*, 13(3), 262-280.
   [DOI: 10.1145/29380.29864](https://doi.org/10.1145/29380.29864)
   *Recommends rt = 0.85 for most problems*

3. **Goffe, W. L., Ferrier, G. D., & Rogers, J. (1994)**
   "Global Optimization of Statistical Functions with Simulated Annealing"
   *Journal of Econometrics*, 60(1-2), 65-99.
   [DOI: 10.1016/0304-4076(94)90038-8](https://doi.org/10.1016/0304-4076(94)90038-8)
   *Practical guidance on choosing rt based on problem characteristics*

**Implementations:**
- Original Netlib SIMANN.F
- SciPy `scipy.optimize.dual_annealing` (default)
- MATLAB Global Optimization Toolbox `simulannealbnd`
- Apache Commons Math `SimulatedAnnealing`
- Most academic and commercial SA implementations

**Pros:**
- Simple and robust
- Well-understood behavior
- Extensive empirical validation
- Single parameter to tune

**Cons:**
- Can be slow to converge
- May require problem-specific tuning of rt

**Recommended Use:**
Default choice for most problems. Use rt ≈ 0.85-0.9 for general optimization.

---

## 2. Fast Annealing (Cauchy) Schedule

**Formula:** `T(k) = T₀ / (1 + k)`

**Parameters:**
- `T₀`: Initial temperature
- `k`: Iteration number

**Description:**
Fast annealing uses an inverse-linear cooling schedule, leading to faster temperature decay than geometric. Originally proposed for use with Cauchy distributions for step generation (hence the name), but can be used with any perturbation distribution.

**References:**

1. **Szu, H., & Hartley, R. (1987)**
   "Fast Simulated Annealing"
   *Physics Letters A*, 122(3-4), 157-162.
   [DOI: 10.1016/0375-9601(87)90796-1](https://doi.org/10.1016/0375-9601(87)90796-1)
   *Original fast annealing paper*

2. **Ingber, L. (1989)**
   "Very Fast Simulated Re-Annealing"
   *Mathematical and Computer Modelling*, 12(8), 967-973.
   [DOI: 10.1016/0895-7177(89)90202-1](https://doi.org/10.1016/0895-7177(89)90202-1)
   *Adaptive version with dynamic adjustments*

3. **Szu, H. H., & Hartley, R. L. (1987)**
   "Nonconvex Optimization by Fast Simulated Annealing"
   *Proceedings of the IEEE*, 75(11), 1538-1540.
   [DOI: 10.1109/PROC.1987.13916](https://doi.org/10.1109/PROC.1987.13916)
   *Application-focused overview*

**Implementations:**
- ASA (Adaptive Simulated Annealing) by Ingber
- Some neural network training algorithms
- Computer vision optimization

**Pros:**
- Faster convergence than geometric
- Good for high-dimensional problems
- No parameter tuning needed (beyond T₀)

**Cons:**
- May converge too quickly for multimodal functions
- Less exploration than geometric
- Can get stuck in local minima if T₀ too low

**Recommended Use:**
Problems with many local minima where faster exploration is beneficial. Good for initial exploration followed by local refinement.

---

## 3. Huang (Generalized) Schedule

**Formula:** `T(k) = T₀ / (1 + c × k)^d`

**Parameters:**
- `T₀`: Initial temperature
- `c`: Cooling parameter (cooling_param, default 1.0)
- `d`: Cooling exponent (cooling_exponent, default 1.0)

**Description:**
A generalized polynomial cooling schedule that includes fast annealing as a special case (d=1, c=1). The exponent d controls the cooling rate, allowing fine-tuning between fast and slow cooling.

**References:**

1. **Huang, M. D., Romeo, F., & Sangiovanni-Vincentelli, A. (1986)**
   "An Efficient General Cooling Schedule for Simulated Annealing"
   *Proceedings of IEEE International Conference on Computer-Aided Design*, 381-384.
   [DOI: 10.1109/ICCAD.1986.1620207](https://doi.org/10.1109/ICCAD.1986.1620207)
   *Original paper on polynomial cooling schedules*

2. **Hajek, B. (1988)**
   "Cooling Schedules for Optimal Annealing"
   *Mathematics of Operations Research*, 13(2), 311-329.
   [DOI: 10.1287/moor.13.2.311](https://doi.org/10.1287/moor.13.2.311)
   *Theoretical analysis of polynomial schedules*

3. **Nourani, Y., & Andresen, B. (1998)**
   "A Comparison of Simulated Annealing Cooling Strategies"
   *Journal of Physics A: Mathematical and General*, 31(41), 8373-8385.
   [DOI: 10.1088/0305-4470/31/41/011](https://doi.org/10.1088/0305-4470/31/41/011)
   *Comprehensive comparison of cooling schedules*

**Implementations:**
- VLSI design tools (original application)
- Some custom SA implementations
- Optimization frameworks with advanced cooling options

**Pros:**
- Flexible: can tune cooling rate with d
- Generalizes fast annealing (d=1)
- Can be faster than geometric with appropriate d

**Cons:**
- Two parameters to tune (c and d)
- Less commonly used, less empirical guidance
- Behavior less intuitive than geometric

**Recommended Use:**
Fine-tuning when geometric is too slow and fast annealing too aggressive. Try d=0.5 to 2.0 range.

---

## 4. Boltzmann (Logarithmic) Schedule

**Formula:** `T(k) = T₀ / log(1 + k + e)`

**Parameters:**
- `T₀`: Initial temperature
- `e = exp(1) ≈ 2.71828`: Base of natural logarithm

**Description:**
Theoretical schedule based on statistical mechanics. Provides very slow cooling with theoretical guarantees of convergence to global optimum (given sufficient time). Named after Ludwig Boltzmann's work on statistical mechanics.

**References:**

1. **Geman, S., & Geman, D. (1984)**
   "Stochastic Relaxation, Gibbs Distributions, and the Bayesian Restoration of Images"
   *IEEE Transactions on Pattern Analysis and Machine Intelligence*, PAMI-6(6), 721-741.
   [DOI: 10.1109/TPAMI.1984.4767596](https://doi.org/10.1109/TPAMI.1984.4767596)
   *Theoretical foundation for logarithmic cooling*

2. **Aarts, E. H. L., & Korst, J. (1989)**
   "Simulated Annealing and Boltzmann Machines"
   *Wiley*, Chichester, UK.
   ISBN: 978-0471921461
   *Comprehensive textbook including cooling schedule theory*

3. **Azizi, N., & Zolfaghari, S. (2004)**
   "Adaptive Temperature Control for Simulated Annealing: A Comparative Study"
   *Computers & Operations Research*, 31(14), 2439-2451.
   [DOI: 10.1016/S0305-0548(03)00197-7](https://doi.org/10.1016/S0305-0548(03)00197-7)
   *Empirical comparison including Boltzmann schedule*

**Implementations:**
- Some theoretical SA implementations
- Image processing (original application from Geman & Geman)
- Academic research implementations

**Pros:**
- Theoretical convergence guarantees
- Very thorough exploration
- Good for highly multimodal functions

**Cons:**
- Extremely slow convergence
- Impractical for most applications
- Requires very large iteration counts

**Recommended Use:**
Primarily of theoretical interest. Use when guaranteed convergence is more important than computation time, or for benchmarking other schedules.

---

## 5. Logarithmic (Modified) Schedule

**Formula:** `T(k) = T₀ / (1 + c × log(1 + k))`

**Parameters:**
- `T₀`: Initial temperature
- `c`: Cooling parameter (cooling_param, default 1.0)

**Description:**
A practical modification of the Boltzmann schedule that cools faster while retaining some theoretical properties. The parameter c allows tuning the cooling rate.

**References:**

1. **Geman, S., & Geman, D. (1984)**
   "Stochastic Relaxation, Gibbs Distributions, and the Bayesian Restoration of Images"
   *IEEE Transactions on Pattern Analysis and Machine Intelligence*, PAMI-6(6), 721-741.
   [DOI: 10.1109/TPAMI.1984.4767596](https://doi.org/10.1109/TPAMI.1984.4767596)
   *Foundation for logarithmic schedules*

2. **Triki, E., Collette, Y., & Siarry, P. (2005)**
   "A Theoretical Study on the Behavior of Simulated Annealing Leading to a New Cooling Schedule"
   *European Journal of Operational Research*, 166(1), 77-92.
   [DOI: 10.1016/j.ejor.2004.03.035](https://doi.org/10.1016/j.ejor.2004.03.035)
   *Analysis of modified logarithmic schedules*

3. **Abramson, D., Krishnamoorthy, M., & Dang, H. (1999)**
   "Simulated Annealing Cooling Schedules for the School Timetabling Problem"
   *Asia-Pacific Journal of Operational Research*, 16(1), 1-22.
   *Practical evaluation of logarithmic variations*

**Implementations:**
- Some operations research applications
- Combinatorial optimization tools
- Custom SA frameworks

**Pros:**
- Faster than pure Boltzmann
- Still provides good exploration
- Tunable with c parameter

**Cons:**
- Still slower than geometric
- Less commonly used
- Limited empirical guidance

**Recommended Use:**
Middle ground between Boltzmann and fast annealing. Use when geometric cooling gets stuck but need more exploration than fast annealing provides.

---

## Comparative Studies

Several papers have compared these cooling schedules empirically:

1. **Nourani, Y., & Andresen, B. (1998)**
   "A Comparison of Simulated Annealing Cooling Strategies"
   *Journal of Physics A: Mathematical and General*, 31(41), 8373-8385.
   [DOI: 10.1088/0305-4470/31/41/011](https://doi.org/10.1088/0305-4470/31/41/011)

2. **Salamon, P., Sibani, P., & Frost, R. (2002)**
   "Facts, Conjectures, and Improvements for Simulated Annealing"
   *SIAM Monographs on Mathematical Modeling and Computation*, Vol. 7.
   ISBN: 978-0898715088

3. **Boussaïd, I., Lepagnot, J., & Siarry, P. (2013)**
   "A Survey on Optimization Metaheuristics"
   *Information Sciences*, 237, 82-117.
   [DOI: 10.1016/j.ins.2013.02.041](https://doi.org/10.1016/j.ins.2013.02.041)

---

## Usage Recommendations

### Quick Selection Guide

| Problem Type | Recommended Schedule | Parameters |
|-------------|---------------------|------------|
| General optimization | Geometric (1) | rt = 0.85 |
| Highly multimodal | Geometric (1) | rt = 0.9-0.95 (slower) |
| High-dimensional | Fast annealing (2) | T₀ = 5-10 |
| Real-time/fast | Fast annealing (2) | T₀ = 1-5 |
| Thorough exploration | Logarithmic (5) | c = 1.0 |
| Fine-tuning | Huang (3) | d = 0.5-1.5 |
| Theoretical guarantee | Boltzmann (4) | (academic use) |

### Implementation in Other SA Codes

**Using geometric cooling (Schedule 1):**
- SciPy's `dual_annealing` and `basinhopping`
- MATLAB's `simulannealbnd`
- Apache Commons Math
- GNU Scientific Library (GSL) SA
- Original SIMANN.F from Netlib

**Using fast annealing (Schedule 2):**
- ASA (Adaptive Simulated Annealing) by Ingber
- Some GPU-accelerated SA implementations
- Fast SA variants in computer vision

**Using adaptive/custom schedules:**
- `simanneal` Python package (allows custom schedules)
- PyGMO/PaGMO optimization framework
- DEAP (Distributed Evolutionary Algorithms in Python)

---

## Implementation Details

### Theoretical Convergence

For asymptotic convergence to global optimum:
- **Logarithmic schedules** (Boltzmann, Modified): Guaranteed convergence
- **Polynomial schedules** (Huang with d ≥ 1): Guaranteed if d large enough
- **Fast annealing**: No guarantee (but empirically good)
- **Geometric**: No guarantee (but most practical)

### Practical Performance

Empirical studies consistently show:
1. **Geometric** (rt ≈ 0.85-0.9) gives best balance for most problems
2. **Fast annealing** is 2-5× faster but may miss global optimum
3. **Logarithmic** schedules too slow for practical use
4. **Huang** with d ≈ 1 similar to fast annealing

### Combining with Step Size Adaptation

All schedules work with the adaptive step size (VM) adjustment. The combination of:
- Temperature schedule (exploration control)
- Step size adaptation (search range control)

provides robust optimization across problem types.

---

## References by Schedule

### Geometric (Schedule 1)
- Kirkpatrick et al. (1983) - Original SA paper
- Corana et al. (1987) - Continuous variable SA
- Goffe et al. (1994) - Practical guidelines

### Fast Annealing (Schedule 2)
- Szu & Hartley (1987) - Fast SA introduction
- Ingber (1989) - Very fast SA
- Szu & Hartley (1987) - IEEE overview

### Huang (Schedule 3)
- Huang et al. (1986) - Polynomial schedules
- Hajek (1988) - Theoretical analysis
- Nourani & Andresen (1998) - Comparative study

### Boltzmann (Schedule 4)
- Geman & Geman (1984) - Theoretical foundation
- Aarts & Korst (1989) - Textbook treatment
- Azizi & Zolfaghari (2004) - Empirical comparison

### Logarithmic (Schedule 5)
- Geman & Geman (1984) - Foundation
- Triki et al. (2005) - Modified logarithmic study
- Abramson et al. (1999) - Practical evaluation

---

## Additional Resources

### Books
- **Aarts, E., & Korst, J. (1989)**: *Simulated Annealing and Boltzmann Machines*. Comprehensive coverage of cooling theory.
- **Van Laarhoven, P. J., & Aarts, E. H. (1987)**: *Simulated Annealing: Theory and Applications*. Mathematical foundations.

### Review Papers
- **Delahaye, D., Chaimatanan, S., & Mongeau, M. (2019)**: "Simulated Annealing: From Basics to Applications", *Handbook of Metaheuristics*, 1-35.
- **Henderson, D., Jacobson, S. H., & Johnson, A. W. (2003)**: "The Theory and Practice of Simulated Annealing", *Handbook of Metaheuristics*, 287-319.

### Online Resources
- Original SIMANN.F: https://www.netlib.org/opt/simann.f
- SciPy documentation: https://docs.scipy.org/doc/scipy/reference/optimize.html
- ASA code: http://www.ingber.com/#ASA

---

**Document Version**: 1.0
**Last Updated**: 2026-06-22
**Maintainer**: See repository contributors
