![simulated-annealing](media/logo.png)
============

A modern Fortran simulated annealing optimization method. ***A work in progress.***
# Status

![Build Status](https://github.com/jacobwilliams/simulated-annealing/actions/workflows/CI.yml/badge.svg)

### Building

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library and tests cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```FoBis.py rule --execute makedoc -f simulated-annealing.fobis```

### See also

  *  https://www.netlib.org/opt/simann.f

### Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/simulated-annealing/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### References

  * Corana et al., "[Minimizing multimodal functions of continuous variables
    with the "simulated annealing" algorithm](https://dl.acm.org/doi/10.1145/29380.29864)", september 1987
    (vol. 13, no. 3, pp. 262-280),
    acm transactions on mathematical software.
  * Goffe, Ferrier and Rogers, "[Global optimization of statistical functions
    with simulated annealing](https://www.sciencedirect.com/science/article/abs/pii/0304407694900388)", journal of econometrics, vol. 60, no. 1/2,
    jan./feb. 1994, pp. 65-100.
  * S. Kirkpatrick, C. D. Gelatt Jr., M. P. Vecchi, "[Optimization by Simulated Annealing](https://pdfs.semanticscholar.org/e893/4a942f06ee91940ab57732953ec6a24b3f00.pdf)", Science 13 May 1983, Vol. 220, Issue 4598, pp. 671-680
  * W. L. Goffe, [SIMANN: A Global Optimization Algorithm using Simulated Annealing](https://www.researchgate.net/publication/24015773_SIMANN_A_Global_Optimization_Algorithm_using_Simulated_Annealing), Studies in Nonlinear Dynamics & Econometrics, De Gruyter, vol. 1(3), pages 1-9, October 1996.