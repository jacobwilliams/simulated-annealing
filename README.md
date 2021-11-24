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

  * Corana et al., "Minimizing multimodal functions of continuous variables
    with the "simulated annealing" algorithm", september 1987
    (vol. 13, no. 3, pp. 262-280),
    acm transactions on mathematical software.
  * Goffe, Ferrier and Rogers, "Global optimization of statistical functions
    with simulated annealing", journal of econometrics, vol. 60, no. 1/2,
    jan./feb. 1994, pp. 65-100.
