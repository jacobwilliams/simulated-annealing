A modern Fortran simulated annealing optimization method. A work in progress.

### Building

A [FoBiS](https://github.com/szaghi/FoBiS) configuration file (`simulated-annealing.fobis`) is also provided that can also build the library and examples. Use the `mode` flag to indicate what to build. For example:

  * To build all the examples using gfortran: `FoBiS.py build -f simulated-annealing.fobis -mode tests-gnu`
  * To build all the examples using ifort: `FoBiS.py build -f simulated-annealing.fobis -mode tests-intel`
  * To build a static library using gfortran: `FoBiS.py build -f simulated-annealing.fobis -mode static-gnu`
  * To build a static library using ifort: `FoBiS.py build -f simulated-annealing.fobis -mode static-intel`

  The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`

  To generate the documentation using [ford](https://github.com/cmacmackin/ford), run: ```FoBis.py rule --execute makedoc -f simulated-annealing.fobis```

### See also

  *  https://www.netlib.org/opt/simann.f

### References

  * Corana et al., "Minimizing multimodal functions of continuous variables
    with the "simulated annealing" algorithm", september 1987
    (vol. 13, no. 3, pp. 262-280),
    acm transactions on mathematical software.
  * Goffe, Ferrier and Rogers, "Global optimization of statistical functions
    with simulated annealing", journal of econometrics, vol. 60, no. 1/2,
    jan./feb. 1994, pp. 65-100.
