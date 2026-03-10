#!/bin/bash
set -e

# Build Fortran library with fpm
fpm install --prefix ./install --profile release

# Build f2py extension
python -m numpy.f2py -c \
    --include-paths $(pwd)/install/include \
    $(pwd)/src/simulated_annealing.pyf \
    $(pwd)/install/lib/libsimulated-annealing.a

mv $(pwd)/simulated_annealing_f2py.*.so $(pwd)/sa_fortran