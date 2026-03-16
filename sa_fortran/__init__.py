"""Python interface to the Fortran simulated annealing library using ctypes.

Note: the shared lib must be built using:
`fpm install --prefix ./sa_fortran/lib --profile release`
"""

from .interface import (sa_fortran,
                        CALLBACK_FUNC,
                        CALLBACK_N_INPUTS,
                        CALLBACK_PARALLEL_INPUT,
                        CALLBACK_PARALLEL_OUTPUT,
                        CALLBACK_REPORT)

__version__ = '0.0.1'

__all__ = ['sa_fortran',
           'CALLBACK_FUNC',
           'CALLBACK_N_INPUTS',
           'CALLBACK_PARALLEL_INPUT',
           'CALLBACK_PARALLEL_OUTPUT',
           'CALLBACK_REPORT']
