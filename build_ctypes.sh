
#!/bin/bash

# build the python/ctypes interface

set -e

# build the fortran shared lib:
fpm install --prefix ./sa_fortran/lib --profile release

# build the python package:
hatch build