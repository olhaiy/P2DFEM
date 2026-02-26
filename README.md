# Batree: An MFEM-based SPM, SPMe and P2D solver

_Nuno Nobre, Karthikeyan Chockalingam, Daniel Ward and Olha Yaman, STFC Hartree Centre_

[![Integration Tests](https://github.com/karthichockalingam/P2DFEM/actions/workflows/tests.yml/badge.svg)](https://github.com/karthichockalingam/P2DFEM/actions/workflows/tests.yml)

###

Compile with `make batree`.

Sample runs:
```
mpirun -np 1 ./batree -m SPM
mpirun -np 3 ./batree -m SPMe
mpirun -np 4 ./batree -m P2D
```

Under active development. No explicit time integration methods are supported
at this time. Use `-m` or `--method` to select from the three electrochemical
models (SPM by default). At the moment, the program will only perform a single
CC discharge cycle until the time specified with `-tf` or `--t-final` (3600s by
default). Run with `-h` or `--help` for all available options.