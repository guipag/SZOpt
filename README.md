# SZOpt
Fortran source code for calculating sound zone filters in fast time.

# Compilation

The algorithm is implemented in Fortran 90 in the file `solve_f_blas.f90`, using three BLAS routines: 'dgemv', 'dsymv' and 'ddot'. To optimize computation time, it is possible to use a BLAS library using parallelism (e.g. openBLAS). In the support paper the file was compiled using f2py with numpy 2.0.2, scipy 1.13.1, meson 1.5.2 and openBLAS with the comand line:

```
f2py -c solve_f_blas.f90 -m solve_f_openblas -lopenblas -llapack --opt="-O3 -march=native -ffast-math" --backend=meson
```
