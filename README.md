# SZOpt
Fortran source code for calculating sound zone filters in fast time.

# Compilation
```
f2py -c solve_f_blas.f90 -m solve_f_openblas -lopenblas -llapack --opt="-O3 -march=native -ffast-math" --backend=meson
```
