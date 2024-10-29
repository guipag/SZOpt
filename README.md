# SZOpt
Fortran source code for calculating sound zone filters in fast time.

# Compilation

The algorithm is implemented in Fortran 90 in the file `solve_f_blas.f90`, using three BLAS routines: 'dgemv', 'dsymv' and 'ddot'. To optimize computation time, it is possible to use a BLAS library using parallelism (e.g. openBLAS). In the support paper the file was compiled using f2py with numpy 2.0.2, scipy 1.13.1, meson 1.5.2 and openBLAS with the comand line:

```
f2py -c solve_f_blas.f90 -m solve_f_openblas -lopenblas -llapack --opt="-O3 -march=native -ffast-math" --backend=meson
```

# Usage

The algorithm solve the same equation than (Simon Galvez et al., 2015) [[1]](#1). The `HB` matrix is the convolution matrix containing the impulse responses of all the microphones in the bright zone, and similary for the `HD` matrix with the dark zones RIs. The `d` vector is the vector containing the desired pressure in the bright zone. The $\lambda$ and $\beta$ parameters are the same than in (Simon Galvez et al., 2015) [[1]](#1).

To use simple gradient descent the function to use is `solvegd`:
```
solvegd(HB, HD, d, beta, lambda, [nbIter])
```
To use conjugate gradient descent the function to use is `solvegdc`:
```
solvegdc(HB, HD, d, beta, lambda, [nbIter])
```
The parameter nbIter is optional, by default it is equal to the number of columns in the HB (or HD) matrix.

# References

<a id="1">[1]</a> Simon Galvez, M., Elliott, S., and Cheer, J. (2015). “Time domain optimization of filters used in a loudspeaker array for personal audio,” IEEE/ACM Transactions on Audio, Speech, and Language Processing 23(11), 1869–1878.
