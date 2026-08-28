# Validation

The release was compiled with gfortran 14.2 using:

```
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

The permanent test covers:

- upstream PCG32 raw words
- ordinary and Owen-scrambled Sobol, including maximum supported dimensions
- Faure Halton, including dimension 256
- deterministic randomized-Halton reproducibility and range checks
- PJ, PMJ, PMJ blue-noise, PMJ(0,2), and PMJ(0,2) blue-noise

A separately compiled C++ reference program using the exact headers from the
supplied `spacefillr` source was used for extended differential validation.
The comparison exercised:

- sequence indices 0, 1, 2, 3, 7, 31, 255, 65535, 123456789, and 2^32-1
- ordinary Sobol dimensions 0, 1, 2, 31, 511, and 1023
- Owen-Sobol dimensions 0, 1, 2, 31, 1023, 4095, and 21200
- Faure Halton dimensions 0, 1, 2, 3, 10, 100, and 255
- 256 points from each of PJ, PMJ, PMJ blue-noise, PMJ(0,2), and PMJ(0,2)
  blue-noise

All 2760 compared floating-point values matched the C++ output exactly when
parsed as double precision values.

Randomized Halton is not included in the bitwise C++ comparison because the
upstream algorithm delegates permutation generation to implementation-dependent
`std::shuffle`. Its Fortran implementation is instead tested for repeatability,
seed sensitivity, and [0,1) range validity.

The example and permanent test both pass in the optimized strict build.
