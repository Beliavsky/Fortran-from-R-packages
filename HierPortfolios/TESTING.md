# Testing

## FPM

```sh
fpm test
fpm run
fpm run --example example_hrp_herc
fpm run --example example_hcaa_gap
fpm run --example example_dhrp_bounds
```

## GNU Fortran scripts

Runtime-checked build:

```sh
./scripts/test_gfortran.sh
```

Optimized warning-as-error build:

```sh
./scripts/test_gfortran_optimized.sh
```

The strict configuration enables Fortran 2018 conformance, warnings as errors,
bounds and runtime checks, and traps for invalid, divide-by-zero, and overflow
floating-point exceptions.

The four tests cover HRP, HCAA with fixed and automatic cluster counts, HERC,
DHRP bounds and `tau`, deterministic gap simulation, link validation, and
invalid covariance/bound inputs.
