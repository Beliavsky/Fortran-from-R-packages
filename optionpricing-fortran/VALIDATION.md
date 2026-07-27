# Validation report

## Environment

- GNU Fortran 14.2.0
- Python 3.13.5 for independent deterministic reference calculations
- Fortran 2018 free-form source
- No external numerical libraries

The `fpm` executable was not installed in the validation environment. The
manifest was parsed with Python's TOML parser, and the same automatic `src`,
`test`, `app`, and `example` layout used by FPM was compiled directly with GNU
Fortran.

## Build modes

The exact release tree passed two independent builds:

1. Debug: `-O0 -g -fcheck=all -fbacktrace`
2. Optimized: `-O2`

Both used:

```text
-std=f2018
-Wall
-Wextra
-Wconversion-extra
-Wimplicit-interface
-Werror
```

## Test results

```text
test_analytic: PASS
test_european: PASS
test_monte_carlo: PASS
test_qmc: PASS
debug: PASS
optimized: PASS
manifest_audit: PASS
source_audit: PASS
release_tree_audit: PASS
checksum_audit: PASS
validation: PASS
```

The demo and both examples also compiled and ran in both build modes.

## Deterministic references

Independent Python calculations give:

```text
Black-Scholes call price: 4.614997129602855
Black-Scholes call delta: 0.5694601832076737
Black-Scholes call gamma: 0.03928800094473793

Curran ECV price: 6.154908593215828
Curran ECV delta: 0.5938761220547435
Curran ECV gamma: 0.03056755530266579
```

The Fortran values agree within approximately 3e-13. The analytic test also
checks:

- Curran lower-bound control expectation
- Quadratic control expectation and its first two spot derivatives
- Lord approximation at the package's standard parameter set
- Finite-difference delta and gamma of the Curran expectation
- Symmetry and nonnegative diagonal of the conditional covariance matrix

## Monte Carlo validation

The stochastic tests check:

- Repeatability with a fixed seed on the same compiler
- Naive price, pathwise delta, and mixed LR/PW gamma estimators
- New-control-variate likelihood-ratio estimates
- Conditional Monte Carlo with six control quantities
- Pilot-run multivariate regression
- Strong error reduction from naive MC to NCV/LR and full CMC/QCV
- Fixed caller-supplied-normal reference values independent of the RNG

At the standard Asian option parameters, the full MC estimate is close to:

```text
price = 6.1560406
delta = 0.5938539
gamma = 0.0305780
```

## QMC validation

The QMC test checks:

- Exact initial Korobov lattice points
- Random shifting and Baker transformation execution
- PCA factor covariance identity
- `pca`, `pcamain`, `lt`, and `ltpca` matrix construction
- Naive randomized-QMC estimates
- Full leave-one-out CMC/QCV estimates
- Error reduction relative to naive QMC
- The high-level `asian_call` dispatcher

## Source and provenance audit

The audit verifies:

- Valid TOML and expected SPDX license expression
- ASCII-only translated release text
- SPDX headers and `implicit none` in every Fortran file
- Maximum free-form line length of 132 characters
- Absence of object files, module files, archives, and executables
- SHA-256 manifests for the supplied ZIP, original package files, and
  translated files

## Linker note

GNU ld reports that the object containing internal callback procedures requires
an executable stack. This is caused by GNU Fortran's trampoline implementation
for nested procedure callbacks used by quadrature and root finding. It does not
affect the numerical tests. Hardened environments may eliminate the warning by
refactoring those callbacks into module-level procedures with explicit context
objects.

## Windows compiler portability regression

A Windows FPM run produced
`1.3997081976970311e-4` instead of the GNU/Linux reference
`1.3997081977210567e-4` for one fixed conditional-Monte-Carlo value. The
absolute difference is about `2.4e-15`; it arises from compiler and math-library
rounding in `erfc`, exponentials, cancellation-prone differences, and Newton
iteration. The original absolute tolerance of `2e-16` was below a practical
cross-platform rounding bound.

`test_monte_carlo` now uses a combined absolute/relative comparison with an
absolute tolerance of `1e-14` and a relative tolerance of `1e-10`. This still
requires approximately 10 to 11 significant decimal digits and does not alter
any library routine or reference value.
