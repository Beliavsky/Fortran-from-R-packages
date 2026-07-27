# Validation

The translated source was validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0
```

Results:

```text
test_densities: PASS
test_fitting: PASS
test_objectives: PASS
test_pricing: PASS
```

The application and both examples also compiled and ran successfully.

The tests cover:

- BSM put-call parity and nonnegative prices
- equivalence of a one-component lognormal mixture and BSM
- Edgeworth reduction to its lognormal baseline
- constant-volatility Shimko equivalence to BSM
- generalized-beta CDF and density behavior
- mixture-density identities
- recovery of synthetic BSM distribution parameters
- recovery of synthetic implied volatility
- quadratic Shimko smile fitting
- extraction of risk-free and dividend rates from put-call parity
- point-density finite differences
- objective-function values at synthetic generating parameters

The final release archive is extracted and rebuilt during packaging validation.
FPM was not installed in the validation environment, so direct GNU Fortran
compilation was used. `fpm.toml` was parsed as TOML and follows FPM automatic
source, application, example, and test discovery.

GNU ld reports that `rnd_fitting.o` requests an executable stack. This is caused
by GNU Fortran trampolines for internal procedure callbacks passed to the
Nelder-Mead optimizer. It does not affect the tested results. Compilers that
implement procedure closures without trampolines may not emit this warning.
