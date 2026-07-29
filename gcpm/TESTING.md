# Validation

The project was compiled on 2026-07-27 with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
```

The automated test program checks:

- standard normal CDF/inverse-CDF accuracy;
- portfolio validation and discretization;
- analytical CreditRisk+ PDF and CDF properties;
- analytical expected-loss agreement;
- VaR monotonicity and expected-shortfall consistency;
- standard-deviation contribution add-up;
- finite analytical EC/VaR/ES contributions;
- original-format CSV parsing, including decimal commas;
- CreditMetrics Monte Carlo execution;
- empirical CDF normalization; and
- finite simulation EC/VaR/ES contributions.

All checks passed. The analytical and CRP Monte Carlo demonstration also ran
successfully and produced closely matching expected loss and tail-risk results.

The container used for translation did not include the `fpm` executable, so the
manifest was syntax-checked as TOML and the same source/test targets were built
directly with GNU Fortran. The project follows the standard FPM `src/`,
`example/`, and `test/` layout.
