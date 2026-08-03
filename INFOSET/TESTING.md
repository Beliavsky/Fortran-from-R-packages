# Testing

Four tests are included:

- `test_core`: gross returns and rolling windows.
- `test_mixture`: two-lognormal fitting and information-set extraction.
- `test_left_risk`: multi-asset rolling Left Risk.
- `test_portfolio`: all four portfolio strategies and summary statistics.

Validation configurations:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

and

```text
-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror
```
