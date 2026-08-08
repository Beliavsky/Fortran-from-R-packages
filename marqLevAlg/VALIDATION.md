# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. upstream two-parameter quadratic with numerical derivatives, analytical
   gradient, and supplied information matrix;
2. upstream four-parameter Powell-style objective;
3. maximization sign handling;
4. `deriva` and `deriva_grad` formulas;
5. random-intercept LMM analytical gradient checked against central finite
   differences.

Both examples are compiled and executed by the strict validation script.
Translated source files are also checked for the standard 132-column free-form
line limit.
