# Validation

Compiler used: GNU Fortran 14.2.0.

Strict flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. Kaplan-Meier / Nelson-Aalen counts and estimates.
2. Cox PH Efron fitting and likelihood improvement.
3. Start/stop counting-process KM and Cox code paths.
4. Concordance and log-rank statistics.
5. Parametric lognormal AFT location/scale recovery.
6. Aalen-Johansen probability conservation.
7. Fine-Gray low-level interval expansion and `surv_split`.
8. `pspline_basis` construction using the supplied splines translation.

All translated library sources and the supplied spline dependency were compiled
with the strict flags above. The final release is also rebuilt from a fresh
archive extraction before publication.
