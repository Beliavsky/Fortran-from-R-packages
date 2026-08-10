# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Five regression programs pass:

1. Lawson-Hanson published LSI/LDP/LSEI examples.
2. Direct numerical reference values from the original fixed-form
   NNLS/PNNLS/HFTI routines.
3. QP and partial-nonnegative QP transformations.
4. Box bounds and infeasibility status.
5. `indx` and `matMaxs` utilities.

The native-reference case agrees with the original source to approximately
machine precision for coefficients, HFTI rank, and residual norms.
