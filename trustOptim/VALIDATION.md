# Validation

The Fortran translation is compiled with GNU Fortran using

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

All **6/6 test programs pass** under the strict flags. Regression coverage:

1. Paired Rosenbrock with BFGS, SR1, and sparse analytic Hessian.
2. BFGS Cholesky preconditioner and sparse modified-Cholesky preconditioner.
3. Maximization through `function_scale_factor = -1`.
4. Sparse lower-triangle storage, dense round trip, and symmetric matvec.
5. Binary-choice gradient finite-difference check and sparse maximization.
6. Explicit negative-curvature Steihaug-CG branch and max-iteration status.

The final release is also unpacked into a fresh directory and rebuilt/tested
from only the archive contents.
