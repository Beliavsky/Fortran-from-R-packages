# Validation

The translation is validated with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. The package's documented six-variable QP with all four methods.
2. Square and low-rank unconstrained box QPs with known analytic solutions.
3. Objective evaluation under both square and low-rank Hessian conventions.
4. Input validation and iteration-limit status.
5. KKT stationarity, equality feasibility, and dual nonnegativity.

The documented QP solution is approximately:

```text
alpha = 0.4761904762 1.0476190476 2.0952380952 2.9523809524 0 0
objective = -2.3809523810
```

The final release archive is also extracted into a clean directory and rebuilt
from only its contents before release.
