# Validation

Validated with GNU Fortran 14.2.0.

Library sources, tests, and examples were compiled using strict settings:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -Wimplicit-interface \
  -fcheck=all ...
```

No nonstandard free-form line-length option is required.

## Tests

1. `test_architecture_v020.f90`
   - Stiefel ParamSet=2 constructed retraction
   - Sphere ParamSet=4 locking/Beta scaling
   - LowRank horizontal projection and metric-preserving scaled transport
   - Armijo, weak Wolfe, strong Wolfe, exact, and custom line-search paths
   - non-default Broyden-family Phi
2. `test_euclidean_methods.f90`
   - all eleven solver names on an analytical quadratic
3. `test_sphere.f90`
   - LRBFGS on a sphere linear objective
4. `test_stiefel.f90`
   - Stiefel optimization and orthogonality
5. `test_remaining_manifolds.f90`
   - Grassmann, OrthGroup, and SPD optimization
6. `test_spd_and_lowrank.f90`
   - SPD and LowRank tangent/retraction feasibility
7. `test_product_numeric.f90`
   - product manifold plus objective-only numerical differentiation

All seven tests pass with the strict flags above.

Both examples compile and run. `sphere_example` reaches the normalized exact
direction for its linear objective, and `brockett_stiefel` reaches objective 4
for its diagonal Brockett problem.

FPM itself is not installed in this validation container. `fpm.toml` is parsed
with Python's TOML parser, and the standard FPM `src`/`test`/`example` layout is
used. The final archive is also unpacked into a fresh directory and rebuilt
from only its contents before release.
