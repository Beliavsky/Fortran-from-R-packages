# Validation

The translated library was compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

The test suite covers:

1. SPSO2007 on a bounded quadratic.
2. SPSO2011 on Rastrigin.
3. Both synchronous/vectorized SPSO2007 and SPSO2011 paths.
4. Hybrid refinement with an analytical gradient and with numerical gradients.
5. Maximization through negative `fnscale`.
6. Restart-limit termination and convergence code.
7. Scalar bounds and supplied initial point.
8. Upstream benchmark objective values and benchmark metadata.

The example programs are also compiled and executed by the direct compiler
validation scripts.

FPM is not required by those scripts; on an FPM installation use:

```text
fpm build
fpm test
```
