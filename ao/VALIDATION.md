# Validation

The translated library is compiled with:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. sequential Rosenbrock AO;
2. random partitions with an analytical gradient;
3. custom overlapping blocks and bounds;
4. maximization;
5. analytical-gradient/Hessian Newton blocks;
6. random-partition and estimate-splitting helpers;
7. multiple-process Cartesian-product optimization.

The final release is also extracted into a fresh directory and rebuilt from
only the archive contents. See `scripts/test_gfortran.sh` and
`scripts/test_gfortran.bat` for the compiler commands.
