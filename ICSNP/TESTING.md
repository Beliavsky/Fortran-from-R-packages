# Testing

The project ships four test programs:

- `test_pairs_location`: pair operations and location estimators;
- `test_shapes`: robust shape/scatter and spatial-sign estimators;
- `test_location_tests`: Hotelling and marginal location tests;
- `test_independence_hp`: marginal/IC independence and Tyler-angle tests.

Validation configuration used during translation:

```text
gfortran 14.2.0
-std=f2018 -Wall -Wextra -Werror -pedantic
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -g
```

An optimized configuration using `-O3 -Werror` is also tested.

The environment used for this translation did not contain the `fpm`
executable, so the included GNU Fortran scripts compile the same source, test,
example, and application targets directly. The `fpm.toml` manifest follows FPM
conventions and uses a plain numeric version.
