# Validation

The release was compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

All seven permanent tests and the example pass:

- `test_gpava`
- `test_repeated`
- `test_active`
- `test_partial_order`
- `test_solvers`
- `test_mregnn`
- `test_chain_stress`

`test_chain_stress` checks 100 randomized weighted chain problems and requires
the independently implemented active-set LS fit to agree with PAVA.

Additional validation against independent SciPy constrained optimizers was
performed outside the permanent Fortran tests:

- 120 randomized weighted LS problems with partial-order constraints:
  0 failures, worst fitted-coordinate discrepancy about `4.00e-7`, and worst
  objective discrepancy about `2.87e-12`;
- 80 randomized `mregnn` problems:
  0 failures, worst fitted-coordinate discrepancy about `2.31e-7`, and worst
  objective discrepancy about `7.64e-13`.

The coordinate differences are dominated by the stopping tolerances of the
independent SLSQP reference; objective values agree much more tightly.
