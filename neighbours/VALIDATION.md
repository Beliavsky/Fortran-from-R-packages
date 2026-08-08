# Validation

The translated core is validated with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Test programs cover:

1. numeric zero-sum, budget-range, no-budget, active and `A*x` update behavior;
2. logical unconstrained, fixed-cardinality, bounded-cardinality and active behavior;
3. real/integer/character permutation moves and `next_subset` enumeration;
4. random logical/numeric vector generation and vector comparison;
5. the 5/10/40 portfolio neighbourhood.

Both FPM examples are compiled/run by the strict validation script.

The optional NMOF adapter is additionally compiled against the supplied NMOF
translation and used with `nmof_optimization::local_search`; the demonstration
converges the test simplex problem to approximately equal 0.2 weights.

FPM itself was not installed in the validation container. `fpm.toml` files are
syntax-checked with a TOML parser, and the source/test/example layout is rebuilt
directly with gfortran using the strict flags above.
