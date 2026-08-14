# Validation

Validation was performed with GNU Fortran 14.2.0.

## Strict build

All source, tests, and the example compile with:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace
```

The checked test executables report:

```text
test_core: PASS
test_population: PASS
test_optimizer: PASS
```

`test_optimizer` exercises the full optimizer on the Schaffer two-objective problem, verifies parameter bounds, checks that the final archive is Pareto-only, and enables sensitivity derivatives.

`test_population` performs 200 randomized archive/population reductions and checks archive/population caps, index validity, uniqueness, and disjointness.

## Delaunay differential test

A development-only driver was compared against `scipy.spatial.Delaunay` (SciPy 1.17.0) on 150 independent nondegenerate random point clouds:

- 80 cases in 2 dimensions;
- 50 cases in 3 dimensions;
- 20 cases in 4 dimensions.

After canonicalizing vertex order, the complete simplex sets matched exactly in all 150 cases (0 mismatches).

This test found an earlier super-simplex sizing weakness: 16/150 cases initially omitted one hull simplex. Enlarging the super-simplex resolved all of those discrepancies; the corrected implementation is the one included in this archive.

## Pareto differential test

The Fortran `pareto` and `dominate` routines were compared with an independent brute-force Python reference on 300 randomized matrices with 2-5 objectives and 3-24 rows. Values were rounded deliberately to create ties, and many cases included exact duplicate rows.

Result: 0 mismatches for both Pareto membership and complete onion-peel ranks.

## End-to-end example

The included Schaffer example completed under runtime checking and produced a non-empty Pareto front while respecting the specified bounds. Since the optimizer is stochastic, front size and exact objective values depend on the RNG seed and runtime implementation.

## FPM note

The sandbox used for this translation did not have an `fpm` executable installed. `fpm.toml` was parsed successfully as TOML, and the exact FPM source/test/example layout was compiled and run directly with gfortran using the strict flags above.
