# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

## Permanent regression programs

1. `test_continuous` — dense LP/equality, SOC, and exponential cone.
2. `test_mixed_integer` — boolean/general-integer branch-and-bound examples.
3. `test_lifecycle_csc` — lifecycle and CSC compatibility/status behavior.
4. `test_status_validation` — primal infeasible, unbounded, and invalid dimensions.
5. `test_mixed_cones` — simultaneous linear and SOC blocks.
6. `test_sparse_backend` — native sparse SOC, sparse exponential cone, and a 2,000-variable diagonal sparse LP.
7. `test_sparse_ldl` — direct sparse quasi-definite LDL solve and residual check.
8. `test_sparse_update` — sparse workspace updates remain on the sparse backend.
9. `test_v04_features` — v0.4 architecture regression covering:
   - persistent symbolic cache reuse;
   - persistent warm-start reuse;
   - same-pattern sparse matrix-value symbolic reuse;
   - warm invalidation after matrix-value changes;
   - cone-preserving equilibration over a `1e16` coefficient range;
   - AMD-style fill reduction versus identity ordering;
   - sparse primal-unbounded certificate;
   - sparse primal-infeasible dual certificate;
   - dual exponential-cone certificate mapping;
   - sparse ECOS_BB and symbolic node reuse;
   - inaccurate-optimal status offset.

All **9/9** tests pass.

## Examples

All four examples compile and run under the same strict flags:

- `soc_example`
- `mixed_integer_example`
- `sparse_large_lp`
- `workspace_reuse`

The workspace example reports:

```text
first symbolic analyses: 1
second symbolic analyses: 0
second cached symbolic reuses: 1
second cached warm starts: 1
```

## 10,000-variable sparse benchmark

The strict-build `sparse_large_lp` executable reports:

```text
variables: 10000
max |x-1|:             6.250001e-10
KKT upper nnz:         30000
LDL off-diagonal nnz:  10000
symbolic analyses:         1
numeric factorizations:    4
LDL/KKT fill ratio:     0.333
ordering CPU seconds:   about 0.24
factor CPU seconds:     about 0.023
refinement CPU seconds: about 0.0085
```

Using `/usr/bin/time` on the same validation container and the same strict `-O0 -fcheck=all` build:

```text
v0.4.0: elapsed ~= 0.30 s, peak RSS ~= 11.1 MB
v0.2.0: elapsed ~= 0.67 s, peak RSS ~= 9.4 MB
```

The extra memory is primarily the dynamic adjacency/heap state used by the AMD-style ordering. These timings are environment-specific and are included as a regression sanity check, not as a general performance guarantee.

## Ordering regression

For a 25-by-25 star-pattern symmetric matrix with the hub numbered first:

```text
identity symbolic L off-diagonal nnz: 300
AMD-style symbolic L off-diagonal nnz: 24
RCM symbolic L off-diagonal nnz:       24
```

This verifies that the minimum-degree implementation avoids the intentionally poor natural ordering.

## Equilibration regression

A two-variable LP uses diagonal constraint coefficients `1e-8` and `1e8`, a ratio of `1e16`. With default equilibration the solution satisfies `x ~= (1,1)` and the reported row/column scale ranges are nontrivial by more than five orders of magnitude.

## Sparse certificate regression

Unbounded LP:

```text
minimize -x
subject to x >= 0
```

returns `ECOS_DINF`, `dual_certificate_valid=.true.`, and a normalized ray approximately `d=1`.

Infeasible LP:

```text
x >= 1
x <= 0
```

returns `ECOS_PINF`, `primal_certificate_valid=.true.`, with a normalized dual ray.

An infeasible exponential-cone model separately validates the dual exponential mapping used in the certificate problem.

## MatrixExtra adapter

The optional MatrixExtra adapter is built after the root core. Its example constructs a sparse LP through MatrixExtra/Matrix sparse types, converts sparse arrays only, and verifies that the ECOS problem has no allocated dense `G` matrix.

## Release hygiene

Before packaging:

- translated Fortran files are scanned for lines over the standard 132-column free-form limit;
- root and adapter FPM TOML manifests are parsed;
- generated object/module/build directories are removed;
- the ZIP is extracted into a fresh directory;
- the complete strict core suite and examples are rebuilt/rerun from only the archive contents;
- the optional MatrixExtra adapter is rebuilt against the freshly extracted root package.

FPM itself is not installed in the validation container, so manifests are parsed independently and the identical FPM source/test/example trees are compiled directly with `gfortran`.
