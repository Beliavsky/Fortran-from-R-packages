# API

## Types

- `cone_result`: projected/fitted vector, coefficients, active face, active matrix, degrees of freedom, steps, status.
- `qprog_result`: QP solution, objective, active face, active matrix, degrees of freedom, steps, status.
- `constreg_result`: constrained/unconstrained fits, coefficients, covariance, standard errors, t statistics, p-values, confidence intervals, optional cone-test p-value.
- `shapereg_result`: parametric coefficients, constrained/linear fits, standard errors, p-values, SSE values, optional cone-test p-value.
- `qr_result`: orthonormal column basis and numerical rank.

## Core procedures

`cone_a(y, amat, result [, weights, start_face, max_iter])`

`cone_b(y, delta, result [, vmat, weights, start_face, max_iter])`

`qprog(q, c, amat, b, result [, start_face, max_iter])`

`make_delta(x, shape, delta, status)`

`constreg_fit(y, xmat, amat, result [, weights, test, nloop, nsim_cov])`

`shapereg_fit(y, t, shape, result [, xmat, weights, test, nloop])`

`qr_decomp(x, result)`

`check_irreducible(edges, keep, reducible, equal_edges, status)`

`seed_rng(seed)` makes Monte-Carlo inference reproducible.
