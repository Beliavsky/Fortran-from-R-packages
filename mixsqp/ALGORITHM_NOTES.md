# Algorithm notes

## Objective

For nonnegative likelihood matrix `L` and normalized nonnegative row weights
`w`, mix-SQP minimizes

    f(x) = -sum_i w(i) log((L x)(i))

over the probability simplex.

The high-level routine follows upstream preprocessing: normalize `w` and `x0`,
remove all-zero columns, optionally normalize rows, optionally replace `L` by a
truncated SVD reconstruction, and augment the likelihoods by the numerical
safeguard `eps` during optimization.

## EM update

The E step multiplies each likelihood column by its current mixture weight,
normalizes each row first by its maximum and then by its sum, and uses a
`1e-15` additive safeguard exactly as the C++ code does. The M step is the
weighted column sum of posterior assignment probabilities.

## SQP step

At each outer iteration the code performs one EM update, thresholds very small
mixture weights, evaluates the objective and its gradient/Hessian, and checks
the dual-residual convergence criterion. The quadratic subproblem is solved by
an active-set method based on Nocedal and Wright. If the reduced Hessian is not
positive definite, a multiple of the identity is increased geometrically until
Cholesky factorization succeeds. A feasibility-aware backtracking line search
then accepts the SQP step.

## SVD backend difference

Upstream `tsvd` repeatedly calls R `irlba`. The Fortran v0.1.0 implementation
uses LAPACK `DGESDD` and truncates the full decomposition according to the same
`tol_svd` rule. This is numerically compatible but can be slower for very large
low-rank matrices. It deliberately avoids copying GPL-licensed irlba code into
the MIT mixsqp distribution.
