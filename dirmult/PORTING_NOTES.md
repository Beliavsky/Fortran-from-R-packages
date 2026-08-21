# Porting notes

## Numerical model

The port retains the upstream gamma parameterization

`gamma_plus = (1-theta)/theta`, `gamma_j = pi_j * gamma_plus`.

As in the R code, `dirmult_loglik` contains only the parameter-dependent part
of the Dirichlet-multinomial log-likelihood and omits the multinomial
coefficient. This constant omission has no effect on parameter estimates or
profile-likelihood differences and preserves the values returned by `dirmult`.

The score and observed Fisher-information routines follow the upstream finite
reciprocal sums. Rising-log products are evaluated with `log_gamma` for better
clarity and numerical range. The expected-information calculation retains the
upstream beta-binomial tail-probability formula, evaluating beta-binomial
masses in log-gamma form.

## Linear algebra

The R package uses `solve`. This port contains a small partial-pivoting dense
linear solver and matrix inverse, avoiding an external BLAS/LAPACK dependency.
`fit_equal_theta` solves the full KKT Newton system directly rather than
constructing the block inverse used by the R implementation. The systems are
algebraically equivalent, while the direct solve avoids unnecessary explicit
inversion.

## Profile-score compatibility

Upstream `profU` subtracts the fixed `gamma_plus` value from every gamma-score
component, rather than the current Lagrange multiplier. Under the fixed-sum
constraint a common shift of all gamma scores is absorbed by the multiplier,
so the gamma/pi Newton step and profile log-likelihood are unchanged; only the
reported multiplier is shifted. `estimate_profile_loglik` deliberately keeps
this convention for compatibility with dirmult 0.1.3-5.

The upstream profile code also temporarily clamps negative gamma values without
writing those clamps back into its combined gamma/lambda state vector. The port
retains this behavior in the profile path. Normal convergent cases do not hit
this branch.

## Zero rows and columns

`dirmult` itself removes zero rows and zero columns. The help page for
`weirMoM` says the same even though the R function body does not explicitly do
so. The Fortran `weir_mom` follows the documented behavior and removes both,
preventing division by zero on empty rows.

## Iteration safety

The R `dirmult` routine has no maximum-iteration argument. `fit_dirmult` adds
an optional `maxit` with default 1000 so a pathological fit cannot loop
forever. The ordinary convergence rule and the upstream negative-gamma clamp
(`0.01`) are preserved.

## Random generation

R calls `rgamma`, `rnorm`, and `rmultinom`. This port uses Fortran
`random_number`, Box-Muller normal generation, Marsaglia-Tsang gamma generation,
and categorical multinomial draws. Distributional behavior is preserved, but
random streams are not bit-for-bit compatible with R. Use `seed_rng` or the
optional `seed` arguments for reproducible Fortran runs.

`null_test` can retain all simulated tables in one `J x K x m` integer array.
Set `store_data=.false.` when only the statistics are needed.

## Data structures

R lists/data frames/names are represented by typed Fortran derived types:
`dirmult_fit_type`, `dirmult_summary_type`, `profile_fit_type`,
`profile_grid_type`, `equal_theta_fit_type`, `sim_pop_result_type`, and
`null_test_result_type`. For a common-theta fit, input tables use
`count_table_type` because different tables may have different category counts.
