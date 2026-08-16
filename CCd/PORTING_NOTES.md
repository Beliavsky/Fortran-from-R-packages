# Porting notes

## Scope

The supplied CCd 1.1 package contains seven computational functions and no
compiled or plotting code. All seven are represented in the Fortran port.
R-specific list/names/model-matrix infrastructure is replaced by derived
result types and an explicit numeric matrix interface.

## Upstream location/CDF limitation

The upstream `dcc` implementation uses

`lambda*tanh(pi*lambda) / (pi*(lambda^2 + (y-mu)^2))`.

That normalization is the exact lattice normalization when `mu` is an integer;
for noninteger `mu`, the normalizing constant depends on `mu`. Nevertheless,
upstream `cc.mle` and `cc.reg` optimize `mu` or linear predictors as continuous
real values while retaining the integer-location normalizer.

There is a second issue in upstream `pcc`: it computes the CDF using symmetry
about zero, even when `mu` is nonzero. Therefore `pcc(y, mu, lambda)` and
`qcc(..., mu, lambda)` should only be interpreted as genuine CDF/quantile
operations for the zero-location model. The Fortran implementation preserves
that upstream formula for compatibility; tests of CDF/quantile identities use
`mu=0`.

## Quantiles

The R implementation minimizes `abs(p - pcc(y,...))` over a large continuous
interval and then takes `ceiling()`. The Fortran implementation computes the
left quantile directly: the smallest integer `q` with `pcc(q) >= p`. This is
faster and deterministic while agreeing with the intended left-tail quantile
semantics for the valid zero-location CDF.

## Optimization

* `cc_mle0`: bounded golden-section maximization over approximately the same
  `(0,1000)` interval used by R `optimize`.
* `cc_mle` and `cc_reg`: standalone Nelder-Mead minimization, matching the
  default derivative-free character of R `optim`.
* `cc_reg` starts from an ordinary least-squares fit and the half-IQR scale,
  mirroring the upstream use of `Rfast::lmfit` and `Rfast::nth`.

No Rfast runtime dependency is required.
