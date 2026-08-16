# Porting notes

## Scope

This is a complete computational port of the supplied R package.  The source
contains no plotting routines or compiled native code.

## Estimation

`estdlaplace2(method="M")` uses bounded optimization in R.  The Fortran port
uses a logistic transformation and Nelder-Mead, enforcing the same
`[err,1-err]` bounds.

Upstream `estdlaplace2(method="ML")` invokes unconstrained `optim`, although a
valid model requires `0 < p,q < 1`.  The Fortran port optimizes the logits of
`p` and `q`, so every trial point is a valid distribution.  The likelihood
being optimized is unchanged.

The `P` and `MM` estimators are direct translations of the R formulas.

## Quantile endpoints

The distributions have unbounded integer support.  For probabilities exactly
0 or 1, Fortran quantile functions return `-huge(int64)` or `huge(int64)` as
sentinels.  For interior probabilities the upstream step convention is
preserved, including its strict comparison at CDF jump points.

## Information matrices

`ifi` is already the inverse Fisher-information matrix in the upstream code.
`ifi2` and `iofi2` form the same 2x2 matrices and use an explicit closed-form
2x2 inverse rather than R's generic `solve`.

## RNG

Random generation uses Fortran's intrinsic uniform generator followed by the
same quantile construction.  Therefore distributions agree with R, but seeded
streams are not intended to match R's RNG bit-for-bit.
