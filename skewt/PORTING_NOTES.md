# Porting notes

## Scope

The upstream package has exactly four computational entry points: `dskt`,
`pskt`, `qskt`, and `rskt`.  All four are translated.  There is no fitting,
plotting implementation, compiled code, or class infrastructure in the
supplied package.

## Student-t dependency

Upstream delegates ordinary Student-t d/p/q operations to R's `stats` package.
The Fortran port is self-contained:

- density uses the log-gamma closed form;
- CDF uses the regularized incomplete beta identity;
- quantiles use monotone bisection with adaptive bracketing;
- a Student-t gamma/normal RNG is provided internally, although public `rskt`
  deliberately follows upstream and uses inverse transform through `qskt`.

## Parameter validation

The Fernandez-Steel skew parameter must satisfy `gamma > 0`; `df > 0` is also
required.  The R source does not explicitly validate these arguments.  The
Fortran routines return IEEE NaN for invalid parameters and for probabilities
outside `[0,1]`, and return infinities for exact endpoint quantiles.

## Distribution identity

For `gamma = 1`, the implementation reduces to the ordinary Student-t law.
For general `gamma`,

    P(X < 0) = 1 / (1 + gamma**2),

which is used as a direct parity test of the Fernandez-Steel construction.
