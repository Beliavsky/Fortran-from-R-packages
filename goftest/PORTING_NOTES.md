# Porting notes

## Scope

The computational exports of `goftest` 1.2-3 are translated:

- `pAD`, `qAD`
- `pCvM`, `qCvM`
- `ad.test`
- `cvm.test`
- `recogniseCdf`

The internal `simpleADtest`, `simpleCvMtest`, and Braun composite-null procedure
are also translated because they are essential to the exported tests.

R-specific `htest` object construction, function-name lookup through R
namespaces, expression deparsing, and printed method strings are intentionally
not reproduced. Their numerical counterparts are exposed through a Fortran
`gof_result` type and CDF callbacks.

## Anderson-Darling

The Marsaglia `adinf`, finite-sample correction, and high-accuracy `ADinf`
algorithms were translated from the package's C sources. The original comments
explicitly allow replacing the hand-coded complementary normal integral with
`erfc`; the Fortran port does this using the standard intrinsic.

The Marsaglia finite-sample correction can return a slightly negative number in
the extreme lower tail for small `n` (for example, near `n=10, q=0.1`). Because
a public CDF must lie in `[0,1]`, the Fortran public routine clamps this numerical
artifact to zero. This is a deliberate numerical correction and can therefore
differ from the literal R/C result in that tiny lower-tail region.

## Cramer-von Mises

The Csorgo-Faraway first-order finite-sample expansion and asymptotic series are
ported directly. R obtains fractional-order `besselK` from its math library;
the Fortran port evaluates

    K_nu(z) = integral_0^infinity exp(-z*cosh(t))*cosh(nu*t) dt

with composite 16-point Gauss-Legendre quadrature, switching to a standard
large-argument asymptotic expansion. The recurrence relation is used for
`K_(5/4)`.

Independent reference tests against SciPy's `scipy.special.kv` and a direct
Python transcription of the Csorgo-Faraway formulas agree to approximately
1e-10 or better over representative probability ranges.

## Quantiles

The R package uses `uniroot`. The Fortran routines use a bracket expansion plus
100 bisection iterations. This is slower than a safeguarded Newton method but is
simple, deterministic, and robust for the one-dimensional monotone CDFs.

## Braun adjustment

The R implementation randomly permutes nearly balanced group labels. The
Fortran routine does the same with `random_number` and Fisher-Yates shuffling.
For reproducible tests or applications, callers may instead supply explicit
`groups(:)` labels.
