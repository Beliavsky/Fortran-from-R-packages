# Porting notes

## Scope

The upstream package has no compiled native source; its computational code is
in `R/dpqr-stable.R` and `R/dist-stableMode.R`. Those numerical routines were
translated. `R/utils.R` contains an R-specific root wrapper/check helper rather
than a statistical algorithm and is represented only by the small internal
bracketed root solver needed by `qstable` and the Nolan peak searches.

## Numerical integration

Upstream `dstable`/`pstable` call R's adaptive `integrate()` (QUADPACK-backed).
The supplied `r_mod.f90` already provides the requested R-compatible
`integrate` helper, so this port uses it rather than introducing another
quadrature implementation. Difficult Nolan intervals are split into several
finite panels before calling `r_mod::integrate`; density integrals retain the
upstream peak-splitting idea around `g(theta)=1`.

This means internal subdivision behavior is not bit-for-bit identical to R's
QUADPACK implementation even when the statistical result agrees closely.

## Quantiles

Upstream `qstable` builds Normal/Cauchy-based initial brackets and uses
`uniroot`. The port uses monotone CDF bracketing with interval expansion and
bisection. This avoids adding a second general optimizer/root package and gives
stable p/q inversion while preserving the public parameterization semantics.

One historical upstream test records
`qstable(.6, alpha=.5, beta=1) = 2.636426573120147`. For this Levy case, the
closed-form identity gives approximately `2.6364178821` in S0. The Fortran port
uses the exact Levy identity already supplied in `inst/xtraR/Levy.R`, so it
follows the closed form rather than reproducing that numerical-integration
error.

## Random generation at alpha=1

The R source evaluates the general Chambers-Mallows-Stuck expression even at
`alpha=1` (except the symmetric Cauchy case), which involves a numerically huge
`tan(pi/2)` and cancellation. The port evaluates the mathematically equivalent
`alpha=1` limiting expression explicitly. Monte Carlo tests verify the result
against the translated S0/S1 CDFs for non-unit scale and nonzero skew.

## S2 and stable mode

`pm=2` uses the same upstream rescaling by `alpha**(-1/alpha)` and translates
the distribution so that `delta` is the mode. The supplied helper module has
no bounded scalar `optimize` equivalent; the port adds a small golden-section
maximizer solely for `stable_mode`.

## Known upstream numerical difficulty

The upstream manual explicitly warns that `alpha` close to 1 or 0 can be
considerably harder numerically and may lose accuracy. The port preserves the
same Nolan representation and does not claim uniform precision in those
regions. Exact/special-case reductions are used where available.

## r_mod

`upstream/r_mod-original.f90` is the exact supplied helper. `src/r_mod.F90` is
a line-wrapped build copy. After removing whitespace and continuation markers,
the two files are lexically identical; no helper implementation was changed.
