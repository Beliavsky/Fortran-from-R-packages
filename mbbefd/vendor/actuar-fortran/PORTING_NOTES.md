# Porting notes

## Design goals

The Fortran port prioritizes:

1. self-contained numerical algorithms;
2. typed, explicit interfaces;
3. deterministic optional random seeds;
4. standard Fortran 2018 portability;
5. clear failure diagnostics for structured results;
6. preservation of upstream licensing and attribution.

The R package is much broader than this release. The code under `original/` is
retained so users can inspect exact upstream behavior when a routine has not
been translated.

## Scalar APIs instead of R vectorization

R probability routines are vectorized and support recycling, `log`,
`lower.tail`, and related conventions. The Fortran functions generally operate
on scalar values and return ordinary lower-tail probabilities. Array operations
can be expressed naturally with loops, array constructors, or `elemental`
wrappers in caller code.

Invalid scalar parameters usually return IEEE NaN. Algorithms with several
outputs return a typed result containing `ok` and `message` fields.

## Distribution parameterization

The names and primary parameterizations follow the upstream actuarial families,
but Fortran argument names are explicit (`shape`, `scale`, `rate`, and so on).
Users should not assume another package's similarly named distribution uses the
same convention.

Boundary points are handled through explicit branches. The validation scripts
therefore suppress GNU Fortran's `-Wcompare-reals` warning for intentional
comparisons against mathematical endpoints such as zero and one.

## Special functions

The upstream package uses mature R math-library and C implementations and links
to `expint`. This port instead contains native implementations of:

- normal distribution functions;
- regularized incomplete gamma and beta functions;
- beta/gamma quantile inversion;
- adaptive Simpson quadrature;
- a direct integral representation of the modified Bessel K function.

These are sufficient for the implemented families and tests, but they are not a
replacement for a full arbitrary-precision special-functions library. Extreme
parameter combinations may require tighter tolerances or a specialist library.

## Random generation

Random generation uses a small self-contained xorshift generator and standard
transformation/rejection algorithms. Seeded Fortran sequences are deterministic
within this project but do not match R's RNG streams.

## Poisson-inverse Gaussian

The PIG random generator uses inverse-Gaussian mixing followed by a Poisson
draw. PMF calculations use the Bessel-K integral supplied by `actuar_special`.
This avoids external dependencies but is slower than a specialized Bessel
implementation for large-scale evaluation.

## Aggregate distributions

The aggregate module works with probability vectors on an equally spaced
nonnegative grid. The `step` value maps an integer index to a monetary loss
amount. Truncation is controlled by the caller through maximum support sizes.
Probability mass remaining above the truncated support is not silently
redistributed.

`discretize_cdf` accepts a procedure callback rather than an R expression or
function object. `compound_simulation` similarly accepts typed frequency and
severity callbacks.

## Panjer recursion

The generic recursion follows the `(a,b,0)` class. Convenience wrappers construct
frequency parameters for Poisson, binomial, and negative-binomial counts. The
result reports failure rather than continuing when the starting denominator is
singular or probability inputs are invalid.

## Aggregate CTE

For a discrete distribution, CTE is computed from the quantile threshold and
includes the threshold atom consistently with the returned discrete VaR. This
may differ from definitions that interpolate within a probability atom.

## Phase-type calculations

A scaling-and-squaring matrix exponential with a Taylor kernel is used. Raw
moments solve repeated linear systems involving the subintensity matrix. This is
appropriate for modest phase dimensions but is not optimized for very large or
highly ill-conditioned matrices.

## Credibility

Buhlmann-Straub uses weighted class means and method-of-moments process and
structural variance estimates. Negative estimated structural variance is
truncated at zero before credibility factors are formed.

The conjugate Bayesian routines expose the posterior credibility combination
directly rather than reproducing R list structures and printing methods.

## Coverage transformation

`apply_coverage` treats `limit` as the maximum gross covered loss before
coinsurance. With an ordinary deductible `d` and limit `u`, the maximum payment
before coinsurance is therefore `max(u-d,0)`. A franchise deductible pays the
covered amount in full once the threshold is exceeded. Inflation is applied
before deductible and limit calculations.

## Ruin calculations

The explicit exponential-claim formula is fully tested. The phase-type route is
a focused Cramer-Lundberg implementation intended for modest transient phase
systems. It should not be represented as the complete upstream `ruin` framework,
which supports additional renewal structures and numerical strategies.

## Excluded estimation frameworks

The current release does not translate the generic minimum-distance (`mde`) and
empirical-moment (`emm`) estimation frameworks. Those routines depend heavily on
R function dispatch, optimization conventions, and distribution metadata. They
can be added later as typed Fortran optimization APIs without changing the
current distribution modules.

## Thread safety

The probability functions are pure where practical. The shared RNG state is not
thread-local; concurrent seeded simulation should use external synchronization
or separate process-level streams.
