# Porting notes

## Scope

This release translates the computational surface of `sadists` 0.2.6. The
only exported function omitted is `runExample`, which launches a Shiny user
interface and contains no distribution mathematics.

All plotting, R documentation machinery, and package/UI plumbing remain only
under `upstream/` for provenance.

## v0.2 dependency refactor

Version 0.1.0 carried private translations of the PDQutils routines needed by
sadists. Now that PDQutils itself has a complete standalone Fortran port,
v0.2.0 removes that compiled duplication.

`sadists_approximations` contains no expansion algorithm. It imports and
renames the canonical PDQutils-fortran routines:

- `dapx_edgeworth` -> `edgeworth_pdf`
- `papx_edgeworth` -> `edgeworth_cdf`
- `qapx_cf` -> `cornish_fisher_quantile`
- `as269` -> `as269`

The old internal names are retained only to keep the rest of the sadists port
and existing v0.1 client code source-compatible.

`sadists_special` also reuses PDQutils' normal PDF/CDF/quantile functions and
its canonical moment/cumulant conversion. The sadists-facing conversion
routines remain as thin subroutine wrappers because v0.1 exposed output-array
interfaces whereas PDQutils returns allocatable function results.

The dependency is vendored at `vendor/pdqutils-fortran` and referenced through
`fpm.toml`, allowing the archive to build without a network fetch.

## Approximation layer

The R package delegates most density, CDF, and quantile calculations to
PDQutils. The Fortran port now delegates these to PDQutils-fortran as well:

- Edgeworth density approximation
- Edgeworth CDF approximation
- Lee and Lin Algorithm AS269
- Cornish-Fisher quantile approximation
- raw moments <-> raw cumulants
- shared normal probability/quantile support

The approximation order defaults remain those of sadists: 6 for most
families, 5 for products of chi-square powers and products of normals, and 4
for products of doubly noncentral F variables.

The approximate density is clipped at zero and the approximate CDF to [0,1],
matching PDQutils. This does not make the CDF necessarily monotone.

## Remaining standalone replacements

The R package also imports `hypergeo` and `orthopolynom`; they remain absent as
runtime dependencies.

- Normal raw moments use a direct recurrence rather than orthopolynom objects.
- Noncentral chi-square power moments use the Poisson-mixture representation
  and gamma ratios rather than a confluent-hypergeometric call.

The small PDQutils R-source snapshot retained in `upstream/PDQutils/` is
provenance material only and is not compiled. The canonical compiled
implementation is the vendored PDQutils-fortran dependency.

## Noncentral log-chi-square moments

Upstream `lnc_moments` sums Poisson-mixture terms for indices 0:50. This port
uses an adaptive interval centered on the Poisson mean and renormalizes the
retained weights. This is more reliable for large noncentrality parameters
while converging to the same calculation in ordinary cases.

## Tail and endpoint handling

Stable Cornish-Fisher endpoint and normal-tail behavior now comes from
PDQutils-fortran. Consequently transformed beta and eta quantiles return their
finite support endpoints without evaluating high-order corrections at infinite
normal quantiles.

## RNG implementation

Sadists-specific simulation remains local. Normal random numbers use
Box-Muller. Gamma random numbers use Marsaglia-Tsang, with shape raising below
one. Noncentral chi-square generation uses a Poisson mixture of central
chi-squares. Large-Poisson generation uses transformed rejection.

The sequence is not intended to match R's RNG bit-for-bit. Distributional
behavior is tested instead.

## Vectorization

R automatically vectorizes scalar `d/p/q` arguments. The Fortran public
functions are scalar. `_vec` subroutines are supplied for doubly noncentral F,
t, beta, and eta. Random generators fill caller-provided arrays.

Multi-component parameter vectors preserve R-style recycling among parameter
arrays.

## Regression protection

`test_pdqutils_integration` constructs a representative weighted chi-square
sum, obtains its sadists cumulants, and verifies that the sadists density, CDF,
and quantile equal direct calls to PDQutils-fortran with those cumulants. It
also compares both directions of the v0.1 moment/cumulant compatibility API
against PDQutils.
