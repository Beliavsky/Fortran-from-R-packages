# Porting notes

## Source mapping

- `R/pearson0.r` -> type 0 d/p/q/r
- `R/pearsonI.R` -> type I d/p/q/r
- `R/pearsonII.r` -> type II d/p/q/r
- `R/pearsonIII.r` -> type III d/p/q/r
- `R/pearsonIV.r` + `src/pearsonIV.c` -> type IV d/p/q/r, normalization, rejection RNG
- `R/pearsonV.r` -> type V d/p/q/r
- `R/pearsonVI.r` -> type VI d/p/q/r
- `R/pearsonVII.r` -> type VII d/p/q/r
- `R/pearson.r` -> generic dispatch and generic moments
- `R/moments.r` -> analytical moment functions
- `R/fit.r` -> empirical moments, moment matching, family MLEs, all-family MLE, MSC
- `R/matchMoments.R` -> `match_moments`
- `R/plots.r` -> intentionally skipped
- `R/zzz.r` and `src/pearson.c` -> R package registration/UI glue, intentionally skipped

## Pearson IV

The public R function `ppearsonIV` selects `.ppearsonIVint` whenever GSL is not
available. This translation deliberately uses that no-GSL numerical-integration
path, with `r_mod`'s integration helper, so no GSL dependency is needed.

The type-IV log normalization and rejection sampler are direct formula/algorithm
ports of `src/pearsonIV.c`, including its Joel Heinrich provenance.

The upstream C tree also contains double-double and quad-double arithmetic used
to improve the internal `F21` hypergeometric series in some GSL-enabled CDF
branches. Those extended-precision accelerator internals are not on the chosen
public CDF path and are not ported in v0.1.0. A normal double-complex
`hypergeom_2f1` series is included for ordinary |z|<1 use.

## Maximum likelihood

PearsonDS uses R's `nlminb` on the direct family parameters. A literal use of a
finite-difference bound optimizer is fragile for families whose support depends
on fitted location/scale parameters. The Fortran port therefore optimizes smooth
transformed parameters that enforce:

- positive shape parameters,
- positive Pearson-IV scale,
- bounded-support endpoints outside the sample for types I/II,
- one-sided support outside the sample for types III/V/VI,
- positive canonical scale for the symmetric type VII.

The decoded likelihood and distribution are the same family likelihoods. Type I
and II MLEs can be mathematically unbounded for some data/configurations as a
support endpoint approaches an observation; this is an inherent likelihood
issue rather than a Fortran-specific error.

## Source behavior retained explicitly

`matchMoments.R` in PearsonDS 1.3.2 has a type-V branch for missing skewness that
references `sss` while `sss` is `NA`. Rather than silently inventing a corrected
formula, the Fortran `match_moments` returns a nonzero status for that branch.
