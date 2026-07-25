# bayesGARCH Modern Fortran

A modern Fortran translation of the computational core of the R package
`bayesGARCH` 2.1.10. It implements Bayesian estimation of a GARCH(1,1) model
with variance-standardized Student-t innovations using the package's tailored
Metropolis-Hastings construction.

This is a source translation, not an R compatibility layer. It does not claim
bit-for-bit random-number equivalence or numerical identity of complete Markov
chains with R.

## Implemented and tested

- GARCH(1,1) conditional-variance filtering and Student-t log likelihood.
- General GARCH(p,q) and sign-dependent threshold-GARCH filtering.
- General GARCH and threshold-GARCH simulation from supplied innovations.
- Variance-standardized Student-t GARCH(1,1) simulation.
- Latent inverse-gamma scale draws for the Student-t mixture representation.
- Tailored Metropolis-Hastings updates for `(alpha0, alpha1)` and `beta`.
- Acceptance-rejection update for Student-t degrees of freedom `nu`.
- Multiple chains, independent starting-value columns, optional stationarity,
  and an optional user procedure for additional prior constraints.
- Original symmetric and asymmetric alpha and ARMA-style W filters.
- Quasi-differencing and one- or two-parameter weighted Gaussian posterior
  calculations used by the proposal distributions.
- Posterior burn-in removal, thinning, chain concatenation, means, and standard
  deviations without R's `coda` classes.
- Normal, gamma, exponential, and standardized Student-t random generation used
  by the sampler.

The regression tests include hard-coded independently computed values for the
GARCH recursion, likelihood, proposal regressions, augmented posterior, and
`nu` proposal equation. Fixed-seed stochastic tests check random-generator
moments, simulation/filter consistency, movement of every MCMC block,
stationarity and custom constraints, and posterior-sample indexing.

## Not implemented

- R `mcmc` and `mcmc.list` classes or `coda` diagnostics.
- S3 methods, R object attributes, R callback/list semantics, or R package
  registration.
- Plotting.
- Exact R RNG streams or bit-for-bit R chain reproduction.
- Parallel-chain execution; chains are currently run sequentially.

No claim of full statistical equivalence to the R package is made. R was not
available in the validation environment, so complete chain-by-chain comparisons
against the original package were not performed.

## Build and test

GNU Fortran:

```sh
make
make check
make optimized-check
```

The checked build uses Fortran 2018, warnings as errors, bounds/runtime checks,
and backtraces. `fpm.toml` is also supplied:

```sh
fpm test
fpm run demo_bayesgarch
fpm run fit_csv -- data/example_returns.csv 2000 500 5 1
```

`fpm` was not available in the validation environment, so only the GNU Make
build is claimed as tested.

## Minimal use

```fortran
use bayesgarch_kinds, only : dp
use bayesgarch_sampler, only : bayesgarch_control, bayesgarch_result, run_bayesgarch

real(dp) :: y(1000)
type(bayesgarch_control) :: control
type(bayesgarch_result) :: fit

control = bayesgarch_control()
control%n_chains = 2
control%n_iter = 10000
control%enforce_stationarity = .true.
call run_bayesgarch(y, fit, control=control)
```

Draws are stored as `fit%draws(iteration, parameter, chain)` in parameter order
`alpha0`, `alpha1`, `beta`, `nu`.

## Command-line CSV fitting

`fit_csv` reads the first field of each line as a real value, skips lines whose
first field is nonnumeric, and accepts:

```text
fit_csv FILE [N_ITER=2000] [BURN=500] [THIN=5] [N_CHAINS=1]
```

## License and citation

The original package declares `GPL (>= 2)`. This derivative is therefore
licensed under `GPL-2.0-or-later`; see `LICENSE`, `COPYRIGHTS`, and `ORIGIN.md`.

Users should cite the original work:

- David Ardia and Lennart F. Hoogerheide (2010), "Bayesian Estimation of the
  GARCH(1,1) Model with Student-t Innovations," *The R Journal*, 2(2), 41-47.
- David Ardia (2008), *Financial Risk Management with Bayesian Estimation of
  GARCH Models*, Springer.
