# Porting notes

## License

The upstream DESCRIPTION declares `GPL-3`; the Fortran project uses the SPDX
identifier `GPL-3.0-only`. The complete GPLv3 text is included in `LICENSE`.
Every translated Fortran source carries a matching SPDX header.

## Numerical stability

- HMM and HHMM likelihoods use log-sum-exp recursions.
- Forward-backward probabilities use per-period scaling.
- Probabilities entering logarithms are bounded below by `tiny_prob`.
- The stationary distribution is solved from the linear stationarity equations
  and falls back to a uniform distribution if the system is singular.
- Gamma and incomplete-beta calculations use convergent series and continued
  fractions.
- Quantile calculations use closed forms where possible and safeguarded
  bisection otherwise.
- Transition rows are constructed through a multinomial-logit representation.

## Parameter ordering

The upstream package extracts off-diagonal transition probabilities in R's
column-major order. The Fortran pack/unpack routines preserve that ordering:
columns are traversed first and diagonal entries are skipped.

The packed ordinary-HMM vector is:

1. transition logits,
2. means, log-linked for gamma and Poisson models,
3. log standard deviations except for Poisson models,
4. log degrees of freedom for Student-t models.

An HHMM vector contains the coarse HMM block followed by one fine-HMM block for
each coarse state.

## Forecast correction

The upstream prediction routine calculates weighted averages of component
quantiles. A mixture quantile must instead solve the mixture CDF. The Fortran
routine uses the mathematically correct mixture quantile by default and exposes
`upstream_quantiles=.true.` for compatibility.

## Hierarchical decoding

Coarse Viterbi emissions combine the coarse observation density with the entire
fine-chunk HMM likelihood conditional on each coarse state. Fine paths are then
decoded under the selected coarse state. This follows the package likelihood
factorization.

## Optimization

`fit_hmm` and `fit_hhmm` replace R's `nlm`, `pracma::hessian`, and
`MASS::ginv` with:

- native multi-start Nelder-Mead,
- central finite-difference gradients and Hessians,
- pivoted linear-system inversion.

The objective context is stored at module scope during fitting, so concurrent
fit calls must be serialized. Other library routines are reentrant.

## Discrete residuals

As upstream, Poisson pseudo-residuals transform the ordinary CDF through the
normal quantile. Randomized probability-integral-transform residuals are not
used. Users requiring randomized discrete residuals can add a uniform draw
inside the probability mass at each observation.

## Missing hierarchical cells

Fortran uses `chunk_lengths` to identify active fine observations. Padding
outside the active range may be NaN or any other value because it is ignored.
This is safer than relying on missing-value scanning inside the likelihood.
