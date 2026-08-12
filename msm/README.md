# msm-fortran

`msm-fortran` is a modern Fortran translation of the computational core of the
R package **msm 1.8.2: Multi-State Markov and Hidden Markov Models in Continuous
Time**.

The project is intended as a numerical library rather than an emulation of R's
formula, data-frame, S3, and plotting interfaces.  It uses allocatable arrays,
explicit derived types for emission models, 1-based state numbers, and a small
public convenience module named `msm`.

## Implemented computational functionality

### Continuous-time Markov models

- construction and validation of generator matrices;
- `P(t) = exp(t Q)` using scaling, squaring, and Pade approximation;
- panel-observation transition probabilities;
- exact-transition kernels used by `msm(exacttimes=TRUE)`;
- exact-death likelihood contributions;
- Frechet derivatives of transition matrices with respect to arbitrary `dQ`;
- ordinary sequence likelihoods and analytic Q gradients;
- aggregate transition-count likelihoods corresponding to the fast ordinary
  likelihood path in `src/lik.c`;
- non-hidden censoring likelihoods using sets of possible true states;
- piecewise-constant Q matrices;
- log-linear covariate effects on individual transition intensities.

### Hidden Markov models

- normalized forward likelihood recursion;
- backward smoothing probabilities;
- Viterbi most-likely paths;
- multivariate outcomes with one emission model per outcome/state;
- missing outcome components through IEEE NaNs;
- optional known true states and allowed-state masks;
- panel, exact, and exact-death observation kernels;
- analytic score with respect to arbitrary Q perturbations.

All 17 emission families in `src/hmm.c` are represented:

`categorical`, `identity`, `uniform`, `normal`, `lognormal`, `exponential`,
`gamma`, `Weibull`, `Poisson`, `binomial`, `beta-binomial`, `truncated normal`,
`truncated normal plus measurement error`, `uniform plus measurement error`,
`negative binomial`, `beta`, and location/scale Student t.

The analytic emission derivatives from `src/hmmderiv.c` are also translated.
As in the original package, derivatives of the truncated-normal and
measurement-error families are not supplied.

### Multi-state summaries and distributions

The library includes mean sojourn times, next-state probabilities, transient
and absorbing-state detection, finite and eventual first-passage
probabilities, expected first-passage times, finite/infinite expected state
occupancy, expected visit counts, and state prevalence.  It also translates
the piecewise-exponential, two-phase, and truncated-normal helpers used by the
R package.

## Numerical consolidation relative to the R package

`src/analyticp.c` in msm contains hand-derived transition-probability formulas
for a collection of two- through five-state Q-matrix topologies.  These are
performance special cases for the same mathematical object `exp(t Q)`.  This
Fortran release deliberately consolidates them into one general, robust
scaling-and-squaring implementation instead of reproducing about 1,300 lines
of topology-specific formulas.  The regression tests include closed-form
2-state comparisons and derivative checks.

Likewise, the original C derivative path switches between eigendecomposition
and a truncated power series when eigenvalues repeat.  The Fortran port uses
the standard block-exponential Frechet identity for every matrix, which is
stable for repeated and nearly repeated eigenvalues and computes the same
mathematical derivative.

## Not translated

The following are R infrastructure or external services rather than the native
numerical engine and are intentionally omitted:

- formula/data-frame/model-frame parsing and S3 methods;
- plotting, printing, `tidy` methods, and survival-object adapters;
- calls to R's `optim`, `optimHess`, `minqa`, `mvtnorm`, `survival`, and
  parallel/foreach infrastructure;
- bootstrap orchestration and R-specific confidence-interval presentation;
- the hand-coded analytic-P topology dispatcher noted above.

The library exposes likelihoods, gradients, transition matrices, and summary
statistics so a Fortran optimization package can be connected directly.

## Building with FPM

```text
fpm build
fpm test
fpm run --example basic_ctmc
fpm run --example hidden_markov
```

BLAS and LAPACK are linked through the FPM manifest.  The only LAPACK routine
required by the current core is `dgesv`.

## Minimal example

```fortran
program demo
    use msm, only : dp, make_generator, transition_matrix
    implicit none
    real(dp) :: off(3,3), q(3,3)
    real(dp), allocatable :: p(:,:)

    off = 0.0_dp
    off(1,2) = 0.15_dp
    off(1,3) = 0.05_dp
    off(2,3) = 0.12_dp
    q = make_generator(off)
    p = transition_matrix(q, 5.0_dp)
    print *, p
end program demo
```

## Validation

The seven-test regression suite checks closed-form two-state transition probabilities,
Frechet derivatives against finite differences (including a repeated-eigenvalue generator),
panel/exact/death likelihoods,
aggregate and censored likelihoods, HMM likelihoods against explicit path
enumeration, Viterbi paths, Q-score derivatives, all emission families,
distribution identities, simulation, multi-state summary quantities (including
competing absorbing classes and non-certain passage), and the delta method.

The release was validated with GNU Fortran 14.2.0 using both bounds-checked
`-O0` and optimized `-O2` builds.  FPM itself was not installed in the build
environment, so the FPM source layout was validated with equivalent direct
`gfortran` compilation and BLAS/LAPACK linking.

## License

GPL-2.0-or-later.  See `LICENSE`, `LICENSES.md`, and `licenses/`.
