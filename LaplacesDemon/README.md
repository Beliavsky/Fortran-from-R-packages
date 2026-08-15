# LaplacesDemon-fortran

Modern Fortran 2018/FPM translation of the computational core of the R
package **LaplacesDemon 16.1.8**.

The upstream package is a large pure-R Bayesian inference environment. This
port exposes its reusable numerical algorithms through explicit procedure
callbacks and numeric arrays. Plotting, printing, R classes/lists, formula
handling, monitor/storage machinery and cluster/HPC orchestration are
intentionally outside the Fortran API.

## v0.3.1 test-portability fix

v0.3.1 contains the same numerical algorithms as v0.3.0. It makes the
`test_v02` stochastic assertions portable across compiler/runtime-specific
Fortran `random_number` sequences and prints the name of any failed assertion.
The revised thresholds were stress-tested over 50 alternate seeds with no
failures.

## v0.3.0 highlights

v0.3.0 is the broad numerical-completion release. In addition to all v0.1 and
v0.2 functionality it adds the remaining fixed-dimensional sampler families,
the full optimizer catalog used by `LaplaceApproximation`, sparse/componentwise
iterative quadrature, additional convergence diagnostics and a substantially
expanded probability/prior library.

### Optimizers and Laplace-approximation backends

The public optimizer layer now contains counterparts for the upstream catalog:

- AGA / genetic search;
- BFGS and BHHH;
- nonlinear conjugate gradient and DFP;
- hit-and-run search and Hooke-Jeeves;
- L-BFGS and Levenberg-Marquardt;
- Nelder-Mead and Newton-Raphson;
- particle swarm optimization;
- resilient backpropagation;
- stochastic gradient descent;
- SOMA;
- spectral projected gradient;
- SR1;
- trust-region optimization.

### MCMC

Together with the earlier RWM, MWG, AM, MALA, HMC, HMCDA, NUTS, DRAM, RAM,
pCN, t-walk, slice, ESS, AIES and DEMC implementations, v0.3 adds numerical
counterparts for:

- ADMG, AFSS and adaptive/griddy Gibbs;
- AHMC;
- adaptive-mixture Metropolis and adaptive MWG;
- CHARM and HARM;
- delayed-rejection Metropolis;
- independence and multiple-try Metropolis;
- in-process INCA adaptation;
- MCMCMC;
- OHSS and UESS;
- RDMH using the upstream `U^s` random-dive proposal and Jacobian;
- refractive sampling;
- RSS;
- SAMWG/SMWG and their updating variants;
- SGLD;
- tempered HMC;
- callback Gibbs sampling;
- fixed-vector reversible-jump variable selection.

The upstream `Experimental` branch is not a stable algorithm and is not
represented as an API routine. Network/cluster transport used by some R
variants is replaced by in-process numeric-array implementations.

### Quadrature and diagnostics

- componentwise iterative quadrature;
- adaptive sparse-grid quadrature;
- BMK diagnostic;
- Heidelberger-Welch diagnostic;
- KS diagnostic;
- the v0.2 Raftery-Lewis and Hangartner diagnostics remain available.

### Expanded probability library

The distribution/prior layer now includes the covariance, precision and
Cholesky parameterizations used by LaplacesDemon for multivariate normal,
Student-t, Cauchy, Laplace and power-exponential laws; Wishart and
inverse-Wishart RNGs and Cholesky forms; matrix gamma/inverse-matrix-gamma;
normal-Laplace and mixture priors; Huang-Wand, horseshoe, LASSO,
normal-Wishart/normal-inverse-Wishart, Yang-Berger, hyper-g and Zellner priors;
continuous-relaxation MRF; generic truncation density/CDF/quantile/RNG plus
truncated moments; and the scalar compatibility CDF/quantile/RNG APIs used by
the upstream package.

Use the umbrella module:

```fortran
use laplacesdemon
```

## Callback interface

Most algorithms use a fixed-dimensional log-target callback:

```fortran
module target_mod
  use laplacesdemon, only: dp
contains
  function log_posterior(theta) result(lp)
    real(dp), intent(in) :: theta(:)
    real(dp) :: lp
    lp = -0.5_dp * sum(theta**2)
  end function
end module
```

For example:

```fortran
call nuts_sample(log_posterior, x0, 5000, 1000, 2, result,
                 adapt_steps=1000, target_accept=0.7_dp)
```

## Fidelity boundary

The numerical algorithms have native Fortran counterparts, but the port does
not reproduce R's dynamic `Model(theta, Data)` list protocol, plotting,
printing, formula evaluation, monitor objects, or parallel socket/cluster
plumbing. BHHH therefore receives an explicit per-observation component-score
callback, Gibbs receives an explicit conditional-draw callback, and the RJ
API represents the package's fixed-vector variable-selection semantics with an
active mask rather than dynamically changing Fortran array rank/size.

## Validation

The source is validated with GNU Fortran using:

```text
-std=f2018 -Werror=implicit-interface -Werror=trampolines -fcheck=all -O0
```

Fortran does not standardize the exact sequence produced by the intrinsic
`random_number` generator. Stochastic recovery tests therefore use statistical
tolerances rather than assuming the seed maps to the same stream on every
compiler.

Permanent tests exercise all historical v0.1/v0.2 functionality plus optimizer
recovery, the completed sampler catalog, sparse/componentwise quadrature,
additional convergence diagnostics, and distribution/parameterization
identities. See `docs/API_MAP.md` for the detailed mapping.
