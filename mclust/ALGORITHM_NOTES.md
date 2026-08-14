# Algorithm notes

## Covariance model family

The port supports the full mclust Gaussian covariance family.  With covariance
decomposition `Sigma_k = lambda_k D_k A_k D_k'`, the model letters encode
whether volume (`lambda`), shape (`A`), and orientation (`D`) are equal or
variable across mixture components.

Supported multivariate models are:

`EII VII EEI VEI EVI VVI EEE EVE VEE VVE EEV VEV EVV VVV`.

The univariate `E` and `V` models are also supported.

## Hybrid modernization

mclust already contains a large, mature legacy Fortran numerical core. In this port the upstream fixed-form sources are converted to free source form while preserving valid legacy language constructs such as arithmetic IF and labeled DO loops.
For the difficult constrained covariance M-steps, v0.1.2 deliberately retains
those algorithms instead of rewriting them into generic full-covariance EM.
Modern free-form modules provide:

- explicit procedure interfaces;
- derived fit/control types;
- allocation and shape checking;
- a common Gaussian-mixture E-step;
- covariance reconstruction as ordinary `Sigma(d,d,G)` arrays;
- hierarchical initialization and transformed data initialization;
- model selection and high-level statistical workflows.

This design makes the application API modern while preserving the specialized
MCLUST numerical algorithms.

## EM

The high-level `fit_model` loop alternates:

1. posterior evaluation using stable log densities and `logsumexp`;
2. the model-specific M-step;
3. relative log-likelihood convergence testing.

The second tolerance and maximum-iteration controls are passed to constrained
inner M-step routines when required.

## Hierarchical initialization

`mclust_select` uses an upstream-style hierarchical initialization by default.
The default is `hc_model='VVV'` and `hc_use='SVD'`.  The port also supports
`VARS`, `STD`, `PCS`, `PCR`, and `SPH` transformations.

## BIC and ICL conventions

mclust uses the convention

`BIC = 2*logLik - nParams*log(n)`

so larger values are better.  The ICL implementation applies the classification
entropy correction to BIC.

## Fortran evaluation-order corrections

Fortran does not guarantee C/R-style short-circuit evaluation.  During strict
bounds-checked testing, index-zero hazards in translated sort/label helpers
were rewritten with explicit nested branches.  This changes no statistical
formula; it makes the intended condition portable Fortran.

## Scope notes

Bayesian-prior fitting, background-noise mixtures, GMMHD, MclustDR subset
selection, and several R convenience wrappers are not implemented in v0.1.2.
See `API_MAPPING.md` for the detailed boundary.
