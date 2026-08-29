# Porting notes

## Scope

This port targets the computational/statistical content of `tweedie` 3.1.0.
Plotting, console/progress output, deprecated dot-name aliases, R formula and
`model.frame` processing, and S3/model-object presentation are not reproduced.

The numeric core of profile likelihood fitting is retained through a
Fortran design-matrix IRLS routine, phi likelihood optimization, and a grid
profile routine. The R wrapper's `splinefun` interpolation, plotting, and
confidence-interval presentation are generic R post-processing and are left
to callers.

## Reused upstream Fortran

The upstream package already contains the Dunn-Smyth Fourier-inversion engine.
The following numerical sources are copied unchanged from `tweedie/src`:

- `00tweedie_params.f90`
- `gaussianData.f90`
- `Calcs_Solvers.f90`
- `Calcs_Real.f90`
- `Calcs_Imag.f90`
- `Integrands.f90`
- `GaussQuadrature.f90`
- `findAccelStart.f90`
- `accelerate.f90`
- `Calcs_K.f90`
- `TweedieIntHelpers.f90`
- `TweedieIntegration.f90`
- `twcomputation_loop.f90`

`R_interfaces.f90` and `rprintf_mod.f90` are standalone replacements for the
R runtime printing glue. The C registration/wrapper files are unnecessary in
a native Fortran library and are not built.

## R-level algorithms translated

### Density

`dtweedie` preserves the package's parameterization and dispatch rules:
closed forms for p=1,2,3; exact mass at zero for 1<p<2; Dunn-Smyth series in
series regions; and two-dimensional interpolation over the stored Fourier
inversion grids in the same p/xix regions as the R implementation.

All ten stored grids were copied from `dtweedie_interpolation.R`. There are
416 coefficients per grid (4,160 total). A programmatic audit against the
attached R source found exact floating-point equality for every value.

A subtle upstream distinction is retained: public `dtweedie(..., power=1)`
uses the scaled-lattice `dpois(y/phi, mu/phi)` convention, whereas the
lower-level `dtweedie_series(..., power=1)` helper divides this value by
`phi`, exactly as its R implementation does.

### Series calculations

Both density expansions from Dunn and Smyth (2005) are translated:

- the nonalternating W-series for 1<p<2;
- the alternating V-series for p>2.

The internal `logw`, `jw`, `logv`, and `kv` calculations are exposed because
the upstream likelihood derivative uses them. Summation is stabilized with a
largest-log-term shift.

`ptweedie_series` implements the 1<p<2 compound-Poisson/gamma CDF series.

### Fourier inversion

`dtweedie_inversion` retains the three scaling representations described by
Dunn and Smyth (2008) and exposes the upstream exit status, relative-error
estimate, and number of integration regions. If no method is supplied, the
same notional multiplier criterion is used to select among the three
representations.

`ptweedie_inversion` calls the upstream CDF inversion directly. The optional
`igexact=.false.` path deliberately sends p=3 through Fourier inversion, which
is useful for parity testing against the exact inverse-Gaussian result.

### Closed forms and DPQR

- p=1: scaled Poisson lattice
- p=2: gamma with shape 1/phi and scale phi*mu
- p=3: inverse Gaussian with variance `phi*mu**3`
- 1<p<2 RNG: direct Poisson-gamma compounding
- p>2 RNG: inverse-CDF generation, matching the R implementation

The inverse-Gaussian density/CDF are included locally because they are needed
for the p=3 special case and are not supplied by `r_mod`.

### Likelihood/profile computations

`tweedie_loglik`, `tweedie_aic`, the phi score, and the internal derivative
of log-density with respect to phi are ported. The p>2 derivative retains the
upstream finite-difference fallback in difficult regions and the phi-scaling
stabilization in `dtweedie_dldphi`.

`tweedie_phi_mle` performs positive one-dimensional likelihood optimization.
`tweedie_glm_fit` implements IRLS for a caller-supplied design matrix and the
power link used by the upstream `statmod::tweedie` family. The common
`link_power=0` log-link path is fully supported. `tweedie_profile_grid` fits
the GLM and profiles phi across caller-supplied Tweedie powers without
requiring R's formula/model-object machinery.

As in the upstream `tweedie_profile`, observation weights affect the GLM fit,
while the profile density likelihood itself is evaluated on the observations
without multiplying the log densities by those weights. The standalone
`tweedie_loglik` additionally accepts optional weights for AIC-style uses.

## r_mod reuse

The user-supplied MIT-licensed `r_mod.f90` is used for R-compatible Poisson,
gamma and normal helpers, RNGs, digamma, and linear solves. No duplicate
versions of those helpers were added. The build copy is a formatting-only
line-wrapped version of the exact supplied file; the original is retained as
`upstream/r_mod-original.f90`.

## Validation

Validation was performed with GNU Fortran 14.2.0 using Fortran 2018,
`-Werror=implicit-interface`, and `-fcheck=all`, linked against BLAS/LAPACK.
The clean validation suite covers:

- p=1/p=2/p=3 special cases;
- inverse-Gaussian closed-form references;
- independent compound-Poisson/gamma references for p=1.5, mu=1, phi=4;
- agreement of series and Fourier inversion on both sides of p=2;
- forced p=3 Fourier inversion against the exact inverse Gaussian;
- all ten stored interpolation-grid regions against direct inversion;
- historical inversion regression inputs copied from the upstream tests;
- analytic phi score versus finite differences;
- phi maximum-likelihood estimation;
- intercept-only Tweedie GLM IRLS and profile-grid fitting.

FPM is not installed in the validation environment, so the project was built
with an equivalent explicit dependency-order `gfortran` command rather than
literally running `fpm test`.
