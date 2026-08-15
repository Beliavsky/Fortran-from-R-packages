# skewunit-fortran

`skewunit-fortran` is a standalone modern Fortran translation of the
computational code in version 1.1 of the R package **skewunit** by Diego
Gallardo, Emilio Gomez-Deniz, Osvaldo Venegas, and Hector W. Gomez.

The upstream package implements estimation, simulation, and model selection
for skew-unit models on `(0,1)`. This port removes the R, `stats`, and `pracma`
runtime dependencies and provides the same statistical calculations through
an FPM library.

## Implemented computational functionality

The five upstream baseline families are included:

- ArcSin (`asin`)
- U-quadratic (`Uquad`)
- symmetric triangular (`triang`)
- Johnson-SB (`JSB`)
- symmetric beta (`sbeta`)

The library also implements:

- `dskewunit`, `pskewunit`, and `rskewunit`
- vector random generation and vector skew-unit CDF evaluation
- maximum-likelihood estimation corresponding to `estimate.skewunit`
- numerical Hessian and covariance/standard-error calculation
- AIC and BIC
- 30-model selection corresponding to `choose.skewunit`
- the exported `cuberoot` utility

There was no plotting code in the supplied upstream package.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example basic
```

The project uses standard Fortran 2018 and has no external library
dependencies.

## Basic use

```fortran
program demo
   use skewunit
   implicit none
   real(dp) :: x(1000)
   type(skewunit_fit_result) :: fit

   call seed_skewunit_rng(12345)
   call rskewunit_vec(size(x),x,lambda=-0.4_dp,delta=1.2_dp, &
      family1=family_triang,family2=family_jsb)

   call estimate_skewunit(x,family_triang,family_jsb,fit)
   print *, fit%coefficients(1:fit%npar)
   print *, fit%loglik, fit%aic, fit%bic
end program demo
```

Family constants are

```fortran
family_none
family_asin
family_uquad
family_triang
family_jsb
family_sbeta
```

`family_id("asin")` and `family_name(family_asin)` are also provided.

## Baseline distribution API

Scalar density/CDF functions are elemental where practical, so they may also
be applied directly to conformable arrays.

```text
dasin,    pasin,    rasin
dtriang,  ptriang,  rtriang
duquad,   puquad,   ruquad
djsb,     pjsb,     rjsb
dsbeta,   psbeta,   rsbeta
```

Vector RNG subroutines use the `_vec` suffix, for example
`rjsb_vec(n,x,delta)` and `rsbeta_vec(n,x,delta)`.

`duquad`, `puquad`, and `ruquad` retain the upstream optional interval
parameters `a` and `b`; skew-unit models use their default interval `[0,1]`.

## Skew-unit parameter convention

The density is the upstream skewing construction

```text
h(x) = 2 f(x) G(lambda*(x - 0.5) + 0.5),   0 < x < 1.
```

`lambda` is constrained to `[-1,1]`.

Johnson-SB and symmetric-beta families have a positive shape parameter. The
upstream convention is retained:

- if only one of `family1` and `family2` needs a shape, it is `delta`;
- if both need shapes, `delta` belongs to `family1` and `delta2` to `family2`.

Setting `family2=family_none` gives the unskewed baseline model.

## Estimation result

`type(skewunit_fit_result)` contains:

```text
family1, family2
npar
convergence
iterations
coefficients(3)
std_error(3)
std_error_available
loglik
aic
bic
```

The coefficient order follows the R package:

- skew model with no shape: `lambda`
- skew model with one shape: `lambda, delta`
- skew model with two shapes: `lambda, delta1, delta2`
- unskewed JSB/sbeta: `delta`

Use `coefficient_name(fit,i)` for labels.

## Model selection

```fortran
use skewunit
real(dp) :: x(100)
type(skewunit_choice_result) :: choice

call choose_skewunit(x,choice,criteria='BIC')
```

`choice%summary` contains all 25 skew models followed by the five baseline
models, sorted by the selected criterion. `choice%best_fit` is a refit of the
winning model.

## Numerical implementation

The R dependencies have standalone replacements:

- regularized incomplete beta by continued fractions
- standard normal CDF through `erfc`
- beta RNG through Gamma variates
- Marsaglia-Tsang Gamma RNG
- Brent minimization for one-parameter transformed fits
- Nelder-Mead for two- and three-parameter transformed fits
- central finite-difference Hessians
- pivoted Gauss-Jordan inversion for the at-most-3-by-3 Hessian
- specialized adaptive Gauss-Kronrod integration for skew-unit CDFs

For CDF integration, the port uses `x = sin(theta)^2`. This removes the
ArcSin endpoint singularity and substantially improves numerical accuracy for
symmetric-beta shapes below one.

## Tests

The test suite checks:

- independently computed baseline density/CDF values
- independently computed skew density/CDF values
- all 25 combinations of baseline/skewing families
- vector CDF evaluation
- simulation moments/direction
- MLE against an independent SciPy optimization
- Hessian standard errors
- AIC model-selection sorting/refitting

See `PORTING_NOTES.md` for deliberate compatibility choices and differences.

## License

The supplied R package declares `GPL (>= 2)`. This translation is therefore
distributed under **GPL-2.0-or-later**. See `COPYING.GPL-2`, `LICENSES.md`,
`NOTICE.md`, and the retained material in `upstream/`.
