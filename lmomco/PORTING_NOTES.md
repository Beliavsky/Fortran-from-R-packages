# Porting notes

## Corrected source archive

An earlier `lmomco-fortran-v0.1.0.zip` accidentally contained only documentation. This rebuilt release contains the actual Fortran sources, tests, example, FPM manifest, and upstream source archive.

## API design

R's list-based `para` objects are represented by `type(lmomco_params)`. Use `make_params(family, values)` to construct parameters and `lmomco_pdf`, `lmomco_cdf`, and `lmomco_quantile` for generic dispatch.

For quantile-defined distributions such as GLD, Govindarajulu, LMRQ, PDQ3/4, and Wakeby, the CDF is obtained by bracketed inversion of the translated quantile function; the PDF is computed from the derivative of the quantile function numerically when no simpler stable closed form is used.

Rice, Eta-Mu, and finite Kappa-Mu use their translated density formulas and numerical integration/inversion. The Kappa-Mu `kappa = Inf` Dirac-mass limit is not implemented in v0.1.0.

## L-moments

Sample L-moments use unbiased probability-weighted moments. Theoretical L-moments use numerical integration of the quantile function against shifted Legendre polynomials.

Direct Hosking-style fits are currently implemented for Normal, Exponential, Gumbel, and GEV. More family-specific `par*` estimators remain worthwhile follow-up targets.

## Numerical notes

The implementation uses self-contained regularized incomplete-gamma and incomplete-beta routines and bracketed quantile inversion. This avoids requiring R, MASS, goftest, or Lmoments at runtime for the currently translated core.
