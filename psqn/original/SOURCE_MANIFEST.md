# Upstream source manifest

Translation input: psqn 0.3.2, CRAN/GitHub source package supplied by the user.

Computational sources inspected and translated:

- inst/include/constant.h
- inst/include/intrapolate.h
- inst/include/lp.h
- inst/include/psqn-bfgs.h
- inst/include/psqn-misc.h
- inst/include/psqn.h
- inst/include/richardson-extrapolation.h
- src/r-api.cpp (used to identify exported computational surfaces)

Interface/support sources inspected but not translated as numerical code:

- inst/include/psqn-Rcpp-wrapper.h
- inst/include/psqn-reporter.h
- src/RcppExports.cpp
- R/RcppExports.R
- R/catch-routine-registration.R

Presentation/vignette/plotting material is outside the standalone numerical port.
