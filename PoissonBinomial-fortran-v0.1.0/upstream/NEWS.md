# PoissonBinomial 1.2.8

* Addresses future `Rcpp` changes that would cause compiler warnings
* Minor documentation updates


# PoissonBinomial 1.2.7

* Updated `Makevars` file, because it contained hard-coded paths which could
  cause problems on some exotic configurations, e.g. MacPorts
* Minor performance improvement for DFT-CF procedure
* Documentation updates


# PoissonBinomial 1.2.6

* Addresses `[-Wformat-security]` compiler warning on newer clang configuration.
  See https://github.com/RcppCore/Rcpp/issues/1287
* change of maintainer email to personal address
* GitHub repo moved to personal account


# PoissonBinomial 1.2.5

* Minor performance improvements for exact methods of `dgpbinom` and `pgpbinom`.
* Added GitHub.


# PoissonBinomial 1.2.4

* Performance improvements and minor bug fix for quantile functions.


# PoissonBinomial 1.2.3

* Further optimizations of determining the number of splits for "DivideFFT"
  procedures.


# PoissonBinomial 1.2.2

* Performance improvements for "Convolve" (and subsequently "DivideFFT")
  procedures.
* GCD optimizations for generalized Poisson distributions have been moved to
  the respective C++ functions, so that packages that import them may benefit
  from them as well.
* Removed dependence on `BH` package, as it was only needed for one constant in
  the C++ code, which has now been defined manually.
* Adjustment to vignettes. Benchmarks (for performance comparisons) are now
  done with 51 runs (instead of 100) to accelerate building of the package.
* Fixed a minor bug with `qgpbinom`.


# PoissonBinomial 1.2.1

* Fixed a minor code issue that prevented compilation on Solaris systems.


# PoissonBinomial 1.2.0

* Performance improvements for exact methods of `dgpbinom` and `pgpbinom`.
* Input variable of the Rcpp implementations of all methods are made `const`
  to prevent inadvertent changes to them. Any package that imports headers must
  be updated. The 'Imports' field of the DESCRIPTION file should include a
  version requirement, i.e. PoissonBinomial (>= 1.2.0).
* Added new random generation methods for `rpbinom` and `rgpbinom` (see
  function documentation). They are much faster than the old quantile-based
  inversion method, which has been removed.


# PoissonBinomial 1.1.3

* Improved numerical accuracy of normal approximations of `dpbinom` and 
  `dgpbinom`.
  

# PoissonBinomial 1.1.2

* Bug fixes and performance improvements for `qpbinom` and `qgpbinom` that also
  affect `rpbinom` and `rgpbinom`. Quantiles were off be one; all code that uses
  the quantile functions should be reviewed!
* When requesting cumulative probabilities, the respective C++ implementations
  are now capable of computing these values for `lower.tail = FALSE` on their
  own, which improves accuracy.


# PoissonBinomial 1.1.1

* Bug fixes in `ppbinom` and `pgpbinom` that caused incorrect calculation of
  logarithms and cumulative upper-tail probabilities.


# PoissonBinomial 1.1.0

* Added exact and approximate algorithms for the generalized Poisson binomial
  distribution described in Zhang, Hong & Balakrishnan (2018). The
  non-generalized distribution is now referred to as the 'ordinary' Poisson
  binomial distribution.
* Restructured vignettes. Added tables of content and fixed smaller issues.
* Minor bug fixes for `dbinom`, `ppbinom` and `qpbinom` functions.


# PoissonBinomial 1.0.2-1

* Fixes and improvements of the vignettes; no code changes.


# PoissonBinomial 1.0.2

* Improvements of C++ helper function `norm_dpb` to achieve better
  normalization.
* Bug fix of DFT-CF method ("Characteristic") so that negative probabilities
  are no longer possible.
* Reworked vignette structure.
* Added author acknowledgements to the Makevars.win file (original author was
  Geoff99 (https://github.com/Geoff99)).
  

# PoissonBinomial 1.0.1

* Fixed a bug in the C++ helper function "norm_dpb" that could cause infinite
  loops (the function is invisible to the user, since it is only used in the
  C++ domain).
  

# PoissonBinomial 1.0.0

* Initial release.
* Added a `NEWS.md` file to track changes to the package.
