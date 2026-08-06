# RPESE 1.0.0
* Initial stable release of package.

# RPESE 1.0.1
* Fix bug for SemiSD.
* Add functionality of package for frequency decimation and truncation.

# RPESE 1.0.2
* Add the option to return the coefficients from the Exponential/Gamma fit.

# RPESE 1.0.3
* Fix error for the StdDev.SE function.
* Fix error for the SortinoRatio.SE function.
* New option to return the correlations of the returns and the IF transformed returns.

# RPESE 1.0.4
* Fix error for prewhitening.

# RPESE 1.0.5
* Another error for prewhitening. Now the residuals are passed on for the GLM fit.

# RPESE 1.0.6
* Typo for prewhitening (variable name).

# RPESE 1.1.0
* Adapt RPESE to changes to RPEIF 1.1.0

# RPESE 1.2.0
* Add the robust m-estimator of location.
* Update vignette for new robust estimator and fix old typos.

# RPESE 1.2.1
* Update the name for robust mean function.
* Change default SE methods for performance measures (Adapt -> PW)

# RPESE 1.2.2
* Update RPEIF function name for mean estimator.
* Remove edhec data.

# RPESE 1.2.3
* Fix code errors in the vignette.

# RPESE 1.2.4
* Update Description file reference to technical article.

# RPESE 1.2.5
* Fix error in point estimate computation for Value-at-Risk.

# RPESE 1.2.6
* Update author field in DESCRIPTION file.

# RPESE 1.2.7
* Removed redundant Maintainer field from DESCRIPTION (already specified in Authors@R).
* Added Encoding: UTF-8 to DESCRIPTION.
* Updated all examples to conditionally check for PerformanceAnalytics availability 
  using requireNamespace(), ensuring compatibility with noSuggests CRAN checks.
* Package now compatible with RPEIF 1.2.6 and RPEGLMEN 1.1.5.