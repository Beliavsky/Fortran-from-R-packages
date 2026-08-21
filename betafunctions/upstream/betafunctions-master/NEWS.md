# betafunctions v. 1.9.0

- Fixed a bug with fitting the two-parameter beta-distribution model in the `HB.CA()`, `HB.CA.MC()`, and `HB.ROC()` functions. 

- Significantly improved computation time for consistency in the `HB.CA()` and `HB.CA.MC()` functions.

- To reduce the verbosity of the output, the `LL.`and `HB.` functions no longer include the positive and negative likelihood ratios as part of their output.

- The `AMS()` and `BMS()` functions have had some arguments removed. Now only allows for specifying mean, variance (or standard deviation), and location parameters.

- The `LABMSU()` and `UABMSL()` functions have had arguments removed (`skewness` and `kurtosis`) due to arguments being removed from `AMS()` and `BMS()`.

- The `Beta.tp.fit()` function is significantly simplified and some arguments (`alpha` and`beta`) are removed.
 
- The `d/p/q/rBetaMS()` functions now includes arguments for specifying lower- and upper-bound location parameters. This allows the mean to be specified outside the [0, 1] interval.

---

# betafunctions v. 1.8.1

- Changes to model fit argument in the `LL.` and `HB.` functions.

  - The classification functions `LL.CA()`, `LL.CA.MC()`, `HB.CA`, and `HB.CA.MC()` functions now allows for specifying not to perform model-fit testing by passing `NULL` to the `modelfit` argument.
  
  - Fixed a bug that affected the binning of of observed values. Due to some strange behavior in R that resulted from scaling and re-scaling values, observations would sometimes be assigned to the wrong bins. This should now be fixed.

- Added arguments to the `mdlfit.gfx()` function that allows for manually controlling a number of plot aesthetics.

---

# betafunctions v. 1.8.0

- Added an implementation of the Hanson-Brennan approach to classification accuracy and consistency.

  - The `HB.CA()` and `HB.MC.CA()` functions for binary and multiple classifications, respectively.
  
  - The `HB.ROC()` function for ROC analysis.
  
- Added d/p/r functions for Lord's two-term approximation to the compound binomial distribution (`xcBinom()`).

- Added d/p/r functions for Beta compound-Binomial distributions where the compound Binomial is Lord's two-term approximation (`xBetacBinom()`).

- Added the `mdlfit.gfx()` as a visual aid to gauge the severity of model misfit for the LL and HB approaches.

---

# betafunctions v. 1.7.0

- Added d/p/r functions for Beta-Binomial distributions (`d/p/rBetaBinom()`, no `q` function as of yet).

- Added the `betabinomialmoments()` function for calculating raw-, central-, and standardized moments of Beta-Binomial distributions.

- Added the `R.ETL()` function for calculating the model-implied reliability of test-scores given the "Effective Test Length" of Livingston and Lewis (1995), and the mean, variance, and minimum and maximum possible scores of the observed-score distribution.

- Added the `mdo()` function which allows for estimating McDonald's Omega reliability coefficient. Requires that there are no negative item-covariances.

- Added the `MC.out.tabular()` function which can be used to organize the accuracy and consistency output from the `LL.CA.MC()` function in tabular format. The `MC.out.tabular()` function takes the output of the `LL.CA.MC()` function as input.

- Changes to the model-fit test procedure for the `LL.CA()`and `LL.CA.MC()` functions. The initial number of bins are now set to 100 and the minimum bin-size is set to 10. The bins are now only grouped into the expected number of observations rather than both expected and observed. The behaviour of this model-fit test-procedure is still under scrutiny and should as of yet be considered an experimental approximation procedure.

---

# betafunctions v. 1.6.1

- Added `betamode()` and `betamedian()` functions that allow for calculating the mode and median (respectively) of two- and four-parameter Beta distributions.

- Added `modelfit` arguments to the `LL.CA()` and `LL.CA.MC()` functions that allow for controlling the maximum number of- and minimum size of the bins that are to be used for the chi-square test of model fit. Tuning the maximum number and minimum size of the bins should facilitate model-fit testing if the default settings does not result in sufficient degrees of freedom to perform the significance-test.

- Fixed bug in the `LL.CA.MC()` function regarding the assembly of the consistency matrix when more than two categories were included which resulted in some entries having slightly too large values. This change results in slight changes to the estimated consistency-indices (now more precise).

- The `LL.CA()` and `LL.CA.MC()` functions no longer issue warnings if parameter estimates are out of bounds if the `true.model` argument is specified to `"2P"` and the `failsafe` argument is set to `TRUE` (the default).

---

# betafunctions v. 1.6.0

- Major updates to the classification accuracy and consistency functionalities:

  - Added new `LL.ROC()` functionality, providing the option to give the ROC-curve a classic, "staircase" look. This option is made the default.

  - Added the `LL.CA.MC()` function extending the Livingston and Lewis approach to using multiple cut-points. The output of this function can grow quite verbose when operating with several cut-points, as diagnostic performance and consistency indices are estimated and reported for each group separately.
  
  - Removed the `error.model` argument of the `LL.CA()` function. In essence, this means that one can no longer specify a beta error model.
  
  - Added model-fit functionality to the `LL.CA()` function. The model-fit is examined by comparing observed and model-expected frequencies by binning frequencies into bins of at least 10 expected and observed frequencies. According to Lord, this model fit procedure has N-bins - 4 degrees of freedom, meaning that model-fit cannot be examined if a grouping of the observed-score distribution reduces to fewer than 5 bins.

---

# betafunctions v. 1.5.0

- Added the `gchoose()` function generalizing the base-R `choose()` function to work with non-integers and positive integers by calculating the factorials of the Binomial coefficient by drawing on the Gamma distribution.

- Added new set of d/p/q/r functions for a new distribution: The "Gamma-Binomial" distribution. This distribution extends the Binomial distribution to all positive real numbers. Not defined for negative integers.
  
  - `dGammaBinom()`: Probability Density Distribution for the Gamma-Binomial distribution.
  
  - `pGammaBinom()`: Cumulative Probability Density Distribution for the Gamma-Binomial distribution.
  
  - `qGammaBinom()`: Quantile function for the Gamma-Binomial distribution. Calls the `pGammaBinom()` function and utilizes a bisecting search-algorithm to find the number of "successful trials" corresponding to the quantile in question. This algorithm is rather slow so the function might take longer to find the appropriate quantile than what one might be used to.
  
  - `rGammaBinom()`: Random number generation for the Gamma-Binomial distribution. Calls the `qGammaBinom()` function. Since the `qGammaBinom()` function searches for the appropriate values using a rather inefficient search-algorithm (bisection), this random-number generation is somewhat slow. 

- Added the `binomialmoments()` function, which allows for calculating the raw, central, and standardized moments of Binomial distributions (for which the Beta distribution is the conjugate prior).

- Changes to the `LL.ROC()` function. 

  - Changed the internal behavior of the function to estimate the true-score distribution only once. This should greatly improve the time required to produce the plots.

  - Added the `locate` argument where it is possible to ask the function to locate the operational cut-point at which the values of sensitivity or NPV are greater than or equal to some value, or specificity or PPV are lesser than or equal to some value.
  
  - Added the `maxAcc` argument to locate the cut-point at which the Accuracy statistic is maximized.
  
  - The raw-output print-out now contains the cut-point specific Accuracy, PPV, and NPV statistics as well.
  
- Changes to the `afac` and `dfac` functions.

  - Due to the problems associated with calculating descending and ascending factorials for low-valued integers by means of the gamma function (particularly for descending factorials), a direct-arithmetic solution is implemented and set as the default method. In order to use the gamma function rather than direct arithmetic, specify any value other than `"product"` as part of the `method` argument.
  
- Changes to the `tsm` function.

  - The `tsm` function now calls the `dfac` function with the direct-arithmetic method for calculating descending factorials as default. In order to use the gamma function rather than direct arithmetic, specify any value other than "product" as part of the `method` argument.

- Correction to the `LL.CA()` function. The Binomial error distribution evaluated up to and including the cut-point. The intended behaviour was to evaluate up to but NOT including the cut-point. Prior to this correction, it is expected that the `LL.CA()` function will have underestimated accuracy somewhat.

  - This new behavior applies to the `dBeta.pBinom()` function as well, which is contrary to the default behavior of base-R's `pbinom()` function.

---

# betafunctions v. 1.4.4

- Fixed a bug with the true-score distribution fitting procedure in the `Beta.tp.fit()` that could occur for low integer values.

- Minor changes to the documentation for various functions.

---

# betafunctions v. 1.4.3

- Added possibility of calculating descending (falling) and ascending (rising) factorials by means of the `dfac()` and `afac()`functions.

- Added `tsm()` argument for calculating raw moments of the true-score distribution under the Livingston and Lewis approach.

- Added `confmat()` function for organizing supplied values of true and false positives and negatives into a confusion matrix.

---

# betafunctions v. 1.4.2

- Changes to the `LL.ROC()` function. Now allows for specifying the lower- and upper-bound parameters of the true-score distribution should "2P" be specified.

- A Shiny application providing a GUI for the Livingston and Lewis approach functionality of the package has been developed and can be found at https://hthaa.shinyapps.io/shinybeta/

- The "3P" functionality of `Beta.tp.fit()`, `LL.CA()` and `LL.ROC` has been removed, as it did not perform satisfactorily. 

---

# betafunctions v. 1.4.1

- The `true.model` argument of the `Beta.tp.fit()` function now includes a `"3P"` option, allowing for the specification of one location-parameter (l or u) and estimating the remaining location-parameter and the shape-parameters (alpha and beta) so as to make the resulting distribution have the same skewness and kurtosis as the estimated true-score distribution.

- The `AMS()` and `BMS()` functions now issue warnings if there was not enough information supplied to calculate the target parameter.

- The `ETL()` function: For the sake of argument consistency, the `l` and `u` arguments are not renamed `min` and `max`, respectively.

---

# betafunctions v. 1.4.0

- Added a number of new functions for working with Beta distributions.

	- The new `LABMSU()` function allows for finding the lower-bound parameter for the four-parameter Beta distribution by supplying the shape-parameters and moments of the resulting distribution, and (optionally) the upper-bound location parameter.
	
	- The new `UABMSU()` function allows for finding the upper-bound parameter for the four-parameter Beta distribution by supplying the shape-parameters and moments of the resulting distribution, and (optionally) the lower-bound location parameter.

- Added additional functionality to some existing functions, allowing for specifying the lower- and upper-bounds, and/or moments of the resulting distribution.

	- The `AMS()` function now includes `l` and `u` arguments, finding the Alpha shape parameter necessary to produce a Beta distribution with target moments and specified lower and upper bounds of the resulting distribution.

	- The `BMS()` function now includes `l` and `u` arguments, finding the Beta shape parameter necessary to produce a Beta distribution with target moments and specified lower and upper bounds of the resulting distribution.

	- The `Beta.2p.fit()` function, essentially a wrapper-function for `AMS()` and `BMS()`, now includes `l` and `u` arguments for specifying the bounds for the resulting distribution.

- Functions relating to estimating classification accuracy and consistency:
	
	- `Beta.tp.fit()`, `LL.CA()`, and `LL.ROC()`: 
	
		- Added `true.model` argument allowing for greater control over the true-score estimation procedure in.

		- The `true.model` allows for specifying whether to fit four- or two-parameter Beta distribution to the estimated moments of the true-score distribution.

		- Currently, the options for the argument are "4P" and "2P". Further options might be added in the future.

	- The `Beta.tp.fit()` function: 

		- Now allows for specifying the lower- and/or upper bounds of the two-parameter solutions.

		- Now allows for specifying that the fitted distribution should estimate two parameters rather than the default of four. Estimates then the parameters necessary to produce a distribution with the same mean and variance as the estimated true-score distribution, given specified moments, shape- and location parameters (default is `l = 0` and `u = 1`).

		- Added the option of specifying the true-score distribution moments (sans a functional form) as output rather than the estimated parameters of the Beta true-score distribution in the new `output` argument. The default is to retrieve the estimated parameters of the Beta true-score distribution.

		- When `parameters` is specified in the `output` argument, the effective test length is included as part of the output. As such, the default output of the `Beta.tp.fit()` function is now complete for the purposes of being used as input for the `LL.CA()` function.

	- The `LL.CA()` function:

		- The `override` and `failsafe` arguments:

			- The `override` argument is rendered inert, to be removed in a later version.

			- The `failsafe` argument is introduced to replace the `override` argument. It is essentially an inversion of the `override` argument.

			- The above change is to streamline and make consistent arguments structures across functions.


- General code cleaning:

	- Changed some argument names so as to have the names be consistent across functions. Argument locations in the function calls remain unchanged.

		- Function arguments previously named `lt` now named `lower.tail`.

		- Function arguments previously named `a` or `b` now named `alpha` and `beta` (respectively).

		- Function arguments previously named `var` now named `variance`.


- Added the possibility of specifying a Beta error model for the `LL.ROC()` function by way of the `true.model` argument.

- Additional correction of typographical errors in documentation.

---

# betafunctions v. 1.3.0

- Fixed some typographical errors in documentation.

- Minor modification to `LL.CA()` where the function now won't automatically terminate if true-score distribution fitting procedure produces `NA` or `NaN` estimates.

	- True-score distribution fitting now performed using the new `Beta.tp.fit()` function.

	- The `grainsize` argument is now fully removed from `LL.CA()`.

	- Added additional checks and diagnostics in `LL.CA()`, issuing warnings for aberrant events.

	- If a list of parameter values are supplied in place of a score-vector, it is now necessary to supply an `etl` (effective test length) parameter as well. See documentation for the `ETL()` function for more information.

- Fixed an error in the fitting procedure of `Beta.4P.fit()` which occurred during positive skew.

- Added the `Beta.tp.fit()` function for estimating four-parameter Beta distribution parameters for an underlying true-score distribution, assuming that the observations are generated from a Beta-Binomial model.

---

# betafunctions v. 1.2.1

- Fixed the `LL.ROC()` function which stopped working after the previous update.

	- Included possibility of specifying the number of points for which the ROC curve is drawn and the AUC statistic is calculated by way of the new `grainsize` argument.

---

# betafunctions v. 1.2.0

- Added functionality for calculating consistency statistics by way of the new `ccStats()` function.

- Substantial developments concerning primarily the `LL.CA()` function, adding functionality and improving performance.

	- Added calculation of consistency statistics by calling the new `ccStats()` function.

	- Added `output` argument indicating which statistics to calculate. Default is to compute both accuracy and consistency.

	- Calculation of distribution-based output now utilizes the `integrate()` function.

		- Subsequently rendered the `grainsize` argument inert, as it was used for the previous method.

	- Added possibility of supplying list of custom parameter values for the four-parameter Beta true-score distribution, forgoing the need to estimate one.

	- Added `override` argument providing the possibility of overriding the default fail-safe reverting to a two-parameter Beta true-score distribution, should the fitting procedure produce impermissible parameter estimates.

---

# betafunctions v. 1.1.0

- Fixed small typographical errors in documentation.

- Added the `cba()` function for computing Cronbach's alpha.
