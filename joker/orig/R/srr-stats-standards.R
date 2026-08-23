#' NA_standards
#'
#' @srrstatsNA {G2.4e, G2.5} No inputs are converted from `factor` to another
#' type.
#' @srrstatsNA {G2.6} The S4 structure of the package ensures that unidimensional
#' input is of an appropriate class.
#' @srrstatsNA {G2.7, G2.8, G2.9, G2.10, G2.11, G2.12} The package does not use
#' data.frames, tibbles, or similar structures. The package's S4 methods make
#' sure that the only part of the package that uses a `data.frame` is a slot of
#' the `SmallMetrics` and `LargeMetrics` S4 classes, whose `setValidity()`
#' function imposes a strict structure. Therefore, no conversion between
#' `matrix`/`array` and `data.frame` (or similar structure) that could
#' potentially lose or add information is used in the package.
#' @srrstatsNA {G3.1, G3.1a} The package does not include ambiguous covariance
#' calculations. All variance-covariance matrices are appropriately defined in
#' the relevant literature.
#' @srrstatsNA {G5.0} Tests rely on simulated data so that the true parameter
#' values are known. Using real data sets is applicable but not practicable.
#' @srrstatsNA {G5.4c} Where applicable, stored values may be drawn from
#' published paper outputs when applicable and where code from original
#' implementations is not available
#' @srrstatsNA {G5.11, G5.11a} No data download or relevant procedure is
#' included in the package tests.
#'
#' The following are TODO, not NA
#' @srrstatsNA {G5.4} Correctness tests are included in the
#' package, testing that statistical algorithms produce expected results.
#' some fixed test data sets (potentially
#' through comparisons using binding frameworks such as
#' [RStata](https://github.com/lbraglia/RStata)).*
#' @srrstatsNA {G5.4a} *For new methods, it can be difficult to separate out
#' correctness of the method from the correctness of the implementation, as
#' there may not be reference for comparison. In this case, testing may be
#' implemented against simple, trivial cases or against multiple implementations
#' such as an initial R implementation compared with results from a C/C++
#' implementation.*
#' @srrstatsNA {G5.4b} *For new implementations of existing methods,
#' correctness tests should include tests against previous implementations. Such
#' testing may explicitly call those implementations in testing, preferably from
#' fixed-versions of other software, or use stored outputs from those where that
#' is not possible.*
#' @srrstatsNA {G5.5} *Correctness tests should be run with a fixed random
#' seed*
#' @srrstatsNA {G5.8, G5.8a, G5.8b, G5.8c, G5.8d} Edge condition tests are
#' included in
#' the package, testing that these conditions produce expected behavior such as
#' clear warnings or errors when confronted with zero-length data, data of
#' unsupported types, data with all-`NA` fields, or data outside the scope of
#' an algorithm.
#' @srrstatsNA {G5.9, G5.9a, G5.9b} Noise susceptibility tests are included in
#' the package, verifying that trivial noise, seeds, and initial conditions do
#' not meaningfully change results
#' @noRd
NULL
