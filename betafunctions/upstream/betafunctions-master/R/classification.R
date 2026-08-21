#' Livingston and Lewis' "Effective Test Length".
#'
#' @description  According to Livingston and Lewis (1995), "The effective test length corresponding to a test score is the number of discrete, dichotomously scored, locally independent, equally difficult items required to produce a total score of the same reliability."
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param min The lower-bound (minimum possible value) of the observed-score distribution. Default is 0 (assuming observed scores represent proportions).
#' @param max The upper-bound (maximum possible value) of the observed-score distribution. Default is 1 (assuming observed scores represent proportions).
#' @param reliability The reliability of the observed scores (proportion of observed-score distribution variance shared with true-score distribution).
#' @return An estimate of the effective length of a test, given the stability of the observations it produces.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the reliability of this test was estimated to 0.7. To estimate and
#' # retrieve the effective test length using ETL():
#' ETL(mean = mean(testdata), variance = var(testdata), min = 0, max = 100,
#' reliability = .7)
#' @export
ETL <- function(mean, variance, min = 0, max = 1, reliability) {
  ((mean - min) * (max - mean) - (reliability * variance)) / (variance * (1 - reliability))
}

#' Model Implied Reliability from Livingston and Lewis' "Effective Test Length".
#'
#' @description  Calculate model-implied reliability given mean, variance, the minimum and maximum possible scores, and the effective test length.
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param min The lower-bound (minimum possible value) of the observed-score distribution. Default is 0 (assuming observed scores represent proportions).
#' @param max The upper-bound (maximum possible value) of the observed-score distribution. Default is 1 (assuming observed scores represent proportions).
#' @param ETL The effective test length as defined by Livingston and Lewis (1995).
#' @return An estimate of the reliability of a test, given the effective test length, mean, variance, and minimum and maximum possible scores of the observed-score distribution..
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # From the data-generating script above, the effective test length is 100.
#' # To estimate and retrieve the model-implied reliability using R.ETL():
#' R.ETL(mean = mean(testdata), variance = var(testdata), min = 0, max = 100,
#' ETL = 100)
#' @export
R.ETL <- function(mean, variance, min = 0, max = 1, ETL) {
  (ETL * variance + min * (max - mean) + mean * (mean - max)) / ((ETL - 1) * variance)
}

#' An Implementation of the Livingston and Lewis (1995) Approach to Estimate Classification Consistency and Accuracy based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{etl} for the effective test length (see documentation for the \code{ETL} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param min The minimum value possible to attain on the test. Default is 0.
#' @param max The maximum value possible to attain on the test. Default is 1 (assumes that the values contained in \code{x} represents proportions of maximum credit).
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves (see documentation for the \code{LL.ROC()} function).
#' @param output Character vector indicating which types of statistics (i.e, accuracy and/or consistency) are to be computed and included in the output. Permissible values are \code{"accuracy"} and \code{"consistency"}.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit. The argument takes either a vector of two values, or \code{NULL}. If set to \code{NULL}, the model-fit test is not executed. If a vector of values is supplied, the first value is to represent the initial number of bins the distribution of scores is to be divided in to. This value is set to a default of 100. If this default results in too few bins to conduct the chi-square test, this value can be made larger. The second value represents the minimum expected number of observations that the bins should consist of. In accordance with standard recommendations for chi-square tests, the default value is set to 10.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false pass/fail categorizations for a test, diagnostic performance statistics, and / or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note It should be noted that this implementation differs from the original articulation of Livingston and Lewis (1995) in some respects. First, the procedure includes a number of diagnostic performance (accuracy) indices which the original procedure enables but that were not included. Second, the way consistency is calculated differs substantially from the original articulation of the procedure, which made use of a split-half approach. Rather, this implementation uses the approach to estimating classification consistency outlined by Hanson (1991).
#' @note A shiny application providing a GUI for this method is available at https://hthaa.shinyapps.io/shinybeta/ .
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To estimate and
#' # retrieve the estimated parameters, confusion matrix, consistency and
#' # accuracy statistics using LL.CA():
#' LL.CA(x = testdata, reliability = .7, cut = 50, min = 0, max = 100)
#'
#' # Suppose the true-score parameter estimation procedure arrived at
#' # impermissible parameter estimates (i.e., l < 0, u > 1, alpha < 0, or
#' # beta < 0). For example:
#' set.seed(9)
#' testdata <- rbinom(100, 25, rBeta.4P(100, 0.25, 1, 5, 3))
#' Beta.tp.fit(testdata, 0, 25, 25, failsafe = TRUE)
#'
#' # Suppose further that you have good grounds for assuming that the lower-
#' # bound parameter is equal to 0.25 (e.g., the test consists of multiple-
#' # choice questions with four response options, leading to a 25% probability
#' # of guessing the correct answer per question), and good reason to believe
#' # that the upper-bound parameter is equal to 1 (i.e., there is no reason to
#' # believe that there are no members of the population who will attain a
#' # perfect score across all possible test-forms.) To set these lower and
#' # upper bounds for the fitting procedure in the LL.CA() function, set
#' # the argument true.model = "2p", and specify the location parameters
#' # l = 0.25 and u = 1:
#' LL.CA(testdata, 0.6287713, 12, 0, 25, true.model = "2p", l = 0.25, u = 1)
#'
#' # Alternatively to supplying scores to which a true-score distribution is
#' # to be fit, a list with true-score distribution parameter values can be
#' # supplied manually along with the effective test length (see documentation
#' # for the ETL() function), foregoing the need for actual data. The list
#' # entries must be named. "l" is the lower-bound and "u" the upper-bound
#' # location parameters of the true-score distribution, "alpha" and "beta" for
#' # the shape parameters, and "etl" for the effective test-length..
#' trueparams <- list("l" = 0.25, "u" = 0.75, "alpha" = 5, "beta" = 3, "etl" = 50)
#' LL.CA(x = trueparams, cut = 50, min = 0, max = 100)
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
LL.CA <- function(x = NULL, reliability, cut, min = 0, max = 1, true.model = "4P", truecut = NULL, output = c("accuracy", "consistency"), failsafe = TRUE, l = 0, u = 1, modelfit = c("nbins" = 100, "minbin" = 10)) {
  out <- base::list()
  if (!is.list(x)) {
    if ((base::min(x) < min) | (base::max(x) > max)) {
      warning(paste("Observed values not within the specified [", min, ", ", max, "] bounds (observed min = ",
                    base::min(x), ", observed max = ", base::max(x), ").", sep = ""))
    }
    N <- ETL(base::mean(x), stats::var(x), min = min, max = max, reliability = reliability)
    if (startsWith(as.character(true.model), "2")) {
      failsafe <- FALSE
    }
    params <- Beta.tp.fit(x, min = min, max = max, etl = N, true.model = true.model, failsafe = failsafe, l = l, u = u)
    if (params$l < 0 | params$u > 1) {
      warning(paste("Parameter out of bounds: l = ", round(params$l, 4), ", u = ", round(params$u, 4), ", alpha = ", round(params$alpha, 4), ", beta = ", round(params$beta, 4),
                    ". Consider constraining the fitting procedure further (e.g., set the location-parameters).", sep = ""))
    }
    x <- (x - min) / (max - min)
  } else {
    params <- x
    N <- params$etl
  }
  if (base::is.null(truecut)) {
    truecut <- cut
  }
  N <- base::round(N)
  cut <- (cut - min) / (max - min)
  truecut <- (truecut - min) / (max - min)
  params[["etl_rounded"]] <- N
  out[["parameters"]] <- params
  if (!is.list(x) & !is.null(modelfit)) {
    x <- round(x * N)
    tcut <- seq(0, N, N / modelfit[1])
    mdlfit <- matrix(nrow = 2, ncol = length(tcut) - 1)
    for (j in 1:(length(tcut) - 1)) {
      mdlfit[1, j] <- stats::integrate(function(x) {
        if (j == 1) {
          dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(tcut[j + 1] - 1, N, x))
        }
        else {
          if (j != 1 & j != (length(tcut) - 1)) {
            dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(tcut[j + 1] - 1, N, x) - stats::pbinom(tcut[j] - 1, N, x))
          } else {
            dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(tcut[j + 1], N, x) - stats::pbinom(tcut[j] - 1, N, x))
          }
          }
        }, lower = 0, upper = 1)$value
      if (j == ncol(mdlfit)) {
        mdlfit[2, j] <- length(x[x >= tcut[j]])
      } else {
        mdlfit[2, j] <- length(x[x < tcut[j + 1] & x >= tcut[j]])
      }
      }
    mdlfit[1, ] <- (mdlfit[1, ]/sum(mdlfit[1, ])) * length(x)
    for (i in 1:ncol(mdlfit)) {
      if (i < ncol(mdlfit)) {
        if (any(mdlfit[1, i] < ncol(mdlfit))) {
          if (any(mdlfit[1, i] < modelfit[2])) {
            mdlfit[, i + 1] <- mdlfit[, i + 1] + mdlfit[, i]
            mdlfit[, i] <- NA
          }
        }
      }
    }
    mdlfit <- mdlfit[, apply(mdlfit, 2, function(x) { any(!is.na(x)) })]
    if (any(mdlfit[1, ncol(mdlfit)] < modelfit[2])) {
      mdlfit[, ncol(mdlfit) - 1] <- mdlfit[, ncol(mdlfit) - 1] + mdlfit[, ncol(mdlfit)]
      mdlfit <- mdlfit[, -ncol(mdlfit)]
    }
    chisquared <- sum(apply(mdlfit, 2, function(x) {
      (x[2] - x[1])^2/x[1]
    }))
    out[["modelfit"]] <- list()
    out[["modelfit"]][["contingencytable"]] <- mdlfit
    out[["modelfit"]][["chisquared"]] <- chisquared
    if (startsWith(as.character(true.model), "2")) {
      out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
    }
    else {
      if ((startsWith(as.character(true.model), "4") & failsafe == TRUE) & (out[["parameters"]]$l == l & out[["parameters"]]$u == u)) {
        out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
      } else {
        out[["modelfit"]][["df"]] <- ncol(mdlfit) - 4
      }
      }
    out[["modelfit"]][["pvalue"]] <- stats::pchisq(chisquared, ncol(mdlfit) - 4, lower.tail = FALSE)
    }
  if (any(startsWith(tolower(output), "a"))) {
    p.tp <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut)
      }, lower = truecut, upper = 1)$value
    p.fp <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut)
      }, lower = 0, upper = truecut)$value
    p.ff <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE)
      }, lower = truecut, upper = 1)$value
    p.tf <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE)
      }, lower = 0, upper = truecut)$value
    camat <- confmat(p.tf, p.tp, p.ff, p.fp, "prop")
    out[["confusionmatrix"]] <- camat
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(startsWith(tolower(output), "c"))) {
    p.ii <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE) *
          stats::pbinom(floor(cut * N) - 1, N, x, lower.tail = TRUE )
      }, lower = 0, upper = 1)$value
    p.ij <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE) *
        (1 - stats::pbinom(floor(cut * N) - 1, N, x, lower.tail = TRUE))
      }, lower = 0, upper = 1)$value
    p.jj <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = FALSE) *
        (1 - stats::pbinom(floor(cut * N) - 1, N, x, lower.tail = TRUE))
      }, lower = 0, upper = 1)$value
    p.ji <- stats::integrate(function(x) {
      dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = FALSE) *
        stats::pbinom(floor(cut * N) - 1, N, x, lower.tail = TRUE)
      }, lower = 0, upper = 1)$value
    ccmat <- base::matrix(nrow = 2, ncol = 2, dimnames = list(c("i", "j"), c("i", "j")))
    ccmat["i", "i"] <- p.ii
    ccmat["i", "j"] <- p.ij
    ccmat["j", "i"] <- p.ji
    ccmat["j", "j"] <- p.jj
    out[["consistencymatrix"]] <- ccmat / sum(ccmat)
    out[["classification.consistency"]] <- ccStats(ccmat["i", "i"], ccmat["i", "j"], ccmat["j", "i"], ccmat["j", "j"])
  }
  base::return(out)
}

#' An Extension of the Livingston and Lewis (1995) Approach to Estimate Classification Consistency and Accuracy for Multiple Classifications based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{etl} for the effective test length (see documentation for the \code{ETL} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param min The minimum value possible to attain on the test. Default is 0.
#' @param max The maximum value possible to attain on the test. Default is 1 (assumes that the values contained in \code{x} represents proportions of maximum credit).
#' @param cut A vector of cut-off values for classifying observations into two or more categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit. The argument takes either a vector of two values, or \code{NULL}. If set to \code{NULL}, the model-fit test is not executed. If a vector of values is supplied, the first value is to represent the initial number of bins the distribution of scores is to be divided in to. This value is set to a default of 100. If this default results in too few bins to conduct the chi-square test, this value can be made larger. The second value represents the minimum expected number of observations that the bins should consist of. In accordance with standard recommendations for chi-square tests, the default value is set to 10.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false positive/negative categorizations for a test, diagnostic performance statistics, and/or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note It should be noted that this implementation differs from the original articulation of Livingston and Lewis (1995) in some respects. First, the procedure includes a number of diagnostic performance (accuracy) indices which the original procedure enables but that were not included. Second, the way consistency is calculated differs substantially from the original articulation of the procedure, which made use of a split-half approach. Rather, this implementation uses the approach to estimating classification consistency outlined by Hanson (1991).
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.1, 0.95, 5, 3)
#' for (i in 1:100) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 1000, ncol = 100)
#'   }
#'   rawdata[, i] <- rbinom(1000, 1, p.success)
#' }
#'
#' # Suppose the cutoff value for being placed in the lower category is a score
#' # below 50, second lowest 60, then 70, 80, and 90. Using the cba() function
#' # to estimate the reliability of this test, to use the LL.CA.MC() function
#' # or estimating diagnostic performance and consistency indices of
#' # classifications when using several cut-points:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)
#'
#' # The output from this function can get quite verbose when operating with
#' # several cut-points. In order to retrieve only model parameter estimates:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$parameters
#'
#' # To retrieve only the model-fit estimate:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$modelfit
#'
#' # To retrieve only the diagnostic performance estimates:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$accuracy
#'
#' # To retrieve only the classification consistency indices:
#' LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100)$consistency
#'
#' # Alternatively, the MC.out.tabular() function can be used to organize the
#' # category-specific indices in a tabular format:
#' MC.out.tabular(LL.CA.MC(rowSums(rawdata), cba(rawdata), c(50, 60, 70, 80, 90), min = 0, max = 100))
#'
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
LL.CA.MC <- function(x = NULL, reliability, cut, min = 0, max = 1, true.model = "4P", failsafe = TRUE, l = 0, u = 1, modelfit = c("nbins" = 100, "minbin" = 10)) {
  out <- base::list()
  if (!is.list(x)) {
    if ((base::min(x) < min) | (base::max(x) > max)) {
      warning(paste("Observed values not within the specified [", min, ", ", max, "] bounds (observed min = ",
                    base::min(x), ", observed max = ", base::max(x), ").", sep = ""))
    }
    N <- ETL(base::mean(x), stats::var(x), min = min, max = max, reliability = reliability)
    if (startsWith(as.character(true.model), "2")) {
      failsafe <- FALSE
    }
    params <- Beta.tp.fit(x, min = min, max = max, etl = N, true.model = true.model, failsafe = failsafe, l = l, u = u)
    if (params$l < 0 | params$u > 1) {
      warning(paste("Parameter out of bounds: l = ", round(params$l, 4), ", u = ", round(params$u, 4), ", alpha = ", round(params$alpha, 4), ", beta = ", round(params$beta, 4),
                    ". Consider constraining the fitting procedure further (e.g., set the location-parameters).", sep = ""))
    }
    x <- (x - min) / (max - min)
  } else {
    params <- x
    N <- params$etl
  }
  cut <- (cut - min) / (max - min)
  N <- base::round(N)
  params[["etl_rounded"]] <- N
  out[["parameters"]] <- params
  pcut <- c(0, cut, 1)
  ocut <- round(pcut * N)
  camat <- matrix(ncol = length(cut) + 1, nrow = length(cut) + 1)
  ccmat <- camat
  if (!is.list(x) & !is.null(modelfit)) {
    x <- round(x * N)
    tcut <- seq(0, N, N / modelfit[1])
    mdlfit <- matrix(nrow = 2, ncol = length(tcut) - 1)
    for (j in 1:(length(tcut) - 1)) {
      mdlfit[1, j] <- stats::integrate(function(x) {
        if (j == 1) {
          dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(tcut[j + 1] - 1, N, x))
        }
        else {
          if (j != 1 & j != (length(tcut) - 1)) {
            dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(tcut[j + 1] - 1, N, x) - stats::pbinom(tcut[j] - 1, N, x))
          } else {
            dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(tcut[j + 1], N, x) - stats::pbinom(tcut[j] - 1, N, x))
          }
        }
      }, lower = 0, upper = 1)$value
      if (j == ncol(mdlfit)) {
        mdlfit[2, j] <- length(x[x >= tcut[j]])
      } else {
        mdlfit[2, j] <- length(x[x < tcut[j + 1] & x >= tcut[j]])
      }
    }
    mdlfit[1, ] <- (mdlfit[1, ]/sum(mdlfit[1, ])) * length(x)
    for (i in 1:ncol(mdlfit)) {
      if (i < ncol(mdlfit)) {
        if (any(mdlfit[1, i] < ncol(mdlfit))) {
          if (any(mdlfit[1, i] < modelfit[2])) {
            mdlfit[, i + 1] <- mdlfit[, i + 1] + mdlfit[, i]
            mdlfit[, i] <- NA
          }
        }
      }
    }
    mdlfit <- mdlfit[, apply(mdlfit, 2, function(x) { any(!is.na(x)) })]
    if (any(mdlfit[1, ncol(mdlfit)] < modelfit[2])) {
      mdlfit[, ncol(mdlfit) - 1] <- mdlfit[, ncol(mdlfit) - 1] + mdlfit[, ncol(mdlfit)]
      mdlfit <- mdlfit[, -ncol(mdlfit)]
    }
    chisquared <- sum(apply(mdlfit, 2, function(x) {
      (x[2] - x[1])^2/x[1]
    }))
    out[["modelfit"]] <- list()
    out[["modelfit"]][["contingencytable"]] <- mdlfit
    out[["modelfit"]][["chisquared"]] <- chisquared
    if (startsWith(as.character(true.model), "2")) {
      out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
    }
    else {
      if ((startsWith(as.character(true.model), "4") & failsafe == TRUE) & (out[["parameters"]]$l == l & out[["parameters"]]$u == u)) {
        out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
      } else {
        out[["modelfit"]][["df"]] <- ncol(mdlfit) - 4
      }
    }
    out[["modelfit"]][["pvalue"]] <- stats::pchisq(chisquared, ncol(mdlfit) - 4, lower.tail = FALSE)
  }
  for (i in 1:(length(cut) + 1)) {
    if (i == 1) {
      for (j in 1:(length(cut) + 1)) {
        if (j == 1) {
          rnam <- NULL
          cnam <- NULL
        }
        if (j != (length(cut) + 1)) {
          if (j == 1) {
            rnam[j] <- paste("Observed <", pcut[j + 1] * max)
            cnam[j] <- paste("True     <", pcut[j + 1] * max)
          } else {
            rnam[j] <- paste(" >=", pcut[j] * max, "& <", pcut[j + 1] * max)
            cnam[j] <- paste(" >=", pcut[j] * max, "& <", pcut[j + 1] * max)
          }
        } else {
          rnam[j] <- paste(" >=", pcut[j] * max)
          cnam[j] <- paste(" >=", pcut[j] * max)
        }
      }
      colnames(camat) <- cnam
      rownames(camat) <- rnam
    }
    for (j in 1:(length(cut) + 1)) {
      camat[j, i] <- stats::integrate(function(x) {
        if(j == 1) {
          dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1] - 1, N, x))
        } else {
          if (j != 1 & j != (length(cut) + 1)) {
            dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1] - 1, N, x) - stats::pbinom(ocut[j] - 1, N, x))
          } else {
            dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1], N, x) - stats::pbinom(ocut[j] - 1, N, x))
          }
        }
      }, lower = pcut[i], upper = pcut[i + 1])$value
    }
  }
  out[["accuracy"]] <- list()
  out[["accuracy"]][["overall"]] <- list()
  out[["accuracy"]][["overall"]][["confusionmatrix"]] <- camat / sum(camat)
  out[["accuracy"]][["overall"]][["accuracy"]] <- sum(diag(camat))
  out[["accuracy"]][["specific"]]
  caout <- list()
  for(i in 1:ncol(camat)) {
    FN <- sum(camat[-i, i])
    FP <- sum(camat[i, -i])
    TP <- camat[i, i]
    TN <- sum(camat[-i, -i])
    caout[[paste("Category.", i, sep = "")]] <- list()
    caout[[i]][["confusionmatrix"]] <- confmat(TP, TN, FP, FN)
    caout[[i]][["statistics"]] <- caStats(TP, TN, FP, FN)
  }
  out[["accuracy"]][["specific"]] <- caout
  for (i in 1:(length(cut) + 1)) {
    if (i == 1) {
      for (j in 1:(length(cut) + 1)) {
        if (j == 1) {
          rnam <- NULL
          cnam <- NULL
        }
        if (j != (length(cut) + 1)) {
          if (j == 1) {
            rnam[j] <- paste("2nd adm. <", pcut[j + 1] * max)
            cnam[j] <- paste("1st adm. <", pcut[j + 1] * max)
          } else {
            rnam[j] <- paste(" >=", pcut[j] * max, "& <", pcut[j + 1] * max)
            cnam[j] <- paste(" >=", pcut[j] * max, "& <", pcut[j + 1] * max)
          }
        } else {
          rnam[j] <- paste(" >=", pcut[j] * max)
          cnam[j] <- paste(" >=", pcut[j] * max)
        }
      }
      colnames(ccmat) <- cnam
      rownames(ccmat) <- rnam
    }
    for (j in 1:(length(cut) + 1)) {
      ccmat[j, i] <- stats::integrate(function(x) {
        if (j == 1) {
          if (i == 1) {
            dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(ocut[j + 1] - 1, N, x) * stats::pbinom(ocut[i + 1] - 1, N, x)
          } else {
            if (i == (length(cut) + 1)) {
              dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(ocut[j + 1] - 1, N, x) *
                (stats::pbinom(ocut[i + 1], N, x) - stats::pbinom(ocut[i] - 1, N, x))
            } else {
              dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(ocut[j + 1] - 1, N, x) *
                (stats::pbinom(ocut[i + 1] - 1, N, x) - stats::pbinom(ocut[i] - 1, N, x))
            }
          }
        } else {
          if (j == (length(cut) + 1)) {
            if (i == 1) {
              dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1], N, x) - stats::pbinom(ocut[j] - 1, N, x)) * stats::pbinom(ocut[i + 1] - 1, N, x)
            } else {
              if (i == (length(cut) + 1)) {
                dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1], N, x) - stats::pbinom(ocut[j] - 1, N, x)) *
                  (stats::pbinom(ocut[i + 1], N, x) - stats::pbinom(ocut[i] - 1, N, x))
              } else {
                dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1], N, x) - stats::pbinom(ocut[j] - 1, N, x)) *
                  (stats::pbinom(ocut[i + 1] - 1, N, x) - stats::pbinom(ocut[i] - 1, N, x))
              }
            }
          } else {
            if (i == 1) {
              dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1] - 1, N, x) - stats::pbinom(ocut[j] - 1, N, x)) * stats::pbinom(ocut[i + 1] - 1, N, x)
            } else {
              if (i == (length(cut) + 1)) {
                dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1] - 1, N, x) - stats::pbinom(ocut[j] - 1, N, x)) *
                  (stats::pbinom(ocut[i + 1], N, x) - stats::pbinom(ocut[i] - 1, N, x))
              } else {
                dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (stats::pbinom(ocut[j + 1] - 1, N, x) - stats::pbinom(ocut[j] - 1, N, x)) *
                  (stats::pbinom(ocut[i + 1] - 1, N, x) - stats::pbinom(ocut[i] - 1, N, x))
              }
            }
          }
        }
      }, lower = 0, upper = 1)$value
    }
  }
  out[["consistency"]] <- list()
  out[["consistency"]][["overall"]] <- list()
  out[["consistency"]][["overall"]][["consistencymatrix"]] <- ccmat / sum(ccmat)
  p <- sum(diag(ccmat))
  p_c <- sum(apply(ccmat, 2, function(x) {
    sum(x)^2
  }))
  Kappa <- (p - p_c) / (1 - p_c)
  out[["consistency"]][["overall"]][["statistics"]] <- list("p" = p, "p_c" = p_c, "Kappa" = Kappa)
  ccout <- list()
  for(i in 1:ncol(ccmat)) {
    p <- ccmat[i, i]
    p_c <- sum(ccmat[i, ])^2
    Kappa <- (p - p_c) / (1 - p_c)
    ccout[[paste("Category.", i, sep = "")]] <- list()
    ccout[[i]][["statistics"]] <- list("p" = p, "p_c" = p_c, "Kappa" = Kappa)
  }
  out[["consistency"]][["specific"]] <- ccout
  base::return(out)
}

#' Confusion matrix
#'
#' @description Organizes supplied values of true and false positives and negatives into a confusion matrix.
#' @param tp The frequency or rate of true-positive classifications.
#' @param tn The frequency or rate of true-negative classifications.
#' @param fp The frequency or rate of false-positive classifications.
#' @param fn The frequency or rate of false-negative classifications.
#' @param output Whether the returned output reflects frequencies or proportions. Defaults to returning frequencies.
#' @return A confusion matrix organizing the input values of true and false positive and negatives.
#' @examples
#' # Generate some true and observed conditions.
#' set.seed(1234)
#' true.ability <- rbeta(50, 4, 4)
#' true.category <- ifelse(true.ability < 0.5, 0, 1)
#' observed.score <- rbinom(50, 50, true.ability)
#' observed.category <- ifelse(observed.score < 25, 0, 1)
#' # Calculate the frequencies of true and false positives and negatives based on the true and
#' # observed conditions.
#' TP <- sum(ifelse(observed.category == 0 & true.category == 0, 1, 0))
#' FP <- sum(ifelse(observed.category == 0 & true.category != 0, 1, 0))
#' TN <- sum(ifelse(observed.category == 1 & true.category == 1, 1, 0))
#' FN <- sum(ifelse(observed.category == 1 & true.category != 1, 1, 0))
#' # Organize the above values in a confusion matrix using the confmat function:
#' confmat(tp = TP, fp = FP, tn = TN, fn = FN)
#' @export
confmat <- function(tp, tn, fp, fn, output = "freq") {
  mat <- base::matrix(nrow = 3, ncol = 3)
  base::rownames(mat) <- c("True", "False", "Total")
  base::colnames(mat) <- c("Positive", "Negative", "Total")
  tot <- base::sum(tp, tn, fp, fn)
  mat[1, 1] <- tp
  mat[1, 2] <- tn
  mat[2, 1] <- fp
  mat[2, 2] <- fn
  mat[, 3] <- base::rowSums(mat[, -3])
  mat[3, ] <- base::colSums(mat[-3, ])
  if (output != "freq") {
    mat <- mat / base::sum(tp, tn, fp, fn)
  }
  mat
}

#' Classification Accuracy Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the certainty of classifications based on tests.
#' @param tp The frequency or rate of true-positive classifications.
#' @param tn The frequency or rate of true-negative classifications.
#' @param fp The frequency or rate of false-positive classifications.
#' @param fn The frequency or rate of false-negative classifications.
#' @return A list of diagnostic performance statistics based on true/false positive/negative statistics. Specifically, the sensitivity, specificity, positive predictive value (PPV), negative predictive value (NPV), Youden's J. (Youden.J), and Accuracy.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. First, compute the
#' # estimated confusion matrix using LL.CA():
#' cmat <- LL.CA(x = testdata, reliability = 0.7, cut = 50, min = 0,
#' max = 100)$confusionmatrix
#'
#' # To estimate and retrieve diagnostic performance statistics using caStats(),
#' # feed it the appropriate entries of the confusion matrix.
#' caStats(tp = cmat["True", "Positive"], tn = cmat["True", "Negative"],
#' fp = cmat["False", "Positive"], fn = cmat["False", "Negative"])
#' @export
caStats <- function(tp, tn, fp, fn) {
sensitivity <-  tp / (tp + fn)
specificity <-  tn / (tn + fp)
ppv <-          tp / (tp + fp)
npv <-          tn / (tn + fn)
accuracy <-     (tp + tn) / (tp + tn + fp + fn)
J <-            (sensitivity + specificity) - 1
base::list("Sensitivity" = sensitivity, "Specificity" = specificity,
           "PPV" = ppv, "NPV" = npv,
           "Youden.J" = J, "Accuracy" = accuracy)
}

#' Classification Consistency Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the consistency of classifications based on tests.
#' @param ii The frequency or rate of consistent classifications into category "i".
#' @param ij The frequency or rate of inconsistent classifications into categories "i" and "j".
#' @param ji The frequency or rate of inconsistent classifications into categories "j" and "i".
#' @param jj The frequency or rate of consistent classifications into category "j".
#' @return A list of classification consistency statistics. Specifically, the coefficient of consistent classification (p), the coefficient of consistent classification by chance (p_c), the proportion of positive classifications due to chance (p_c_pos), the proportion of negative classifications due to chance (p_c_neg), and Cohen's Kappa coefficient.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. First, compute the
#' # estimated consistency matrix using LL.CA():
#' cmat <- LL.CA(x = testdata, reliability = .7, cut = 50, min = 0,
#' max = 100)$consistencymatrix
#'
#' # To estimate and retrieve consistency statistics using ccStats(),
#' # feed it the appropriate entries of the consistency matrix.
#' ccStats(ii = cmat["i", "i"], ij = cmat["i", "j"],
#' ji = cmat["j", "i"], jj = cmat["j", "j"])
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @export
ccStats <- function(ii, ij, ji, jj) {
  p <- (ii + jj) / (ii + ij + ji + jj)
  p_c <- (ii + ij) * (ii + ji) + (ij + jj) * (ji + jj)
  p_c_pos <- (ii + ij) * (ii + ji)
  p_c_neg <- (ij + jj) * (ji + jj)
  Kappa <- (p - p_c) / (1 - p_c)
  base::list("p" = p, "p_c" = p_c, "p_c_pos" = p_c_pos, "p_c_neg" = p_c_neg, "Kappa" = Kappa)
}

#' ROC curves for the Livingston and Lewis approach.
#'
#' @description Generate a ROC curve plotting the false-positive rate against the true-positive rate at different cut-off values across the observed-score scale.
#' @param x A vector of observed results.
#' @param min The minimum possible value to attain on the observed-score scale.
#' @param max The maximum value possible to attain on the test. Default is 1 (assumes that the values contained in \code{x} represents proportions of maximum credit).
#' @param reliability The reliability coefficient of the test.
#' @param truecut The true point along the x-scale that marks the categorization-threshold.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The \code{"4P"} method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the \code{"2P"} method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe If true-model == "4P": Whether to engage a fail-safe reverting to a two-parameter true-score distribution solution should the four-parameter fitting procedure produce impermissible results. Default is TRUE (engage fail-safe in the event of impermissible estimates).
#' @param l If \code{true.model == "2P"} or \code{failsafe == TRUE}: The lower-bound location parameter of the two-parameter true-score distribution solution.
#' @param u If \code{true.model == "2P"} or \code{failsafe == TRUE}: The upper-bound location parameter of the two-parameter true-score distribution solution.
#' @param AUC Logical. Calculate and include the area under the curve? Default is \code{FALSE}.
#' @param maxJ Logical. Mark the point along the curve where Youden's J statistic is maximized? Default is \code{FALSE}.
#' @param maxAcc Logical. Mark the point along the curve where the Accuracy statistic is maximized? Default is \code{FALSE}.
#' @param locate Ask the function to locate the cut-point at which sensitivity or NPV is greater than or equal to some value, or specificity or PPV is lesser than or equal to some value. Take as input a character-vector of length 2, with the first argument being which index is to be found (e.g., "sensitivity"), and the second argument the value to locate (e.g., "0.75"). For example: c("sensitivity", "0.75").
#' @param raw.out Give raw coordinates as output rather than plot? Default is \code{FALSE}
#' @param grainsize Specify the number of cutoff-points for which the ROC curve is to be calculated. The greater this number the greater the accuracy. Default is 100 points.
#' @return A plot tracing the ROC curve for the test, or matrix of coordinates if raw.out is \code{TRUE}.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' hist(testdata / 100, xlim = c(0, 1), freq = FALSE)
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct.
#' # Suppose further that the reliability of the test-scores were estimated to
#' # 0.75. To produce a plot with an ROC curve using LL.ROC(), along with the
#' # AUC statistics and the points at which Youden's J. is maximized:
#' LL.ROC(x = testdata, reliability = 0.7, truecut = 50, min = 0, max = 100,
#' AUC = TRUE, maxJ = TRUE)
#' # Or to locate the point at which accuracy is maximized:
#' LL.ROC(x = testdata, reliability = 0.7, truecut = 50, min = 0, max = 100,
#' maxAcc = TRUE)
#'
#' # Using the example data above, the function can be instructed to locate an
#' # operational cut-point at which sensitivity or specificity is equal to or
#' # greater than some specified value by specifying the "locate" argument with
#' # c("statistic", value). For example, to locate the operational cut-point at
#' # which sensitivity is first equal to or greater than 0.9:
#' LL.ROC(testdata, reliability = 0.7, min = 0, max = 100, truecut = 50,
#' locate = c("sensitivity", 0.9))
#' # For Negative Predictive value, the point at which it is equal or greater:
#' LL.ROC(testdata, reliability = 0.7, min = 0, max = 100, truecut = 50,
#' locate = c("NPV", 0.9))
#' # And so on for other statistics such as Specificity and Positive Predictive
#' # Value.
#' @export
LL.ROC <- function(x = NULL, reliability, min = 0, max = 1, truecut, true.model = "4P", failsafe = TRUE, l = 0, u = 1, AUC = FALSE, maxJ = FALSE, maxAcc = FALSE, locate = NULL, raw.out = FALSE, grainsize = 100) {
  if (!raw.out) {
    oldpar <- graphics::par(no.readonly = TRUE)
    base::on.exit(graphics::par(oldpar))
  }
  if (!is.list(x)) {
    x <- Beta.tp.fit(x, min, max, reliability = reliability, true.model = true.model, failsafe = failsafe, l = l, u = u)
  }
  for (i in 1:(grainsize + 1)) {
    if (i == 1) {
      cuts <- seq(min, max, (max - min) / grainsize)
      outputmatrix <- matrix(nrow = grainsize + 1, ncol = 7)
      outputmatrix[, 4] <- cuts
    }
    axval <- LL.CA(x = x, min = min, max = max, cut = cuts[i],
                   truecut = truecut, true.model = true.model,
                   output = "a", l = l, u = u, modelfit = NULL)$classification.accuracy
    outputmatrix[i, 1] <- 1 - axval$Specificity
    outputmatrix[i, 2] <- axval$Sensitivity
    outputmatrix[i, 3] <- axval$Youden.J
    outputmatrix[i, 5] <- axval$Accuracy
    outputmatrix[i, 6] <- axval$PPV
    outputmatrix[i, 7] <- axval$NPV
    base::colnames(outputmatrix) <- c("FPR", "TPR", "Youden.J", "Cutoff", "Accuracy", "PPV", "NPV")
    outputmatrix[base::which(base::is.na(outputmatrix[, 1])), 1] <- 0
    outputmatrix[base::which(base::is.na(outputmatrix[, 2])), 2] <- 1
    outputmatrix[base::which(base::is.na(outputmatrix[, 6])), 6] <- 1
    outputmatrix[base::which(base::is.na(outputmatrix[, 7])), 7] <- 1
  }
  if (raw.out) {
    base::return(outputmatrix)
  }
  graphics::plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
  graphics::abline(h = seq(0, 1, .1), v = seq(0, 1, .1), col = "lightgrey", lty = "dotted")
  graphics::par(new = TRUE)
  graphics::plot(outputmatrix[, 1], outputmatrix[, 2], type = "s",
       xlab = "False-Positive Rate (1 - Specificity)",
       ylab = "True-Positive Rate (Sensitivity)",
       main = paste("ROC curve for true-cut equal to", truecut), lwd = 2,
       xlim = c(0, 1), ylim = c(0, 1))
  if (AUC) {
    graphics::legend("bottomright", bty = "n", cex = 1.5,
                     legend = paste("AUC =", round(AUC(outputmatrix[, 1], outputmatrix[, 2]), 3)))
  }
  if (maxJ) {
    graphics::points(outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 1],
           outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 2], cex = 1.5, pch = 19)
    graphics::text(outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 1] + .025,
         outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 2] - .025,
         labels = base::paste("Maximum Youden's J. ", "(", base::round(base::max(outputmatrix[, 3]), 3), ") ",  "at cut-off = ",
                        base::round(outputmatrix[which(outputmatrix[, 3] == base::max(outputmatrix[, 3]))[1], 4], 3), ".", sep = ""),
         adj = c(0, 1))
  }
  if (maxAcc) {
    graphics::points(outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 1],
                     outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 2], cex = 1.5, pch = 19)
    graphics::text(outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 1] + .025,
                   outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 2] - .025,
                   labels = base::paste("Maximum Accuracy ", "(", base::round(base::max(outputmatrix[, 5]), 3), ") ",  "at cut-off = ",
                                  base::round(outputmatrix[which(outputmatrix[, 5] == base::max(outputmatrix[, 5]))[1], 4], 3), ".", sep = ""),
                   adj = c(0, 1))
  }
  if (!base::is.null(locate)) {
    if (base::startsWith(base::tolower(locate[1]), "se")) {
      rowloc <- base::which(outputmatrix[, "TPR"] >= base::as.numeric(locate[2]))[1]
      colloc <- "TPR"
      value <- outputmatrix[rowloc, "TPR"]
      statistic <- "Sensitivity >= "
    }
    if (base::startsWith(base::tolower(locate[1]), "p")) {
      rowloc <- base::which(outputmatrix[, "PPV"] <= base::as.numeric(locate[2]))[1]
      colloc <- "PPV"
      value <- outputmatrix[rowloc, "PPV"]
      statistic <- "PPV <= "
    }
    if (base::startsWith(base::tolower(locate[1]), "n")) {
      rowloc <- base::which(outputmatrix[, "NPV"] >= base::as.numeric(locate[2]))[1]
      colloc <- "NPV"
      value <- outputmatrix[rowloc, "NPV"]
      statistic <- "NPV >= "
    }
    if (base::startsWith(base::tolower(locate[1]), "sp")) {
      rowloc <- base::which(outputmatrix[, "FPR"] <= 1 - base::as.numeric(locate[2]))[base::length(base::which(outputmatrix[, "FPR"] <= 1 - base::as.numeric(locate[2])))]
      colloc <- "FPR"
      value <- 1 - outputmatrix[rowloc, "FPR"]
      statistic <- "Specificity <= "
    }
    graphics::points(outputmatrix[rowloc, 1], outputmatrix[rowloc, 2], cex = 1.5, pch = 19)
    graphics::text(outputmatrix[rowloc, 1] + .025,
                   outputmatrix[rowloc, 2] - .025,
                   labels = base::paste(statistic, base::as.numeric(locate[2]) , " (", base::round(value, 3), ")",
                                  " at cut-off = ", outputmatrix[rowloc, 4], ".", sep = ""),
                   adj = if (base::startsWith(base::tolower(locate[1]), "sp")) { c(0, 1) } else { c(0, 1) })
  }
}

#' Area Under the ROC Curve.
#'
#' @description Given a vector of false-positive rates and a vector of true-positive rates, calculate the area under the Receiver Operator Characteristic (ROC) curve.
#' @param FPR Vector of False-Positive Rates.
#' @param TPR Vector of True-Positive Rates.
#' @return A value representing the area under the ROC curve.
#' @note Script originally retrieved and modified from https://blog.revolutionanalytics.com/2016/11/calculating-auc.html.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To calculate the
#' # necessary (x, y) coordinates to compute the area under the curve statistic
#' # one can use the LL.ROC() function with the argument
#' # raw.out = TRUE.
#' coords <- LL.ROC(x = testdata, reliability = .7, truecut = 50, min = 0,
#' max = 100, raw.out = TRUE)
#'
#' # To calculate and retrieve the Area Under the Curve (AUC) with the AUC()
#' # function, feed it the raw coordinates calculated above.
#' AUC(coords[, "FPR"], coords[, "TPR"])
#' @export
AUC <- function(FPR, TPR) {
  dFPR <- base::c(base::diff(FPR), 0)
  dTPR <- base::c(base::diff(TPR), 0)
  base::sum(TPR * dFPR) + base::sum(dTPR * dFPR) / 2
}

#' Calculate Cronbach's Alpha reliability-coefficient from supplied variables.
#'
#' @description Calculates Cronbach's Alpha reliability coefficient of the sum-score.
#' @param x A data-frame or matrix of numerical values where rows represent respondents, and columns represent items.
#' @note Missing values are treated by passing \code{na.rm = TRUE} to the \code{var} function call.
#' @note Be aware that this function does not issue a warning if there are negative correlations between variables in the supplied data-set.
#' @return Cronbach's Alpha for the sum-score of supplied variables.
#' @references Cronbach, L.J. (1951). Coefficient alpha and the internal structure of tests. Psychometrika 16, 297--334. doi: 10.1007/BF02310555
#' @examples
#' # Generate some fictional data. Say 100 students take a 50-item long test
#' # where all items are equally difficult.
#' set.seed(1234)
#' p.success <- rBeta.4P(100, 0.25, 0.75, 5, 3)
#' for (i in 1:50) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 100, ncol = 50)
#'   }
#'   rawdata[, i] <- rbinom(100, 1, p.success)
#' }
#' # To calculate Cronbach's Alpha for this test:
#' cba(rawdata)
#' @export
cba <- function(x) {
  (base::ncol(x) / (base::ncol(x) - 1)) *
    (1 - (base::sum(base::diag(stats::var(x, na.rm = TRUE))) /
            base::sum(stats::var(x, na.rm = TRUE))))
}

#' Calculate McDonald's Omega reliability-coefficient from supplied variables.
#'
#' @description Calculates McDonalds's Omega reliability-coefficient of the sum-score from the Spearman one-factor model using the procedure outlined in McDonald (1999).
#' @param x A data-frame or matrix of numerical values where rows represent respondents, and columns represent items.
#' @param fit Logical. Default is \code{FALSE}. If \code{TRUE}, the output changes from a vector containing the Omega reliability-estimate to a list containing additional detailed information concerning the fitted factor model.
#' @note Missing values are treated by passing \code{na.rm = TRUE} to the \code{var} function call and \code{use = "pairwise.complete.obs"} to the \code{cov} function call.
#' @note The function terminates with an error if there are negative covariance-matrix entries.
#' @return If \code{fit = FALSE}, A vector of length 1 containing the estimated McDonalds's Omega reliability-coefficient for the sum-score of the supplied variables. If \code{fit = TRUE}, a list containing the Omega-coefficient reliability-estimate as the first entry, followed by the goodness-of-fit index (GFI), a two-row matrix containing the estimated factor-loadings and error-variances, and the observed and fitted covariance-matrices and the discrepancy matrix.
#' @references McDonald, R. P. (1999). Test Theory: A Unified Treatment. Routledge.
#' @examples
#' # Generate some fictional data.
#' set.seed(1234)
#' rawdata <- matrix(rnorm(500), ncol = 5)
#' common <- rnorm(100)
#' rawdata <- apply(rawdata, 2, function(x) {x + common})
#'
#' # To estimate McDonald's Omega from this data:
#' mdo(rawdata)
#'
#' # To retrieve additional information such as the GFI fit-index and model-
#' # parameter estimates:
#' mdo(rawdata, fit = TRUE)
#' @export
mdo <- function(x, fit = FALSE) {
  vars <- base::apply(x, 2, stats::var, na.rm = TRUE)
  covs <- base::list()
  if (base::any(stats::cov(x, use = "pairwise.complete.obs") <= 0)) {
    stop("Item-covariance(s) less than or equal to 0 detected. Consider reverse-scoring or excluding variables.")
  }
  for (i in 1:base::ncol(x)) {
    covs[[i]] <- base::list()
    if (i == 1) {
      y <- stats::cov(x, use = "pairwise.complete.obs")
    } else {
      y <- base::rbind(y, y[1, ]); y <- base::cbind(y, y[, 1]); y <- y[-1, ]; y <- y[, -1]
    }
    for (j in 1:(base::ncol(x) - 1)) {
      covs[[i]][[j]] <- y[(j + 1):base::ncol(y), j]
    }
  }
  fload <- base::lapply(covs, function(y) {
    for (i in 1:(base::length(y) - 1)) {
      for (j in 1:base::length(y[[i + 1]])) {
        if (i == 1 & j == 1)  {
          out <- base::sqrt((y[[1]][i] * y[[1]][i + j]) /  y[[i + 1]][j])
        } else {
          out[base::length(out) + 1] <- base::sqrt((y[[1]][i] * y[[1]][i + j]) /  y[[i + 1]][j])
        }
      }
    }
    base::mean(out)
  })
  for (i in 1:base::length(fload)) {
    if (i == 1) {
      floads <- fload[[1]][1]
    } else {
      floads[i] <- fload[[i]][1]
    }
  }
  if (!fit) {
    base::return(base::sum(floads)^2 / (base::sum(vars - floads^2) + base::sum(floads)^2))
  } else {
    Omega <- base::sum(floads)^2 / (base::sum(vars - floads^2) + base::sum(floads)^2)
    omat <- stats::cov(x, use = "pairwise.complete.obs")
    evars <- vars - floads^2
    for (i in 1:base::length(floads)) {
      for (j in 1:base::length(floads)) {
        if (i == 1 & j == 1) {
          imat <- base::matrix(ncol = base::length(floads), nrow = base::length(floads))
        }
        if (i == j) {
          imat[j, i] <- floads[j]^2 + evars[j]
        } else {
          imat[j, i] <- floads[i] * floads[j]
        }
      }
    }
    dmat <- omat - imat
    gfi <- 1 - (base::mean(dmat^2) / base::mean(omat))
    floads <- base::matrix(base::c(floads, evars), nrow = 2, byrow = TRUE, dimnames = base::list(c("Loadings", "Errors")))
    if (!base::is.null(base::colnames(x))) {
      base::colnames(floads) <- base::rownames(imat) <- base::colnames(imat) <- base::colnames(x)
    }
    base::return(base::list("Omega" = Omega,
                            "GFI" = gfi,
                            "Parameters" = floads,
                            "Matrices" = base::list("Observed" = omat,
                                                    "Fitted" = imat,
                                                    "Discrepancy" = dmat)))
  }
}

#' Estimate Beta true-score distribution based on observed-score raw-moments and the effective test length.
#'
#' @description Estimator for the Beta true-score distribution shape-parameters from the observed-score distribution and Livingston and Lewis' effective test length. Returns a list with entries representing the lower- and upper shape parameters (l and u), and the shape parameters (alpha and beta) of the four-parameters beta distribution, and the effective test length.
#' @param x Vector of observed-scores.
#' @param min The minimum possible score to attain on the test.
#' @param max The maximum possible score to attain on the test.
#' @param etl The value of Livingston and Lewis' effective test length. See \code{?ETL()}. Not necessary to specify if reliability is supplied to the \code{reliability} argument.
#' @param reliability Optional specification of the test-score reliability coefficient. If specified, overrides the input of the \code{etl} argument.
#' @param true.model The type of Beta distribution which is to be fit to the moments of the true-score distribution. Options are \code{"4P"} and \code{"2P"}, where \code{"4P"} refers to the four-parameter (with the same mean, variance, skewness, and kurtosis), and \code{"2P"} the two-parameter solution where both location-parameters are specified (with the same mean and variance).
#' @param failsafe Logical. Whether to revert to a fail-safe two-parameter solution should the four-parameter solution contain invalid parameter estimates.
#' @param l If \code{failsafe = TRUE} or \code{true.model = "2P"}: The lower-bound of the Beta distribution. Default is 0 (i.e., the lower-bound of the Standard, two-parameter Beta distribution).
#' @param u If \code{failsafe = TRUE} or \code{true.model = "2P"}: The upper-bound of the Beta distribution. Default is 1 (i.e., the upper-bound of the Standard, two-parameter Beta distribution).
#' @param output Option to specify true-score distribution moments as output if the value of the output argument does not equal \code{"parameters"}.
#' @return A list with the parameter values of a four-parameter Beta distribution. "l" is the lower location-parameter, "u" the upper location-parameter, "alpha" the first shape-parameter, "beta" the second shape-parameter, and "etl" the effective test length.
#' @references Hanson, B. A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series. Retrieved from https://files.eric.ed.gov/fulltext/ED344945.pdf
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @references Rogosa, D. &  Finkelman, M. (2004). How Accurate Are the STAR Scores for Individual Students? An Interpretive Guide. Retrieved from http://statweb.stanford.edu/~rag/accguide/guide04.pdf
#' @examples
#' # Generate some fictional data. Say 1000 individuals take a 100-item test
#' # where all items are equally difficult, and the true-score distribution
#' # is a four-parameter Beta distribution with location parameters l = 0.25,
#' # u = 0.75, alpha = 5, and beta = 3:
#' set.seed(12)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#'
#' # Since this test contains items which are all equally difficult, the true
#' # effective test length (etl) is the actual test length. I.e., etl = 100.
#' # To estimate the four-parameter Beta distribution parameters underlying
#' # the draws from the binomial distribution:
#' Beta.tp.fit(testdata, 0, 100, 100)
#'
#' # Imagine a case where the fitting procedure produces an impermissible
#' # estimate (e.g., l < 0 or u > 1).
#' set.seed(1234)
#' testdata <- rbinom(1000, 50, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' Beta.tp.fit(testdata, 0, 50, 50)
#'
#' # This example produced an l-value estimate less than 0. One way of
#' # dealing with such an occurrence is to revert to a two-parameter
#' # model, specifying the l and u parameters and estimating the
#' # alpha and beta parameters necessary to produce a Beta distribution
#' # with the same mean and variance as the estimated true-score distribution.
#'
#' # Suppose you have good theoretical reasons to fix the l parameter at a
#' # value of 0.25 (e.g., the test is composed of multiple-choice questions
#' # with four response-options, resulting in a 25% chance of guessing the
#' # correct answer). The l-parameter could be specified to this theoretically
#' # justified value, and the u-parameter could be specified to be equal to the
#' # estimate above (u = 0.7256552) as such:
#' Beta.tp.fit(testdata, 0, 50, 50, true.model = "2P", l = 0.25, u = 0.7256552)
#' @export
Beta.tp.fit <- function(x, min, max, etl = NULL, reliability = NULL, true.model = "4P", failsafe = FALSE, l = 0, u = 1, output = "parameters") {
  if (is.null(etl)) {
    etl <- ETL(base::mean(x), stats::var(x), min, max, reliability)
  }
  x <- (x - min)/(max - min) * etl
  m <- HB.tsm(x, 4, etl, 0)
  s2 <- m[2] - m[1]^2
  s3 <- (m[3] - 3 * m[1] * m[2] + 2 * m[1]^3)
  s4 <- (m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4)
  g3 <- (m[3] - 3 * m[1] * m[2] + 2 * m[1]^3) / (sqrt(s2)^3)
  g4 <- (m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4) / (sqrt(s2)^4)
  if (output == "parameters") {
    if (startsWith(as.character(true.model), "4")) {
      out <- Beta.4p.fit(mean = m[1], variance = s2, skewness = g3, kurtosis = g4)
      if (failsafe == TRUE & (out$l < 0 | out$u > 1)) {
        warning(paste("Fail-safe engaged: l = ", out$l, ", u = ", out$u, ", alpha = ", out$alpha, ", beta = ", out$beta,
                      ". \n  Finding permissible solution for the true-score distribution in accordance with specifications.", sep = ""))
        out <- Beta.2p.fit(mean = m[1], variance = s2, l = l, u = u)
      }
    }
    if (startsWith(as.character(true.model), "2")) {
      out <- Beta.2p.fit(mean = m[1], variance = s2, l = l, u = u)
    }
    out[["etl"]] <- etl
    return(out)
  } else {
    moments <- list()
    moments[["Raw"]] <- base::list(m[1], m[2], m[3], m[4])
    moments[["Central"]] <- base::list(0, s2, s3, s4)
    moments[["Standardized"]] <- base::list(0, 1, g3, g4)
    base::return(moments)
  }
}

#' Estimate Beta True-Score Distribution Based on Observed-Score Raw-Moments and Lord's k.
#'
#' @description Estimator for the Beta true-score distribution shape-parameters from the observed-score distribution and Lord's k. Returns a list with entries representing the lower- and upper shape parameters (l and u), and the shape parameters (alpha and beta) of the four-parameters beta distribution, as well as Lord's k and the test length.
#' @param x Vector of observed-scores.
#' @param N The test length.
#' @param k Lord's k (see documentation for the \code{Lords.k()} function).
#' @param true.model The type of Beta distribution which is to be fit to the moments of the true-score distribution. Options are \code{"4P"} and \code{"2P"}, where \code{"4P"} refers to the four-parameter (with the same mean, variance, skewness, and kurtosis), and \code{"2P"} the two-parameter solution where both location-parameters are specified (with the same mean and variance).
#' @param failsafe Logical. Whether to revert to a fail-safe two-parameter solution should the four-parameter solution contain invalid parameter estimates.
#' @param l If \code{failsafe = TRUE} or \code{true.model = "2P"}: The lower-bound of the Beta distribution. Default is 0 (i.e., the lower-bound of the Standard, two-parameter Beta distribution).
#' @param u If \code{failsafe = TRUE} or \code{true.model = "2P"}: The upper-bound of the Beta distribution. Default is 1 (i.e., the upper-bound of the Standard, two-parameter Beta distribution).
#' @return A list with the parameter values of a four-parameter Beta distribution. "l" is the lower location-parameter, "u" the upper location-parameter, "alpha" the first shape-parameter, and "beta" the second shape-parameter. Also includes Lord's k and the test length.
#' @references Hanson, B. A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series. Retrieved from https://files.eric.ed.gov/fulltext/ED344945.pdf
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @examples
#' # Generate some fictional data. Say 1000 individuals take a 100-item test
#' # where all items are equally difficult, and the true-score distribution
#' # is a four-parameter Beta distribution with location parameters l = 0.25,
#' # u = 0.75, alpha = 5, and beta = 3, and the error distribution is Binomial
#' # with Lord's k = 0:
#' set.seed(12)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#'
#' # To estimate the four-parameter Beta distribution parameters from this
#' # sample of observations:
#' HB.beta.tp.fit(testdata, 100, 0)
#' @export
HB.beta.tp.fit <- function(x, N, k, true.model = "4P", failsafe = FALSE, l = 0, u = 1) {
  m <- HB.tsm(x, 4, N, k)
  s2 <- m[2] - m[1]^2
  g3 <- (m[3] - 3 * m[1] * m[2] + 2 * m[1]^3) / (sqrt(s2)^3)
  g4 <- (m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4) / (sqrt(s2)^4)
  if (startsWith(as.character(true.model), "4")) {
    out <- Beta.4p.fit(mean = m[1], variance = s2, skewness = g3, kurtosis = g4)
    if (failsafe == TRUE & (out$l < 0 | out$u > 1)) {
      warning(paste("Fail-safe engaged: l = ", out$l, ", u = ", out$u, ", alpha = ", out$alpha, ", beta = ", out$beta,
                    ". \n  Finding permissible solution for the true-score distribution in accordance with specifications.", sep = ""))
      out <- Beta.2p.fit(mean = m[1], variance = s2, l = l, u = u)
    }
  }
  if (startsWith(as.character(true.model), "2")) {
    out <- Beta.2p.fit(mean = m[1], variance = s2, l = l, u = u)
  }
  out[["k"]] <- k
  out[["N"]] <- N
  out
}

#' Descending (falling) factorial.
#'
#' @description Calculate the descending (or falling) factorial of a value \code{x} of order \code{r}.
#' @param x A value for which the descending factorial is to be calculated.
#' @param r The power \code{x} is to be raised to.
#' @return The descending factorial of value \code{x} raised to the \code{r}'th power.
#' @param method The method by which the descending factorials are to be calculated. Default is \code{"product"} which uses direct arithmetic. Alternative is \code{"gamma"} which calculates the ascending factorial using the Gamma function. The alternative method might be faster but might fail because the Gamma function is not defined for negative integers (returning Inf).
#' @export
#' @examples
#' # To calculate the 4th descending factorial for a value (e.g., 3.14):
#' dfac(x = 3.14, r = 4)
#'
#' # To calculate the 5th descending factorial for values 3.14, 2.72, and 0.58:
#' dfac(x = c(3.14, 2.72, 0.58), r = 5)
dfac <- function(x, r, method = "product") {
  if (method == "product") {
    x <- base::ifelse(x < r, 0, x)
    if (r <= 1) {
      x^r
    } else {
      mat <- base::matrix(nrow = length(x), ncol = r)
      for (i in 1:r) {
        if (i == 1) {
          mat[, 1] <- x
        } else {
          mat[, i] <- x - i + 1
        }
      }
      base::apply(mat, 1, prod)
    }
  } else {
    base::gamma(x + 1) / base::gamma(x - r + 1)
  }
}

#' Ascending (rising) factorial.
#'
#' @description Calculate the ascending (or rising) factorial of a value \code{x} of order \code{r}.
#' @param x A value for which the ascending factorial is to be calculated.
#' @param r The power \code{x} is to be raised to.
#' @param method The method by which the descending factorials are to be calculated. Default is \code{"product"} which uses direct arithmetic. Alternative is \code{"gamma"} which calculates the descending factorial using the Gamma function. The alternative method might be faster but might fail because the Gamma function is not defined for negative integers (returning Inf).
#' @return The ascending factorial of value \code{x} raised to the \code{r}'th power.
#' @export
#' @examples
#' # To calculate the 4th ascending factorial for a value (e.g., 3.14):
#' afac(x = 3.14, r = 4)
#'
#' # To calculate the 5th ascending factorial for values 3.14, 2.72, and 0.58:
#' afac(x = c(3.14, 2.72, 0.58), r = 5)
afac <- function(x, r, method = "product") {
  if (method == "product") {
    if (r <= 1) {
      x^r
    } else {
      mat <- base::matrix(nrow = length(x), ncol = r)
      for (i in 1:r) {
        if (i == 1) {
          mat[, 1] <- x
        } else {
          mat[, i] <- x + i - 1
        }
      }
      base::apply(mat, 1, prod)
    }
  } else {
    base::gamma(x + r) / base::gamma(x)
  }
}

#' Proportional true-score distribution raw moments from Livingston and Lewis' effective test-score and effective test-length.
#'
#' @description An implementation of Lords (1965, p. 265) equation 37 for estimating the raw moments of the true-score distribution, modified to function for the Livingston and Lewis approach.
#' @param x The effective test-score of test-takers.
#' @param r The moment-order that is to be calculated (where 1 is the mean, 2 is the raw variance, 3 is the raw skewness, etc.).
#' @param n The effective test-length.
#' @param method The method by which the descending factorials are to be calculated. Default is \code{"product"} which uses direct arithmetic. Alternative is "gamma" which calculates the descending factorial using the Gamma function. The alternative method might be faster but might fail because the Gamma function is not defined for negative integers (returning Inf).
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
#' @examples
#' # Examine the raw moments of the underlying Beta distribution that is to provide the basis for
#' # observed-scores:
#' betamoments(alpha = 5, beta = 3, l = 0.25, u = 0.75, types = "raw")
#'
#' # Generate observed-scores from true-scores by passing the true-scores as binomial probabilities
#' # for the rbinom function.
#' set.seed(1234)
#' obs.scores <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
#' # Examine the raw moments of the observed-score distribution.
#' observedmoments(obs.scores, type = "raw")
#'
#' # First four estimated raw moment of the proportional true-score distribution from the observed-
#' # score distribution. As all items are equally difficult, the effective test-length is equal to
#' # the actual test-length.
#' tsm(x = obs.scores, r = 1, n = 100)
#' tsm(x = obs.scores, r = 2, n = 100)
#' tsm(x = obs.scores, r = 3, n = 100)
#' tsm(x = obs.scores, r = 4, n = 100)
#' # Which is fairly close to the true raw moments of the proportional true-score distribution
#' # calculated above.
tsm <- function(x, r, n, method = "product") {
  if (method != "product") {
    base::mean(dfac(x, r, method)) / base::mean(dfac(n - 2, r - 2, method)) / dfac(n, 2, method)
  } else {
    if (r == 1) {
      base::mean(x) / n
    } else {
      (base::mean(dfac(x, r)) / dfac(n - 2, r - 2)) * (1 / dfac(n, 2))
    }
  }
}

#' Proportional True-Score Distribution Raw Moments for the Hanson-Brennan Approach to Classification Accuracy and Consistency.
#'
#' @description An implementation of Lords (1965, p. 265) equation 37 for estimating the raw moments of the true-score distribution.
#' @param x Vector of values representing sum-scores.
#' @param r The number of raw moments to be calculated.
#' @param N The number of test items (i.e., test length).
#' @param k Lord's k (see documentation for the \code{Lords.k()} function.
#' @export
#' @examples
#' # Generate some data under the Beta Compound-Binomial distribution, where the
#' # Compound Binomial distribution has 100 trials and Lord's k = 2, and the
#' # Beta distribution has location parameters l = .15 and u = .85, and shape
#' # parameters alpha = 6 and beta = 4:
#' obs <- rBetacBinom(1000, 100, 2, .15, .85, 6, 4)
#'
#' # To estimate the first four raw moments of the underlying Beta distribution:
#' HB.tsm(x = obs, r = 4, N = 100, k = 2)
HB.tsm <- function(x, r, N, k) {
  m <- vector("numeric", r)
  for(i in 1:r) {
    if (i == 1) {
      m[i] <- mean(x) / N
      } else {
        m[i] <- 1 / (dfac(N, 2) + k*dfac(i, 2)) * ((mean(dfac(x, i)) / dfac(N - 2, i - 2)) + k*dfac(i, 2) * m[i - 1])
      }
    }
  m
}

#' Tabular organization of accuracy and consistency output from the \code{LL.CA.MC()} function.
#'
#' @description Function that takes the output from the \code{LL.CA.MC()} function and organizes it in a table with accuracy and consistency indices represented by columns and categories as rows.
#' @param x The list-output from the \code{LL.CA.MC()} function.
#' @export
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.1, 0.95, 5, 3)
#' for (i in 1:100) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 1000, ncol = 100)
#'   }
#'   rawdata[, i] <- rbinom(1000, 1, p.success)
#' }
#'
#' # Estimate accuracy and consistency where the lowest category are scores
#' # below 50, second lowest 60, then 70, 80, and 90. Using the cba() function
#' # to estimate the reliability of this test, to use the LL.CA.MC() function
#' # or estimating diagnostic performance and consistency indices of
#' # classifications when using several cut-points:
#' output <- LL.CA.MC(rowSums(rawdata), cba(rawdata), seq(50, 90, 10), 0, 100)
#'
#' # As this output can get quite verbose as the number of categories increase,
#' # the MC.out.tabular() function can be used to organize the output more
#' # concisely in a tabular format.
#' MC.out.tabular(output)
MC.out.tabular <- function(x) {
  tab.out <- base::matrix(nrow = length(x$accuracy$specific), ncol = 13)
  base::colnames(tab.out) <- base::c("TP", "FP", "TN", "FN", base::names(x$accuracy[[2]][[2]][[2]]), "p", "p_c", "Kappa")
  nams <- NULL
  for (i in 1:base::length((x$accuracy$specific))) {
    nams[i] <- base::paste("Category.", i, sep = "")
    tab.out[i, 1:4] <- base::as.vector(as.matrix(x$accuracy$specific[[i]][[1]][-3, -3]))
    for (j in 1:base::length(x$accuracy[[2]][[i]][[2]])) {
      tab.out[i, j + 4] <- x$accuracy[[2]][[i]][[2]][[j]]
    }
    for (k in 1:base::length(x$consistency[[2]][[i]][[1]])) {
      tab.out[i, k + 10] <- x$consistency[[2]][[i]][[1]][[k]]
    }
  }
  base::rownames(tab.out) <- nams
  tab.out
}

#' Graphical presentation of model fit for the Beta-Binomial classification accuracy and consistency model.
#'
#' @description Tool for visually gauging the discrepancy between the observed and model-implied frequencies of observed-scores.
#' @param x The output object from the \code{LL.CA()}, \code{LL.MC.CA()}, \code{HB.CA()}, or \code{HB.CA.MC()} functions.
#' @param x.tickat The points along the x-axis that bins are to be labeled. Default is \code{NULL} (places a tick for each of the bins).
#' @param y.tickat The points along the y-axis where frequencies are to be labelled. Default is \code{NULL}.
#' @param y.lim The limits of the y-axis (frequencies). Useful for keeping the scale equal across several plots.
#' @param main.lab The main label (title) of the plot.
#' @param x.lab The label for the x-axis (the bins).
#' @param y.lab The label for the y-axis (the frequencies).
#' @param x.grid Control the vertical grid-lines of the plot. Takes \code{NULL}, \code{NA}, or a vector of values as input. If \code{NULL}, grid-lines are drawn automatically for each bin. If \code{NA}, no grid-lines are drawn. If a vector of values are supplied, lines are drawn at each value provided along the x-axis.
#' @param y.grid Control the horizontal grid-lines of the plot. Takes \code{NULL}, \code{NA}, or a vector of values as input. If \code{NULL}, grid-lines are drawn automatically for each frequency (i.e., increments of 1). If \code{NA}, no grid-lines are drawn. If a vector of values are supplied, lines are drawn at each value provided along the y-axis.
#' @export
#' @examples
#' # Generate some data. 1000 respondents taking 100 item test:
#' set.seed(060121)
#' p.success <- rBeta.4P(1000, 0.25, 0.75, 5, 3)
#' for (i in 1:100) {
#'   if (i == 1) {
#'    rawdata <- matrix(nrow = 1000, ncol = 100)
#'  }
#'  rawdata[, i] <- rbinom(1000, 1, p.success)
#' }
#'
#' # Analyse the accuracy and consistency of the test and store the object:
#' out <- LL.CA(x = rowSums(rawdata), reliability = cba(rawdata), cut = 50,
#' min = 0, max = 100, modelfit = c(nbins = 20, minbin = 1))
#'
#' # Feed the object to the mdlfit.gfx() function:
#' mdlfit.gfx(out)
#'
#' # Given the number of observations, the y-axis ticks are a bit crowded. We
#' # can make it look less crowded by changing the number of ticks, labels, and
#' # the grid-lines:
#' mdlfit.gfx(out, y.tickat = seq(0, 250, 25), y.lim = c(0, 250),
#' y.grid = seq(0, 250, 12.5))
mdlfit.gfx <- function(x, x.tickat = NULL, y.tickat = NULL, y.lim = NULL, main.lab = "Observed vs. Expected Frequencies",  x.lab = "Bins", y.lab = "Frequency", x.grid = NULL, y.grid = NULL) {
  if (is.null(y.lim)) y.lim <- c(0, base::max(x$modelfit$contingencytable) * 1.3)
  if (is.null(x.tickat)) x.tickat <- 1:ncol(x$modelfit$contingencytable)
  if (is.null(y.tickat)) y.tickat <- seq(0, ceiling(y.lim[2]), ceiling(y.lim[2] / 10))
  plot(NULL, xlim = c(1, ncol(x$modelfit$contingencytable)), ylim = y.lim, ylab = y.lab, xlab = x.lab, axes = FALSE)
    if (is.null(x.grid)) {
      graphics::abline(v = seq(1:ncol(x$modelfit$contingencytable)), col = "lightgrey", lty = 3)
      } else {
        if (!is.na(x.grid[1])) {
          graphics::abline(v = x.grid, col = "lightgrey", lty = 3)
        }
      }
  if (is.null(y.grid)) {
    graphics::abline(h = seq(0, ceiling(y.lim[2]), ceiling(y.lim[2] / 10)), col = "lightgrey", lty = 3)
    } else {
      if (!is.na(y.grid[1])) {
        graphics::abline(h = y.grid, col = "lightgrey", lty = 3)
      }
    }
  graphics::par(new = TRUE)
    base::plot(1:ncol(x$modelfit$contingencytable), x$modelfit$contingencytable[1, ],
             type = "o", col = "grey", ylim = y.lim, ylab = "", xlab = "", lwd = 3, lty = 1,
             axes = FALSE, pch = 16, cex = 1.5)
    graphics::par(new = TRUE)
    base::plot(1:ncol(x$modelfit$contingencytable), x$modelfit$contingencytable[2, ],
             type = "o", bg = "black", ylim = y.lim, ylab = "", xlab = x.lab,
             main = main.lab, lwd = 3, axes = FALSE, pch = 1, cex = 1.5)
  graphics::box()
  graphics::axis(1, at = x.tickat, labels = FALSE)
  graphics::axis(2, at = y.tickat)
  graphics::legend("topright", legend = c("Observed",
                                          "Expected"), col = c("black", "darkgrey"),
                   lty = c(1, 1), lwd = c(3, 3), pch = c(1, 16), bty = "n", cex = 1.5)
  if (x$modelfit$p < 0.001) {
    graphics::legend("topleft",
                     legend = c(paste("chi-square = ", round(x$modelfit$chisquared, 2)),
                                paste("df = ", x$modelfit$df), "p < 0.001"),
                     bty = "n", cex = 1.5)
    } else {
      graphics::legend("topleft",
                       legend = c(paste("chi-square = ", round(x$modelfit$chisquared, 2)),
                                  paste("df = ", x$modelfit$df),
                                  paste("p = ", round(x$modelfit$p, 3))),
                       bty = "n", cex = 1.5)
    }
}



#' Function for estimating "Lord's k" for Lord's two-term approximation to the compound binomial distribution.
#'
#' @description Calculates Lord's k.
#' @param x A vector of observed-scores.
#' @param N The test length.
#' @param reliability The test-score reliability coefficient.
#' @return A value representing Lord's k
#' @export
#' @examples
#' # Generate some fictional data. Say 100 students take a 50-item long test
#' # where all items are equally difficult (i.e., where the true Lord's k = 0).
#' set.seed(1234)
#' p.success <- rBeta.4P(100, 0.25, 0.75, 5, 3)
#' for (i in 1:50) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 100, ncol = 50)
#'   }
#'   rawdata[, i] <- rbinom(100, 1, p.success)
#' }
#'
#' # Estimate the reliability of these scores with Cronbach's Alpha:
#' reliability <- cba(rawdata)
#'
#' # Estimate Lord's k using Lords.k():
#' Lords.k(rowSums(rawdata), 50, reliability)
Lords.k <- function(x, N, reliability) {
  mu <- base::mean(x)
  sigma2 <- stats::var(x)
  sigma2e <- sigma2*(1 - reliability)
  num <- N*((N-1) * (sigma2 - sigma2e) - N*sigma2 + mu*(N-mu))
  den <- 2*(mu*(N - mu) - (sigma2 - sigma2e))
  num/den
}


#' An Implementation of the Hanson and Brennan Approach to Estimate Classification Consistency and Accuracy based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Hanson and Brennan approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to a two-term approximation of the compound binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{k} for the "Lord's k" parameter (see documentation for the \code{Lords.k} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param testlength The total number of test items (or maximum possible score). Must be an integer.
#' @param cut The cutoff value for classifying observations into above/below categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves (see documentation for the \code{HB.ROC()} function).
#' @param output Character vector indicating which types of statistics (i.e, accuracy and/or consistency) are to be computed and included in the output. Permissible values are \code{"accuracy"} and \code{"consistency"}.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit by setting the minimum bin-size for expected observations. Can alternatively be set to \code{NULL} to forego model-fit testing (speeding up the function). In accordance with standard recommendations for chi-square tests the default input to this argument is 10.
#' @note This implementation of the Hanson-Brennan approach is much slower than the implementation of the Livingston and Lewis approach, as there is no native implementation of Lord's two-term approximation to the Compound-Binomial distribution in R. This implementation uses a "brute-force" method of computing the cumulative probabilities from the compound-Binomial distribution, which will by necessity be more resource intensive.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false pass/fail categorizations for a test, diagnostic performance statistics, and / or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 50.
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#'  for (i in 1:20) {
#'    if (i == 1) {
#'      rawdata <- matrix(nrow = 1000, ncol = 20)
#'      }
#'    rawdata[, i] <- rbinom(1000, 1, p.success)
#'  }
#'
#' # Suppose the cutoff value for attaining a pass is 10 items correct, and
#' # that the reliability of this test was estimated using the Cronbach's Alpha
#' # estimator. To estimate and retrieve the estimated parameters, confusion and
#' # consistency matrices, and accuracy and consistency indices using HB.CA():
#' HB.CA(x = rowSums(rawdata), reliability = cba(rawdata), cut = 10,
#' testlength = 20)
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
HB.CA <- function(x = NULL, reliability, cut, testlength, true.model = "4P", truecut = NULL, output = c("accuracy", "consistency"), failsafe = TRUE, l = 0, u = 1, modelfit = 10) {
  out <- base::list()
  if (!is.list(x)) {
    k <- Lords.k(x, testlength, reliability)
    if (startsWith(as.character(true.model), "2")) {
      failsafe <- FALSE
    }
    params <- HB.beta.tp.fit(x, testlength, k, true.model = true.model, failsafe = failsafe, l = l, u = u)
    if (params$l < 0 | params$u > 1) {
      warning(paste("Parameter out of bounds: l = ", round(params$l, 4), ", u = ", round(params$u, 4), ", alpha = ", round(params$alpha, 4), ", beta = ", round(params$beta, 4),
                    ". Consider constraining the fitting procedure further (e.g., set the location-parameters).", sep = ""))
    }
  } else {
    params <- x
  }
  if (base::is.null(truecut)) {
    truecut <- cut
  }
  truecut <- truecut / params$N
  out[["parameters"]] <- params
  if (!is.list(x) & !is.null(modelfit)) {
    mdlfit <- matrix(ncol = params$N + 1, nrow = 2)
    rownames(mdlfit) <- c("Expected", "Observed")
    for (i in 0:params$N) {
      mdlfit[1, i + 1] <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x) }, lower = 0, upper = 1)$value * length(x)
      mdlfit[2, i + 1] <- length(x[x == i])
      }
    for (i in 1:ncol(mdlfit)) {
      if (i < ncol(mdlfit)) {
        if (any(mdlfit[1, i] < ncol(mdlfit))) {
          if (any(mdlfit[1, i] < modelfit)) {
            mdlfit[, i + 1] <- mdlfit[, i + 1] + mdlfit[, i]
            mdlfit[, i] <- NA
          }
        }
      }
      }
    mdlfit <- mdlfit[, apply(mdlfit, 2, function(x) {any(!is.na(x))})]
    if (any(mdlfit[1, ncol(mdlfit)] < modelfit)) {
      mdlfit[, ncol(mdlfit) - 1] <- mdlfit[, ncol(mdlfit) - 1] + mdlfit[, ncol(mdlfit)]
      mdlfit <- mdlfit[, -ncol(mdlfit)]
      }
    chisquared <- sum(apply(mdlfit, 2, function(x) {(x[2] - x[1])^2 / x[1]}))
    out[["modelfit"]] <- list()
    out[["modelfit"]][["contingencytable"]] <- mdlfit
    out[["modelfit"]][["chisquared"]] <- chisquared
    if (startsWith(as.character(true.model), "2")) {
      out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
      } else {
        if ((startsWith(as.character(true.model), "4") & failsafe == TRUE) & (out[["parameters"]]$l == l & out[["parameters"]]$u == u)) {
          out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
          } else {
            out[["modelfit"]][["df"]] <- ncol(mdlfit) - 4
          }
        }
    out[["modelfit"]][["pvalue"]] <- stats::pchisq(chisquared, ncol(mdlfit) - 4, lower.tail = FALSE)
  }
  if (any(startsWith(tolower(output), "a"))) {
    p.fp <- NULL
    for (i in 0:(cut - 1)) {
      p.fp[i + 1] <- stats::integrate(function(x) {dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x)}, lower = truecut, upper = 1)$value
    }
    p.fp <- sum(p.fp)
    p.tp <- NULL
    for (i in 0:(cut - 1)) {
      p.tp[i + 1] <- stats::integrate(function(x) {dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x)}, lower = 0, upper = truecut)$value
    }
    p.tp <- sum(p.tp)
    p.tn <- NULL
    for (i in cut:params$N) {
      p.tn[i - cut + 1] <- stats::integrate(function(x) {dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x)}, lower = truecut, upper = 1)$value
    }
    p.tn <- sum(p.tn)
    p.fn <- NULL
    for (i in cut:params$N) {
      p.fn[i - cut + 1] <- stats::integrate(function(x) {dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x)}, lower = 0, upper = truecut)$value
    }
    p.fn <- sum(p.fn)
    camat <- confmat(p.tp, p.tn, p.fp, p.fn, "freq")
    out[["confusionmatrix"]] <- camat
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(startsWith(tolower(output), "c"))) {
    ccmat <- matrix(ncol = params$N + 1, nrow = params$N + 1)
    for (i in 0:params$N) {
      for (j in 0:params$N) {
        if (i <= j) {
          ccmat[i + 1, j + 1] <- stats::integrate(function(x) {dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x) * dcBinom(j, params$N, params$k, x)}, lower = 0, upper = 1)$value
        }
      }
    }
    ccmat[lower.tri(ccmat)] <- t(ccmat)[lower.tri(ccmat)]
    p.ii <- sum(ccmat[1:cut, 1:cut])
    p.jj <- sum(ccmat[(cut + 1):(params$N + 1), (cut + 1):(params$N + 1)])
    p.ij <- sum(ccmat[1:cut, (cut + 1):(params$N + 1)])
    p.ji <- sum(ccmat[(cut + 1):(params$N + 1), 1:cut])
    ccmat <- base::matrix(nrow = 2, ncol = 2, dimnames = list(c("i", "j"), c("i", "j")))
    ccmat["i", "i"] <- p.ii
    ccmat["i", "j"] <- p.ij
    ccmat["j", "i"] <- p.ji
    ccmat["j", "j"] <- p.jj
    ccmat <- ccmat / sum(ccmat)
    out[["consistencymatrix"]] <- ccmat
    out[["classification.consistency"]] <- ccStats(ccmat["i", "i"], ccmat["i", "j"], ccmat["j", "i"], ccmat["j", "j"])
  }
  base::return(out)
}
#' ROC curves for the Hanson and Brennan approach.
#'
#' @description Generate a ROC curve plotting the false-positive rate against the true-positive rate at different cut-off values across the observed-score scale.
#' @param x A vector of observed results (sum scores) or a list of parameter values (see documentation for the \code{HB.beta.tp.fit() function}.
#' @param reliability The reliability coefficient of the test.
#' @param testlength The total number of test items (or maximum possible score). Must be an integer.
#' @param truecut The point along the x-scale that marks true category membership.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The \code{"4P"} method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the \code{"2P"} method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe If true-model == "4P": Whether to engage a fail-safe reverting to a two-parameter true-score distribution solution should the four-parameter fitting procedure produce impermissible results. Default is TRUE (engage fail-safe in the event of impermissible estimates).
#' @param l If \code{true.model == "2P"} or \code{failsafe == TRUE}: The lower-bound location parameter of the two-parameter true-score distribution solution.
#' @param u If \code{true.model == "2P"} or \code{failsafe == TRUE}: The upper-bound location parameter of the two-parameter true-score distribution solution.
#' @param AUC Logical. Calculate and include the area under the curve? Default is \code{FALSE}.
#' @param maxJ Logical. Mark the point along the curve where Youden's J statistic is maximized? Default is \code{FALSE}.
#' @param maxAcc Logical. Mark the point along the curve where the Accuracy statistic is maximized? Default is \code{FALSE}.
#' @param locate Ask the function to locate the cut-point at which sensitivity or NPV is greater than or equal to some value, or specificity or PPV is lesser than or equal to some value. Take as input a character-vector of length 2, with the first argument being which index is to be found (e.g., "sensitivity"), and the second argument the value to locate (e.g., "0.75"). For example: c("sensitivity", "0.75").
#' @param raw.out Give raw coordinates as output rather than plot? Default is \code{FALSE}
#' @param grainsize Specify the number of cutoff-points for which the ROC curve is to be calculated. The greater this number the greater the accuracy. Default is set to the stated test length (N).
#' @note This implementation of the Hanson-Brennan approach is much slower than the implementation of the Livingston and Lewis approach, as there is no native implementation of Lord's two-term approximation to the Compound-Binomial distribution in R. This implementation uses a "brute-force" method of computing the cumulative probabilities from the compound-Binomial distribution, which will by necessity be more resource intensive.
#' @return A plot tracing the ROC curve for the test, or matrix of coordinates if raw.out is \code{TRUE}.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a test with a
#' # maximum score of 50.
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#'  for (i in 1:20) {
#'    if (i == 1) {
#'      rawdata <- matrix(nrow = 1000, ncol = 20)
#'      }
#'    rawdata[, i] <- rbinom(1000, 1, p.success)
#'  }
#'
#' # Suppose the cutoff value for attaining a pass is 10 items correct, and
#' # that the reliability of this test was estimated using the Cronbach's Alpha
#' # estimator. To draw the ROC-graph and locate the points at which Youden's J
#' # and Accuracy are maximized:
#' HB.ROC(rowSums(rawdata), cba(rawdata), 20, 10, maxAcc = TRUE, maxJ = TRUE)
#'
#' # For further examples regarding how to use the locate argument to locate
#' # points at which various criteria are satisfied, see documentation for the
#' # LL.ROC() function.
#' @export
HB.ROC <- function(x = NULL, reliability, testlength, truecut, true.model = "4P", failsafe = TRUE, l = 0, u = 1, AUC = FALSE, maxJ = FALSE, maxAcc = FALSE, locate = NULL, raw.out = FALSE, grainsize = testlength) {
  if (!raw.out) {
    oldpar <- graphics::par(no.readonly = TRUE)
    base::on.exit(graphics::par(oldpar))
  }
  if (!is.list(x)) {
    k <- Lords.k(x, testlength, reliability)
    x <- HB.beta.tp.fit(x, testlength, k, true.model = true.model, failsafe = failsafe, l = l, u = u)
  }
  for (i in 1:(grainsize + 1)) {
    if (i == 1) {
      cuts <- seq(0, testlength, testlength / grainsize)
      outputmatrix <- matrix(nrow = grainsize + 1, ncol = 7)
      outputmatrix[, 4] <- cuts
    }
    axval <- HB.CA(x = x, cut = cuts[i],
                   truecut = truecut, true.model = true.model,
                   output = "a", l = l, u = u, modelfit = NULL)$classification.accuracy
    outputmatrix[i, 1] <- 1 - axval$Specificity
    outputmatrix[i, 2] <- axval$Sensitivity
    outputmatrix[i, 3] <- axval$Youden.J
    outputmatrix[i, 5] <- axval$Accuracy
    outputmatrix[i, 6] <- axval$PPV
    outputmatrix[i, 7] <- axval$NPV
    base::colnames(outputmatrix) <- c("FPR", "TPR", "Youden.J", "Cutoff", "Accuracy", "PPV", "NPV")
    outputmatrix[base::which(base::is.na(outputmatrix[, 1])), 1] <- 0
    outputmatrix[base::which(base::is.na(outputmatrix[, 2])), 2] <- 1
    outputmatrix[base::which(base::is.na(outputmatrix[, 6])), 6] <- 1
    outputmatrix[base::which(base::is.na(outputmatrix[, 7])), 7] <- 1
  }
  if (raw.out) {
    base::return(outputmatrix)
  }
  graphics::plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
  graphics::abline(h = seq(0, 1, .1), v = seq(0, 1, .1), col = "lightgrey", lty = "dotted")
  graphics::par(new = TRUE)
  graphics::plot(outputmatrix[, 1], outputmatrix[, 2], type = "s",
                 xlab = "False-Positive Rate (1 - Specificity)",
                 ylab = "True-Positive Rate (Sensitivity)",
                 main = paste("ROC curve for true-cut equal to", truecut), lwd = 2,
                 xlim = c(0, 1), ylim = c(0, 1))
  if (AUC) {
    graphics::legend("bottomright", bty = "n", cex = 1.5,
                     legend = paste("AUC =", round(AUC(outputmatrix[, 1], outputmatrix[, 2]), 3)))
  }
  if (maxJ) {
    graphics::points(outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 1],
                     outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 2], cex = 1.5, pch = 19)
    graphics::text(outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 1] + .025,
                   outputmatrix[base::which(outputmatrix[, 3] == base::max(outputmatrix[, 3])), 2] - .025,
                   labels = base::paste("Maximum Youden's J. ", "(", base::round(base::max(outputmatrix[, 3]), 3), ") ",  "at cut-off = ",
                                        base::round(outputmatrix[which(outputmatrix[, 3] == base::max(outputmatrix[, 3]))[1], 4], 3), ".", sep = ""),
                   adj = c(0, 1))
  }
  if (maxAcc) {
    graphics::points(outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 1],
                     outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 2], cex = 1.5, pch = 19)
    graphics::text(outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 1] + .025,
                   outputmatrix[base::which(outputmatrix[, 5] == base::max(outputmatrix[, 5])), 2] - .025,
                   labels = base::paste("Maximum Accuracy ", "(", base::round(base::max(outputmatrix[, 5]), 3), ") ",  "at cut-off = ",
                                        base::round(outputmatrix[which(outputmatrix[, 5] == base::max(outputmatrix[, 5]))[1], 4], 3), ".", sep = ""),
                   adj = c(0, 1))
  }
  if (!base::is.null(locate)) {
    if (base::startsWith(base::tolower(locate[1]), "se")) {
      rowloc <- base::which(outputmatrix[, "TPR"] >= base::as.numeric(locate[2]))[1]
      colloc <- "TPR"
      value <- outputmatrix[rowloc, "TPR"]
      statistic <- "Sensitivity >= "
    }
    if (base::startsWith(base::tolower(locate[1]), "p")) {
      rowloc <- base::which(outputmatrix[, "PPV"] <= base::as.numeric(locate[2]))[1]
      colloc <- "PPV"
      value <- outputmatrix[rowloc, "PPV"]
      statistic <- "PPV <= "
    }
    if (base::startsWith(base::tolower(locate[1]), "n")) {
      rowloc <- base::which(outputmatrix[, "NPV"] >= base::as.numeric(locate[2]))[1]
      colloc <- "NPV"
      value <- outputmatrix[rowloc, "NPV"]
      statistic <- "NPV >= "
    }
    if (base::startsWith(base::tolower(locate[1]), "sp")) {
      rowloc <- base::which(outputmatrix[, "FPR"] <= 1 - base::as.numeric(locate[2]))[base::length(base::which(outputmatrix[, "FPR"] <= 1 - base::as.numeric(locate[2])))]
      colloc <- "FPR"
      value <- 1 - outputmatrix[rowloc, "FPR"]
      statistic <- "Specificity <= "
    }
    graphics::points(outputmatrix[rowloc, 1], outputmatrix[rowloc, 2], cex = 1.5, pch = 19)
    graphics::text(outputmatrix[rowloc, 1] + .025,
                   outputmatrix[rowloc, 2] - .025,
                   labels = base::paste(statistic, base::as.numeric(locate[2]) , " (", base::round(value, 3), ")",
                                        " at cut-off = ", outputmatrix[rowloc, 4], ".", sep = ""),
                   adj = if (base::startsWith(base::tolower(locate[1]), "sp")) { c(0, 1) } else { c(0, 1) })
  }
}



#' An Extension of the Hanson and Brennan Approach to Estimate Classification Consistency and Accuracy for Multiple Classifications based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Hanson and Brennan approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores, or a list specifying parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location-, and \code{alpha} and \code{beta} for the shape parameters of the Beta true-score distribution, and \code{k} for the "Lord's k" parameter (see documentation for the \code{Lords.k} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param testlength The total number of test items (or maximum possible score). Must be an integer.
#' @param cut A vector of cut-off values for classifying observations into two or more categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param failsafe Logical value indicating whether to engage the automatic fail-safe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{TRUE} (i.e., the function will engage failsafe if the four-parameter Beta-distribution fitting-procedure produced impermissible estimates).
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param modelfit Allows for controlling the chi-square test for model fit by setting the minimum bin-size for expected observations. Can alternatively be set to \code{NULL} to forego model-fit testing (speeding up the function). In accordance with standard recommendations for chi-square tests the default input to this argument is 10.
#' @return A list containing the estimated parameters necessary for the approach (i.e., Lord's k, test-length, and the true-score Beta distribution parameters), a chi-square test of model-fit, the confusion matrix containing estimated proportions of true/false positive/negative categorizations for a test, diagnostic performance statistics, and/or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note This implementation of the Hanson-Brennan approach is much slower than the implementation of the Livingston and Lewis approach, as there is no native implementation of Lord's two-term approximation to the Compound-Binomial distribution in R. This implementation uses a "brute-force" method of computing the cumulative probabilities from the compound-Binomial distribution, which will by necessity be more resource intensive.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#'  for (i in 1:20) {
#'    if (i == 1) {
#'      rawdata <- matrix(nrow = 1000, ncol = 20)
#'      }
#'    rawdata[, i] <- rbinom(1000, 1, p.success)
#'  }
#'
#' # Suppose the cutoff value for attaining a pass is 10 items correct, and
#' # that the reliability of this test was estimated using the Cronbach's Alpha
#' # estimator. To estimate and retrieve the estimated parameters, confusion and
#' # consistency matrices, and accuracy and consistency indices using HB.CA():
#' (output <- HB.CA.MC(x = rowSums(rawdata), reliability = cba(rawdata),
#' cut = c(8, 12), testlength = 20))
#'
#' # The output for this function can get quite verbose as more categories are
#' # included. The output from the function can be fed to the MC.out.tabular()
#' # function in order to organize the output in a tabular format.
#' MC.out.tabular(output)
#'
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @references Lord. Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @references Lewis, Don and Burke, C. J. (1949). The Use and Misuse of the Chi-Square Test. Psychological Bulletin, 46(6).
#' @export
HB.CA.MC <- function(x = NULL, reliability, cut, testlength, true.model = "4P", failsafe = TRUE, l = 0, u = 1, modelfit = 10) {
  out <- base::list()
  if (!is.list(x)) {
    k <- Lords.k(x, testlength, reliability)
    if (startsWith(as.character(true.model), "2")) {
      failsafe <- FALSE
    }
    params <- HB.beta.tp.fit(x, testlength, k, true.model = true.model, failsafe = failsafe, l = l, u = u)
    if (params$l < 0 | params$u > 1) {
      warning(paste("Parameter out of bounds: l = ", round(params$l, 4), ", u = ", round(params$u, 4), ", alpha = ", round(params$alpha, 4), ", beta = ", round(params$beta, 4),
                    ". Consider constraining the fitting procedure further (e.g., set the location-parameters).", sep = ""))
    }
  } else {
    params <- x
  }
  out[["parameters"]] <- params
  truecut <- c(0, cut / params$N, 1)
  cut <- c(0, cut, params$N)
  camat <- matrix(ncol = length(cut) -1, nrow = length(cut) - 1)
  ccmat <- camat
  if (!is.list(x) & !is.null(modelfit)) {
    mdlfit <- matrix(ncol = params$N + 1, nrow = 2)
    rownames(mdlfit) <- c("Expected", "Observed")
    for (i in 0:params$N) {
      mdlfit[1, i + 1] <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x) }, lower = 0, upper = 1)$value * length(x)
      mdlfit[2, i + 1] <- length(x[x == i])
    }
    for (i in 1:ncol(mdlfit)) {
      if (i < ncol(mdlfit)) {
        if (any(mdlfit[1, i] < ncol(mdlfit))) {
          if (any(mdlfit[1, i] < modelfit)) {
            mdlfit[, i + 1] <- mdlfit[, i + 1] + mdlfit[, i]
            mdlfit[, i] <- NA
          }
        }
      }
    }
    mdlfit <- mdlfit[, apply(mdlfit, 2, function(x) {any(!is.na(x))})]
    if (any(mdlfit[1, ncol(mdlfit)] < modelfit)) {
      mdlfit[, ncol(mdlfit) - 1] <- mdlfit[, ncol(mdlfit) - 1] + mdlfit[, ncol(mdlfit)]
      mdlfit <- mdlfit[, -ncol(mdlfit)]
    }
    chisquared <- sum(apply(mdlfit, 2, function(x) {(x[2] - x[1])^2 / x[1]}))
    out[["modelfit"]] <- list()
    out[["modelfit"]][["contingencytable"]] <- mdlfit
    out[["modelfit"]][["chisquared"]] <- chisquared
    if (startsWith(as.character(true.model), "2")) {
      out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
    } else {
      if ((startsWith(as.character(true.model), "4") & failsafe == TRUE) & (out[["parameters"]]$l == l & out[["parameters"]]$u == u)) {
        out[["modelfit"]][["df"]] <- ncol(mdlfit) - 2
      } else {
        out[["modelfit"]][["df"]] <- ncol(mdlfit) - 4
      }
    }
    out[["modelfit"]][["pvalue"]] <- stats::pchisq(chisquared, ncol(mdlfit) - 4, lower.tail = FALSE)
  }
  for (i in 1:(length(cut) - 1)) {
    if (i == 1) {
      for (j in 1:(length(cut) - 1)) {
        if (j == 1) {
          rnam <- NULL
          cnam <- NULL
        }
        if (j != (length(cut) - 1)) {
          if (j == 1) {
            rnam[j] <- paste("Observed <", cut[j + 1])
            cnam[j] <- paste("True     <", cut[j + 1])
          } else {
            rnam[j] <- paste(" >=", cut[j], "& <", cut[j + 1])
            cnam[j] <- paste(" >=", cut[j], "& <", cut[j + 1])
          }
        } else {
          rnam[j] <- paste(" >=", cut[j])
          cnam[j] <- paste(" >=", cut[j])
        }
      }
      colnames(camat) <- cnam
      rownames(camat) <- rnam
    }
  }

  camat.bc <- matrix(ncol = (length(cut) - 1), nrow = params$N + 1)
  for(i in 0:(nrow(camat.bc) - 1)) {
    for (j in 1:ncol(camat.bc)) {
      camat.bc[i + 1, j] <- stats::integrate(function(x) {
        dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x)
      }, lower = truecut[j], upper = truecut[j + 1])$value
    }
  }
  for (i in 1:ncol(camat)) {
    for (j in 1:ncol(camat)) {
      if (i != ncol(camat)) {
        camat[i, j] <- sum(camat.bc[(cut[i] + 1):(cut[i + 1]), j])
      } else {
        camat[i, j] <- sum(camat.bc[(cut[i] + 1):(cut[i + 1] + 1), j])
      }
    }
  }
  out[["accuracy"]] <- list()
  out[["accuracy"]][["overall"]] <- list()
  out[["accuracy"]][["overall"]][["confusionmatrix"]] <- camat / sum(camat)
  out[["accuracy"]][["overall"]][["accuracy"]] <- sum(diag(camat))
  out[["accuracy"]][["specific"]]
  caout <- list()
  for(i in 1:ncol(camat)) {
    FN <- sum(camat[-i, i])
    FP <- sum(camat[i, -i])
    TP <- camat[i, i]
    TN <- sum(camat[-i, -i])
    caout[[paste("Category.", i, sep = "")]] <- list()
    caout[[i]][["confusionmatrix"]] <- confmat(TP, TN, FP, FN)
    caout[[i]][["statistics"]] <- caStats(TP, TN, FP, FN)
  }
  out[["accuracy"]][["specific"]] <- caout
  for (i in 1:(length(cut) - 1)) {
    if (i == 1) {
      for (j in 1:(length(cut) - 1)) {
        if (j == 1) {
          rnam <- NULL
          cnam <- NULL
        }
        if (j != (length(cut) - 1)) {
          if (j == 1) {
            rnam[j] <- paste("2nd adm. <", cut[j + 1])
            cnam[j] <- paste("1st adm. <", cut[j + 1])
          } else {
            rnam[j] <- paste(" >=", cut[j], "& <", cut[j + 1])
            cnam[j] <- paste(" >=", cut[j], "& <", cut[j + 1])
          }
        } else {
          rnam[j] <- paste(" >=", cut[j])
          cnam[j] <- paste(" >=", cut[j])
        }
      }
      colnames(ccmat) <- cnam
      rownames(ccmat) <- rnam
    }
    ccmat.bc <- matrix(ncol = params$N + 1, nrow = params$N + 1)
    for (i in 0:params$N) {
      for (j in 0:params$N) {
        if (i <= j) {
          ccmat.bc[i + 1, j + 1] <- stats::integrate(function(x) {dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * dcBinom(i, params$N, params$k, x) * dcBinom(j, params$N, params$k, x)}, lower = 0, upper = 1)$value
        }
      }
    }
    ccmat.bc[lower.tri(ccmat.bc)] <- t(ccmat.bc)[lower.tri(ccmat.bc)]
    for (i in 1:ncol(ccmat)) {
      for (j in 1:ncol(ccmat)) {
        if (i == 1 & j == 1) {
          ccmat[i, j] <- sum(ccmat.bc[1:(cut[i + 1]), 1:(cut[j + 1])])
        }
        if (i == 1 & (j != 1 & j != ncol(ccmat))) {
          ccmat[i, j] <- sum(ccmat.bc[1:(cut[i + 1]), (cut[j] + 1):cut[j + 1]])
        }
        if (i == 1 & (j == ncol(ccmat))) {
          ccmat[i, j] <- sum(ccmat.bc[1:(cut[i + 1]), (cut[j] + 1):(cut[j + 1] + 1)])
        }
        if ((i != 1 & i != ncol(ccmat)) & j == 1) {
          ccmat[i, j] <- sum(ccmat.bc[(cut[i] + 1):cut[i + 1], 1:(cut[j + 1])])
        }
        if ((i != 1 & i != ncol(ccmat)) & (j != 1 & j != ncol(ccmat))) {
          ccmat[i, j] <- sum(ccmat.bc[(cut[i] + 1):cut[i + 1], (cut[j] + 1):cut[j + 1]])
        }
        if ((i != 1 & i != ncol(ccmat)) & j == ncol(ccmat)) {
          ccmat[i, j] <- sum(ccmat.bc[(cut[i] + 1):cut[i + 1], (cut[j] + 1):(cut[j + 1] + 1)])
        }
        if (i == ncol(ccmat) & j == 1) {
          ccmat[i, j] <- sum(ccmat.bc[(cut[i] + 1):(cut[i + 1] + 1), 1:(cut[j + 1])])
        }
        if (i == ncol(ccmat) & (j != 1 & j != ncol(ccmat))) {
          ccmat[i, j] <- sum(ccmat.bc[(cut[i] + 1):(cut[i + 1] + 1), (cut[j] + 1):cut[j + 1]])
        }
        if (i == ncol(ccmat) & j == ncol(ccmat)) {
          ccmat[i, j] <- sum(ccmat.bc[(cut[i] + 1):(cut[i + 1] + 1), (cut[j] + 1):(cut[j + 1] + 1)])
        }
      }
    }
  }
  out[["consistency"]] <- list()
  out[["consistency"]][["overall"]] <- list()
  out[["consistency"]][["overall"]][["consistencymatrix"]] <- ccmat / sum(ccmat)
  p <- sum(diag(ccmat))
  p_c <- sum(apply(ccmat, 2, function(x) {
    sum(x)^2
  }))
  Kappa <- (p - p_c) / (1 - p_c)
  out[["consistency"]][["overall"]][["statistics"]] <- list("p" = p, "p_c" = p_c, "Kappa" = Kappa)
  ccout <- list()
  for(i in 1:ncol(ccmat)) {
    p <- ccmat[i, i]
    p_c <- sum(ccmat[i, ])^2
    Kappa <- (p - p_c) / (1 - p_c)
    ccout[[paste("Category.", i, sep = "")]] <- list()
    ccout[[i]][["statistics"]] <- list("p" = p, "p_c" = p_c, "Kappa" = Kappa)
  }
  out[["consistency"]][["specific"]] <- ccout
  base::return(out)
}
