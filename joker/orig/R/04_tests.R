#' @title Estimation and Variance Tests
#' @name tests
#'
#' @description
#' This set of functions employs Monte Carlo simulations to check the
#' consistency of the estimators (i.e. that the estimators are coded correctly)
#' and their asymptotic normality (i.e. that their asymptotic variance is
#' coded correctly).
#'
#' @srrstats {G5.6, G5.6a, G5.6b, G5.7} Parameter recovery tests are extensively
#' used in the package tests, employing Monte Carlo simulations to check the
#' consistency and asymptotic normality of the estimators.
#'
#' @param est character. The estimator to be tested.
#' @param D0 An object of class `Distribution`.
#' @param n integer. The sample size to be simulated.
#' @param m integer. The number of samples to be simulated.
#' @param seed integer. Passed to `set.seed()`.
#' @param bar logical. Should a progress bar be printed?
#' @param ... extra arguments passed to the estimator.
#'
#' @returns A list with the simulation and the expected results so that they
#' can be compared in tests.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' D <- Beta(2, 3)
#' test1 <- test_consistency("me", D)
#' test2 <- test_avar("mle", D)
#' }
test_consistency <- function(est, D0, n = 1e4, seed = 1, ...) {

  # Random Sampling
  set.seed(seed)
  sam <- r(D0)(n)

  # Return
  list(prm_true = get_unknown_params(D0),
       prm_est = do.call(est, list(distr = D0, x = sam, ...)))

}

#' @rdname tests
test_avar <- function(est, D0, n = 1e4, m = 1e3, seed = 1, bar = FALSE, ...) {

  # Preliminaries
  set.seed(seed)
  params <- get_unknown_params(D0, list = FALSE)
  y <- matrix(nrow = m, ncol = length(params))

  # Loading bar
  if (bar) {
    start <- Sys.time()
  }

  # Estimation
  for (i in 1:m) {

    # Progress Bar
    if (bar) {
      progress_bar(i, m, start = start, message = "Computing:")
    }

    # CLT
    sam <- r(D0)(n)
    params_est <- unlist(do.call(est, list(x = sam, distr = D0, ...)))
    y[i, ] <- sqrt(n) * (params - params_est)

  }

  # Calculate avar
  avar_est <- var(y)
  prm_names <- names(get_unknown_params(D0, list = FALSE))

  if (nrow(avar_est) == 1) {
    avar_est <- as.vector(avar_est)
    names(avar_est) <- prm_names
  } else {
    dimnames(avar_est) <- list(prm_names, prm_names)
  }

  # Return
  list(avar_true = do.call(paste0("avar_", est), list(distr = D0)),
       avar_est = avar_est)

}
