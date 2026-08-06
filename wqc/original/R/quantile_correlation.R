#' Apply Quantile Correlation Analysis
#'
#' @importFrom waveslim  mra
#' @importFrom QCSIS    qc
#' @importFrom stats    rnorm sd quantile
#'
#' @param data Data frame containing the time series data. The first column is the reference series; subsequent columns are the target series.
#' @param quantiles Numeric vector of quantiles.
#' @param wf Wavelet family name.
#' @param J Decomposition level.
#' @param n_sim Number of simulations for confidence intervals.
#' @return A combined data.frame of quantile correlation results, with one row per level-quantile-series combination.
#' @examples
#' data <- data.frame(x = rnorm(1000), y = rnorm(1000), z = rnorm(1000))
#' quantiles <- c(0.05, 0.5, 0.95)
#' res_df <- apply_quantile_correlation(data, quantiles,n_sim=10)
#' head(res_df)
#' @export
apply_quantile_correlation <- function(data, quantiles, wf = 'la8', J = 8, n_sim = 1000) {
  x <- data[[1]]
  results_list <- lapply(colnames(data)[-1], function(col_name) {
    y <- data[[col_name]]
    res <- quantile_correlation_analysis(x, y, quantiles, wf, J, n_sim)
    cbind(Series = col_name, res)
  })
  do.call(rbind, results_list)
}

#' Quantile Correlation Analysis
#'
#' @importFrom waveslim  mra
#' @importFrom QCSIS    qc
#' @importFrom stats    rnorm sd quantile
#'
#' @param x Numeric vector for the first time series.
#' @param y Numeric vector for the second time series.
#' @param quantiles Numeric vector of quantiles.
#' @param wf Wavelet family name.
#' @param J Decomposition level.
#' @param n_sim Number of simulations for confidence intervals.
#' @return Data frame with quantile correlation estimates and confidence intervals for one pair of series.
#' @examples
#' data <- data.frame(x = rnorm(1000), y = rnorm(1000))
#' quantiles <- c(0.05, 0.5, 0.95)
#' result <- quantile_correlation_analysis(data$x, data$y, quantiles,n_sim=10)
#' head(result)
#' @export
quantile_correlation_analysis <- function(x, y, quantiles, wf = 'la8', J = 8, n_sim = 1000) {
  decomp_x <- mra(x, wf = wf, J = J, method = 'modwt', boundary = 'periodic')
  decomp_y <- mra(y, wf = wf, J = J, method = 'modwt', boundary = 'periodic')

  details_x <- decomp_x[-(J + 1)]
  details_y <- decomp_y[-(J + 1)]

  qc_results <- lapply(seq_len(J), function(j) {
    sapply(quantiles, function(q) as.numeric(qc(details_x[[j]], details_y[[j]], tau = q)[2]))
  })

  n <- length(x)
  simulated_qc <- array(NA, dim = c(n_sim, J, length(quantiles)))

  for (j in seq_len(J)) {
    for (sim in seq_len(n_sim)) {
      sim_x <- rnorm(n, mean = mean(details_x[[j]]), sd = sd(details_x[[j]]))
      sim_y <- rnorm(n, mean = mean(details_y[[j]]), sd = sd(details_y[[j]]))
      qc_vals <- sapply(quantiles, function(q) as.numeric(qc(sim_x, sim_y, tau = q)[2]))
      simulated_qc[sim, j, ] <- qc_vals
    }
  }

  ci_lower <- apply(simulated_qc, c(2, 3), quantile, probs = 0.025)
  ci_upper <- apply(simulated_qc, c(2, 3), quantile, probs = 0.975)

  data.frame(
    Level        = rep(seq_len(J), each = length(quantiles)),
    Quantile     = rep(quantiles, times = J),
    Estimated_QC = unlist(qc_results),
    CI_Lower     = as.vector(ci_lower),
    CI_Upper     = as.vector(ci_upper),
    stringsAsFactors = FALSE
  )
}
