# Declare global variables for ggplot2 (avoids NOTEs in R CMD check)
utils::globalVariables(c("t", "beta_tmg", "sector"))

#' Plot beta(TMG_t) evolution by sector
#'
#' Creates a line plot showing the evolution of time-varying beta
#' coefficients for selected sectors.
#'
#' @param results List returned by \code{\link{fit_ou_nonlinear_tmg}}
#' @param sectors Character or integer vector. Sectors to plot.
#'   If NULL, plots all sectors.
#'
#' @return A ggplot2 object if ggplot2 is available, otherwise NULL
#'   with a base R plot produced as side effect.
#'
#' @examples
#' \donttest{
#' # 1. Create mock data (T x S matrix)
#' T_obs <- 50
#' S <- 3
#' beta_mat <- matrix(rnorm(T_obs * S), nrow = T_obs, ncol = S)
#' colnames(beta_mat) <- paste0("Sector_", 1:S)
#' 
#' # 2. Wrap in list structure expected by function
#' results_mock <- list(
#'   beta_tmg = list(
#'     beta_point = beta_mat
#'   )
#' )
#' 
#' # 3. Plot
#' plot_beta_tmg(results_mock)
#' }
#'
#' @export
plot_beta_tmg <- function(results, sectors = NULL) {
  
  bt <- results$beta_tmg$beta_point
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    
    if (!is.null(sectors)) {
      if (is.character(sectors)) {
        bt <- bt[, sectors, drop = FALSE]
      } else {
        bt <- bt[, sectors, drop = FALSE]
      }
    }
    
    graphics::matplot(
      bt,
      type = "l",
      lwd = 1.5,
      lty = 1,
      xlab = "Time",
      ylab = expression(beta[TMG](t)),
      main = "Beta(TMG_t) evolution by sector"
    )
    graphics::grid()
    return(invisible(NULL))
  }
  
  df <- as.data.frame(bt)
  df$t <- seq_len(nrow(df))
  dfl <- tidyr::pivot_longer(df, -t, names_to = "sector", values_to = "beta_tmg")
  
  if (!is.null(sectors)) {
    dfl <- dfl[dfl$sector %in% sectors, ]
  }
  
  ggplot2::ggplot(dfl, ggplot2::aes(t, beta_tmg, color = sector)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Time",
      y = expression(beta[TMG](t)),
      title = "Beta(TMG_t) evolution by sector"
    ) +
    ggplot2::theme_minimal()
}


#' Plot cubic OU drift curves
#'
#' Displays the drift function for selected sectors showing mean reversion
#' with cubic nonlinearity.
#'
#' @param results List returned by \code{\link{fit_ou_nonlinear_tmg}}
#' @param sectors Integer vector. Sector indices to plot.
#'   If NULL, plots all sectors.
#'
#' @return NULL invisibly. Produces a base R plot as side effect.
#'
#' @examples
#' \donttest{
#' # 1. Create mock data
#' # z: vector of state deviations
#' z_seq <- seq(-3, 3, length.out = 100)
#' # drift: matrix (rows=z, cols=sectors)
#' drift_mat <- cbind(
#'   -0.5 * z_seq - 0.1 * z_seq^3, # Sector 1
#'   -0.8 * z_seq - 0.05 * z_seq^3 # Sector 2
#' )
#' 
#' # 2. Wrap in list structure
#' results_mock <- list(
#'   nonlinear = list(
#'     drift_decomp = list(
#'       z = z_seq,
#'       drift = drift_mat
#'     )
#'   )
#' )
#' 
#' # 3. Plot
#' plot_drift_curves(results_mock)
#' }
#'
#' @export
plot_drift_curves <- function(results, sectors = NULL) {
  dc <- results$nonlinear$drift_decomp
  z <- dc$z
  M <- dc$drift
  
  if (is.null(sectors)) {
    sectors <- seq_len(ncol(M))
  }
  
  graphics::matplot(
    z,
    M[, sectors, drop = FALSE],
    type = "l",
    lwd = 2,
    lty = 1,
    xlab = "z = Y - theta",
    ylab = "E[dY|z]",
    main = "Cubic OU drift function"
  )
  graphics::abline(h = 0, lty = 2, col = "gray50")
  graphics::grid()
  
  invisible(NULL)
}


#' Plot stochastic volatility evolution
#'
#' Displays the estimated volatility path for a selected sector.
#'
#' @param results List returned by \code{\link{fit_ou_nonlinear_tmg}}
#' @param sector Integer. Sector index to plot. Default 1.
#'
#' @return NULL invisibly. Produces a base R plot as side effect.
#'
#' @examples
#' \donttest{
#' # 1. Create mock data (Volatility must be positive)
#' T_obs <- 50
#' sigma_mat <- matrix(exp(rnorm(T_obs * 2)), nrow = T_obs, ncol = 2)
#' 
#' # 2. Wrap in list structure
#' results_mock <- list(
#'   sv = list(
#'     h_summary = list(
#'       sigma_t = sigma_mat
#'     )
#'   )
#' )
#' 
#' # 3. Plot sector 1
#' plot_sv_evolution(results_mock, sector = 1)
#' }
#'
#' @export
plot_sv_evolution <- function(results, sector = 1) {
  sig <- results$sv$h_summary$sigma_t[, sector]
  
  graphics::plot(
    sig,
    type = "l",
    lwd = 2,
    main = sprintf("Stochastic volatility (sector %d)", sector),
    xlab = "t",
    ylab = expression(sigma[t])
  )
  graphics::grid()
  
  invisible(NULL)
}


#' Plot posterior distributions of key parameters
#'
#' Creates density plots for selected parameters.
#'
#' @param fit Fitted Stan model object
#' @param params Character vector. Parameter names to plot.
#' @param verbose Logical. Print progress messages. Default FALSE.
#'
#' @return NULL invisibly. Produces plots as side effect.
#'
#' @keywords internal
plot_posterior_densities <- function(fit, params = c("beta1", "nu"), 
                                     verbose = FALSE) {
  
  for (p in params) {
    vmsg(sprintf("Plotting density for %s", p), verbose)
    
    if (inherits(fit, "CmdStanMCMC")) {
      draws <- as.vector(fit$draws(p, format = "matrix"))
    } else {
      draws <- as.vector(rstan::extract(fit, p)[[1]])
    }
    
    graphics::hist(
      draws,
      breaks = 50,
      freq = FALSE,
      main = sprintf("Posterior density: %s", p),
      xlab = p,
      col = "lightblue",
      border = "white"
    )
    graphics::lines(stats::density(draws), lwd = 2, col = "darkblue")
  }
  
  invisible(NULL)
}