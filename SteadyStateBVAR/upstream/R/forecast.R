#' Forecast from a fitted steady-state BVAR model
#'
#' Computes and plots forecasts from a fitted steady-state \code{bvar} object. Draws from the 
#' joint predictive distribution are used to construct point forecasts and
#' prediction intervals. Optionally converts monthly or quarterly growth rate forecasts to annual growth
#' rates for variables specified as \code{100*diff(log(x))}. Optionally overlays the posterior steady
#' states.
#'
#' @param x A steady-state \code{bvar} object that has been passed through \code{\link{fit}}.
#' @param pi Numeric. The prediction interval width. Default \code{0.95}, i.e. 95% prediction interval.
#' @param fcst_type Character. Whether to use \code{"mean"} or \code{"median"}
#'   as the point forecast. Default \code{"mean"}.
#' @param growth_rate_idx Integer vector. Indices of variables of which to convert forecasts to
#'   annual growth rates \eqn{100 (\ln x_{t} - \ln x_{t-f})}, where \eqn{f} is
#'   the frequency of the data (4 for quarterly, 12 for monthly).
#'   Only suitable for variables specified as \eqn{100 (\ln x_{t} - \ln x_{t-1})}, i.e. \code{100*diff(log(x))}.
#'   Computed by summing up to \eqn{f} log first differences. Default is \code{NULL}.
#' @param plot_idx Integer vector. Indices of variables to plot. If \code{NULL}
#'   (default), all variables are plotted. Forecasts are always computed and
#'   returned for all variables, regardless of \code{plot_idx}.
#' @param show_all Logical. If \code{TRUE} (default), the full
#'   history is shown. If \code{FALSE}, the last eight years
#'   of history are shown.
#' @param ss Logical. If \code{TRUE}, overlays the posterior steady state
#'   \eqn{\mu_t = \Psi d_t}. Default \code{TRUE}. For variables in \code{growth_rate_idx},
#'   the posterior steady-state is annualized.
#' @param ss_type Character. Whether to use \code{"mean"} or \code{"median"}
#' of the posterior as the steady-state point estimate. Default \code{"mean"}.
#' @param ss_ci Numeric. The posterior credible interval width for the steady state.
#' Default \code{0.95}, i.e. 95% credible interval.
#'
#' @return Invisibly returns a list with three matrices: \code{forecast}, \code{lower}, and
#'   \code{upper}, each of dimension \code{H x k} where \code{H} is the
#'   forecast horizon and \code{k} is the number of variables.
#' @export
#'
#' @examples
#' \donttest{
#' #homoscedastic with Jeffreys prior
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1, deterministic = "constant")
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    Jeffreys = TRUE)
#'                    
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1)
#'                 
#' forecast(bvar_obj,
#'          fcst_type = "mean",
#'          pi = 0.90,
#'          show_all = TRUE,
#'          ss = TRUE,
#'          ss_type = "mean",
#'          ss_ci = 0.95)
#' }
forecast <- function(x, pi = 0.95, fcst_type = c("mean", "median"),
                     growth_rate_idx = NULL, plot_idx = NULL, show_all = TRUE,
                     ss = TRUE, ss_type = c("mean", "median"), ss_ci = 0.95) {
  
  fcst_type <- match.arg(fcst_type)
  ss_type <- match.arg(ss_type)
  
  yt   <- x$data
  freq <- frequency(yt)
  k    <- ncol(yt)
  
  if (is.null(plot_idx)) {plot_idx <- 1:ncol(yt)}

  posterior    <- rstan::extract(x$fit$stan)
  y_pred       <- posterior$y_pred
  S            <- dim(y_pred)[1]
  alpha        <- 1 - pi
  alpha_ss     <- 1 - ss_ci
  N            <- nrow(yt)
  H            <- x$predict$H
  
  if (isTRUE(ss)) {
    Psi_draws <- posterior$Psi
    dt_full   <- x$setup$dt
    d_pred    <- x$predict$d_pred
    dt_all    <- rbind(dt_full, d_pred)
    
    mu_draws <- array(NA, dim = c(S, N+H, k))
    for (s in 1:S) {
      mu_draws[s, , ] <- dt_all %*% t(Psi_draws[s, , ])
    }
  }
  
  forecast_ret <- matrix(NA, H, k)
  lower_ret    <- matrix(NA, H, k)
  upper_ret    <- matrix(NA, H, k)
  colnames(forecast_ret) <- colnames(yt)
  colnames(lower_ret)    <- colnames(yt)
  colnames(upper_ret)    <- colnames(yt)
  
  point_fcst <- apply(y_pred, c(2, 3), fcst_type)
  fcst_lower <- apply(y_pred, c(2, 3), quantile, probs = alpha / 2)
  fcst_upper <- apply(y_pred, c(2, 3), quantile, probs = 1 - alpha / 2)
  
  time_hist <- as.numeric(time(yt))
  time_fore <- seq(tail(time_hist, 1) + 1 / freq, by = 1 / freq, length.out = H)
  mu_time   <- c(time_hist, time_fore)
  
  for (i in 1:k) {
    yt_i <- yt[, i]
    
    if (isTRUE(ss)) {
      
      if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
        
        mu_scaled     <- mu_draws[, , i] * freq
        mu_line       <- apply(mu_scaled, 2, ss_type)
        mu_line_lower <- apply(mu_scaled, 2, quantile, probs = alpha_ss / 2)
        mu_line_upper <- apply(mu_scaled, 2, quantile, probs = 1 - alpha_ss / 2)
        
      } else {
        
        mu_line       <- apply(mu_draws[, , i], 2, ss_type)
        mu_line_lower <- apply(mu_draws[, , i], 2, quantile, probs = alpha_ss / 2)
        mu_line_upper <- apply(mu_draws[, , i], 2, quantile, probs = 1 - alpha_ss / 2)
      }
    }
    if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      
      annual_hist <- rep(NA, length(yt_i))
      for (t in freq:length(yt_i)) {
        annual_hist[t] <- sum(yt_i[(t - (freq - 1)):t])
      }
      annual_hist <- ts(annual_hist, start = start(yt), frequency = freq)
      
      last_obs <- tail(yt_i, (freq - 1))
      
      pred_extended <- matrix(NA, nrow = S, ncol = (freq - 1) + H)
      for (s in 1:S) {
        pred_extended[s, ] <- c(last_obs, y_pred[s, ,i])
      }
      
      annual_pred <- matrix(NA, nrow = S, ncol = H)
      for (s in 1:S) {
        for (h in 1:H) {
          annual_pred[s, h] <- sum(pred_extended[s, h:(h + (freq - 1))])
        }
      }
      
      annual_fcst  <- apply(annual_pred, 2, fcst_type)
      annual_lower <- apply(annual_pred, 2, quantile, probs = alpha / 2)
      annual_upper <- apply(annual_pred, 2, quantile, probs = 1 - alpha / 2)
      
      forecast_ret[, i] <- annual_fcst
      lower_ret[, i]    <- annual_lower
      upper_ret[, i]    <- annual_upper
      
      if (i %in% plot_idx) {
        time_full  <- c(tail(time_hist, 1), time_fore)
        m_full     <- c(tail(annual_hist, 1), annual_fcst)
        lower_full <- c(tail(annual_hist, 1), annual_lower)
        upper_full <- c(tail(annual_hist, 1), annual_upper)
        
        if (isFALSE(show_all)) {
          xlim_vals    <- c(head(time_fore, 1) - 8, tail(time_fore, 1))
          hist_in_plot <- annual_hist[time_hist >= xlim_vals[1] & time_hist <= xlim_vals[2]]
          ylim_vals    <- c(hist_in_plot, annual_lower, annual_upper)
          if (isTRUE(ss)) {
            mu_in_plot_lower <- mu_line_lower[mu_time >= xlim_vals[1] & mu_time <= xlim_vals[2]]
            mu_in_plot_upper <- mu_line_upper[mu_time >= xlim_vals[1] & mu_time <= xlim_vals[2]]
            ylim_vals <- c(ylim_vals, mu_in_plot_lower, mu_in_plot_upper)
          }
          ylim <- range(ylim_vals, na.rm = TRUE)
          
          plot.ts(annual_hist, main = paste(colnames(yt)[i], "(annual)"),
                  xlab = "Time", ylab = NULL,
                  xlim = xlim_vals,
                  ylim = ylim, col = "black", lwd = 2)
          grid(col = "gray", lty = "dotted") 
          points(as.numeric(time_hist), annual_hist, pch = 16, col = "black", cex=0.8)
          points(time_full[-1], m_full[-1], pch = 16, col = "blue", cex=0.8)
        } else {
          ylim_vals <- c(upper_full, lower_full, annual_hist)
          if (isTRUE(ss)) ylim_vals <- c(ylim_vals, mu_line_lower, mu_line_upper)
          plot.ts(annual_hist, main = paste(colnames(yt)[i], "(annual)"),
                  xlab = "Time", ylab = NULL, col = "black", lwd = 2,
                  xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                  ylim = range(ylim_vals, na.rm = TRUE))
          grid(col = "gray", lty = "dotted") 
        }
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full, col = "blue", lwd = 2)
        if (isTRUE(ss)) {
          polygon(c(mu_time, rev(mu_time)), c(mu_line_upper, rev(mu_line_lower)),
                  col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
          lines(mu_time, mu_line, col = "gray40", lwd = 2, lty = "dashed")
          legend("bottomleft",
                 legend = c(paste0("Forecast (", 100 * pi, "% PI)"),
                            paste0("Posterior steady-state (", 100 * ss_ci, "% CI)")),
                 col    = c("blue", "gray40"),
                 lty    = c(1, 2),
                 lwd    = 2,
                 bty    = "o",
                 bg     = "white")
        } else {
          legend("bottomleft",
                 legend = paste0("Forecast (", 100 * pi, "% PI)"),
                 col    = "blue",
                 lty    = 1,
                 lwd    = 2,
                 bty    = "o",
                 bg     = "white")
        }
      }
      
    } else {
      
      forecast_ret[, i] <- point_fcst[,i]
      lower_ret[, i]    <- fcst_lower[,i]
      upper_ret[, i]    <- fcst_upper[,i]
      
      if (i %in% plot_idx) {
        time_full  <- c(tail(time_hist, 1), time_fore)
        m_full     <- c(tail(yt_i, 1), point_fcst[, i])
        lower_full <- c(tail(yt_i, 1), fcst_lower[, i])
        upper_full <- c(tail(yt_i, 1), fcst_upper[, i])
        
        if (isFALSE(show_all)) {
          xlim_vals    <- c(head(time_fore, 1) - 8, tail(time_fore, 1))
          hist_in_plot <- yt_i[time_hist >= xlim_vals[1] & time_hist <= xlim_vals[2]]
          ylim_vals    <- c(hist_in_plot, lower_full, upper_full)
          if (isTRUE(ss)) {
            mu_in_plot_lower <- mu_line_lower[mu_time >= xlim_vals[1] & mu_time <= xlim_vals[2]]
            mu_in_plot_upper <- mu_line_upper[mu_time >= xlim_vals[1] & mu_time <= xlim_vals[2]]
            ylim_vals <- c(ylim_vals, mu_in_plot_lower, mu_in_plot_upper)
          }
          ylim <- range(ylim_vals, na.rm = TRUE)
          
          plot.ts(yt_i, main = colnames(yt)[i], xlab = "Time", ylab = NULL,
                  xlim = xlim_vals,
                  ylim = ylim, col = "black", lwd = 2)
          grid(col = "gray", lty = "dotted") 
          points(as.numeric(time_hist), yt_i, pch = 16, col = "black", cex=0.8)
          points(time_full[-1], m_full[-1], pch = 16, col = "blue", cex=0.8)
        } else {
          ylim_vals <- c(upper_full, lower_full, yt_i)
          if (isTRUE(ss)) ylim_vals <- c(ylim_vals, mu_line_lower, mu_line_upper)
          plot.ts(yt_i, main = colnames(yt)[i], xlab = "Time", ylab = NULL,
                  col = "black", lwd = 2,
                  xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                  ylim = range(ylim_vals, na.rm = TRUE))
          grid(col = "gray", lty = "dotted") 
        }
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full, col = "blue", lwd = 2)
        if (isTRUE(ss)) {
          polygon(c(mu_time, rev(mu_time)), c(mu_line_upper, rev(mu_line_lower)),
                  col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
          lines(mu_time, mu_line, col = "gray40", lwd = 2, lty = "dashed")
          legend("bottomleft",
                 legend = c(paste0("Forecast (", 100 * pi, "% PI)"),
                            paste0("Posterior steady-state (", 100 * ss_ci, "% CI)")),
                 col    = c("blue", "gray40"),
                 lty    = c(1, 2),
                 lwd    = 2,
                 bty    = "o",
                 bg     = "white")
        } else {
          legend("bottomleft",
                 legend = paste0("Forecast (", 100 * pi, "% PI)"),
                 col    = "blue",
                 lty    = 1,
                 lwd    = 2,
                 bty    = "o",
                 bg     = "white")
        }
      }
    }
  }
  
  invisible(list(forecast = forecast_ret, lower = lower_ret, upper = upper_ret))
}