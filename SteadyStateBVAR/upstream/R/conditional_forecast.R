#' Conditional forecasts from a fitted steady-state BVAR model
#'
#' Computes and plots conditional forecasts from a fitted steady-state \code{bvar} object.
#' Conditions are imposed on specific variables at specific horizons using Algorithm 3.3.1
#' from Dieppe, Legrand, and van Roye (2016). Both conditional and unconditional
#' forecasts are plotted for comparison. Please note that for the moment, conditional forecasting
#' is only enabled for the homoscedastic steady-state BVAR, i.e. when \code{SV=FALSE} in \code{priors()}.
#'
#' @param x A steady-state \code{bvar} object that has been passed through
#'   \code{\link{fit}}.
#' @param conditions A data frame with three columns: \code{var} (integer index
#'   of the conditioned variable), \code{horizon} (integer forecast horizon at
#'   which the condition is imposed), and \code{value} (numeric, the imposed
#'   condition). Note that \code{var} and \code{horizon} must be integers, as
#'   they are used for indexing.
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
#'
#' @return Invisibly returns a list with three matrices: \code{forecast}, \code{lower}, and
#'   \code{upper}, each of dimension \code{H x k}.
#'   
#' @details
#' See Section 5.3 of Dieppe, Legrand, and van Roye (2016) for more details.
#' Please note the limitations of this method, see the detailed discussion in Section 5.4 of Dieppe, Legrand, and van Roye (2016).
#' 
#' @export
#'
#' @references
#' Dieppe, A., Legrand, R., and van Roye, B. (2016).
#' The BEAR toolbox. Working Paper Series, No. 1934. European Central Bank.
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
#'                    Jeffreys = TRUE,
#'                    SV = FALSE,
#'                    SV_type = NULL,
#'                    SV_priors = NULL)
#'                    
#' bvar_obj <- fit(bvar_obj,
#'                 H = 8,
#'                 iter = 200,
#'                 warmup = 50,
#'                 chains = 1,
#'                 cores = 1)
#'                 
#' conditions <- data.frame(var = rep(2,8),
#'                          horizon = rep(1:8),
#'                          value   = rep(1,8))
#'                          
#' (cond_fcst <- conditional_forecast(bvar_obj,
#'                                   conditions,
#'                                   pi=0.68,
#'                                   fcst_type = "mean"))
#' }
conditional_forecast <- function(x, conditions, pi = 0.95,
                                 fcst_type = c("mean", "median"),
                                 growth_rate_idx = NULL, plot_idx = NULL) {
  
  if (isTRUE(x$priors$SV)) {
    stop("conditional_forecast is only supported for the homoscedastic steady-state BVAR (SV = FALSE).")
  }
  
  fcst_fun  <- match.arg(fcst_type)
  posterior <- rstan::extract(x$fit$stan)
  n_draws   <- dim(posterior$beta)[1]
  
  p      <- x$setup$p
  k      <- x$setup$k
  Y      <- x$setup$Y
  X      <- x$setup$X
  N      <- x$setup$N
  d_pred <- x$predict$d_pred
  H      <- x$predict$H
  alpha  <- 1 - pi
  s      <- k * H
  
  cond_forecast_array <- array(NA, dim = c(H, k, n_draws))
  
  for (n in 1:n_draws) {
    
    beta_n  <- posterior$beta[n, , ]
    Sigma_n <- posterior$Sigma_u[n, , ]
    Psi_n   <- posterior$Psi[n, , ]
    D_n     <- t(chol(Sigma_n))
    
    Pi_n <- vector("list", p)
    for (i in 1:p) {
      rows_idx <- ((i - 1) * k + 1):(i * k)
      Pi_n[[i]] <- t(beta_n[rows_idx, ])
    }
    
    Y_pred_mat <- matrix(NA, nrow = H, ncol = k)
    for (h in 1:H) {
      ytilde_t <- d_pred[h, ] %*% t(Psi_n)
      if (h > 1) {
        for (i in 1:min(h - 1, p)) {
          term     <- (Y_pred_mat[h - i, ] - d_pred[h - i, ] %*% t(Psi_n)) %*% t(Pi_n[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      if (h <= p) {
        for (i in h:p) {
          term     <- (Y[N + h - i, ] - X[N + h - i, ] %*% t(Psi_n)) %*% t(Pi_n[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      Y_pred_mat[h, ] <- ytilde_t
    }
    
    Pi_n        <- t(beta_n)
    Phi_n       <- MTS::VARpsi(Pi_n, lag = H)$psi
    n_Phi       <- ncol(Phi_n) / k
    Phi_n_tilde <- lapply(1:n_Phi, function(i) {
      cols <- ((i - 1) * k + 1):(i * k)
      Phi_n[, cols] %*% D_n
    })
    
    v <- nrow(conditions)
    R <- matrix(0, v, s)
    r <- rep(0, v)
    
    for (row_idx in 1:v) {
      i     <- conditions$var[row_idx]
      h     <- conditions$horizon[row_idx]
      value <- conditions$value[row_idx]
      row   <- c()
      for (j in 1:h) {
        row <- c(row, Phi_n_tilde[[h - j + 1]][i, ])
      }
      R[row_idx, 1:length(row)] <- row
      r[row_idx] <- value - Y_pred_mat[h, i]
    }
    
    svd_R  <- svd(R, nu = v, nv = s)
    U      <- svd_R$u
    D      <- diag(svd_R$d)
    V_1    <- svd_R$v[, 1:v]
    V_2    <- svd_R$v[, (v + 1):s]
    lambda <- rnorm(s - v)
    eta    <- V_1 %*% solve(D) %*% t(U) %*% r + V_2 %*% lambda
    
    eta_mat    <- matrix(eta, nrow = k, ncol = H)
    Y_cond_mat <- matrix(NA, nrow = H, ncol = k)
    for (h in 1:H) {
      sum_shocks <- matrix(0, nrow = k, ncol = 1)
      for (j in 1:h) {
        sum_shocks <- sum_shocks + Phi_n_tilde[[h - j + 1]] %*% eta_mat[, j]
      }
      Y_cond_mat[h, ] <- Y_pred_mat[h, ] + sum_shocks
    }
    cond_forecast_array[, , n] <- Y_cond_mat
  }
  
  cond_point <- apply(cond_forecast_array, c(1, 2), fcst_fun)
  cond_lower <- apply(cond_forecast_array, c(1, 2), quantile, probs = alpha / 2)
  cond_upper <- apply(cond_forecast_array, c(1, 2), quantile, probs = 1 - alpha / 2)
  
  y_pred_uncond   <- posterior$y_pred
  y_pred_m_uncond <- apply(y_pred_uncond, c(2, 3), fcst_fun)
  
  data_Y    <- x$data
  freq      <- frequency(data_Y)
  m         <- ncol(data_Y)
  if (is.null(plot_idx)) plot_idx <- 1:m
  
  time_hist <- as.numeric(time(data_Y))
  time_fore <- seq(tail(time_hist, 1) + 1 / freq, by = 1 / freq, length.out = H)
  
  forecast_ret <- matrix(NA, H, m)
  lower_ret    <- matrix(NA, H, m)
  upper_ret    <- matrix(NA, H, m)
  colnames(forecast_ret) <- colnames(data_Y)
  colnames(lower_ret)    <- colnames(data_Y)
  colnames(upper_ret)    <- colnames(data_Y)
  
  for (i in 1:m) {
    smply      <- data_Y[, i]
    fcst_m     <- cond_point[, i]
    fcst_lower <- cond_lower[, i]
    fcst_upper <- cond_upper[, i]
    uncond_m   <- y_pred_m_uncond[, i]
    
    if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      annual_hist <- rep(NA, length(smply))
      for (t in freq:length(smply)) {
        annual_hist[t] <- sum(smply[(t - (freq - 1)):t])
      }
      annual_hist <- ts(annual_hist, start = start(data_Y), frequency = freq)
      
      last_obs <- tail(smply, (freq - 1))
      
      cond_draws_full <- matrix(NA, nrow = n_draws, ncol = (freq - 1) + H)
      for (n in 1:n_draws) {
        cond_draws_full[n, ] <- c(last_obs, cond_forecast_array[, i, n])
      }
      
      cond_annual_draws <- matrix(NA, nrow = n_draws, ncol = H)
      for (n in 1:n_draws) {
        for (h in 1:H) {
          cond_annual_draws[n, h] <- sum(cond_draws_full[n, h:(h + (freq - 1))])
        }
      }
      
      annual_fcst  <- apply(cond_annual_draws, 2, fcst_fun)
      annual_lower <- apply(cond_annual_draws, 2, quantile, probs = alpha / 2)
      annual_upper <- apply(cond_annual_draws, 2, quantile, probs = 1 - alpha / 2)
      
      uncond_draws_full <- matrix(NA, nrow = n_draws, ncol = (freq - 1) + H)
      for (n in 1:n_draws) {
        uncond_draws_full[n, ] <- c(last_obs, y_pred_uncond[n, ,i])
      }
      
      uncond_annual_draws <- matrix(NA, nrow = n_draws, ncol = H)
      for (n in 1:n_draws) {
        for (h in 1:H) {
          uncond_annual_draws[n, h] <- sum(uncond_draws_full[n, h:(h + (freq - 1))])
        }
      }
      
      annual_fcst_uncond <- apply(uncond_annual_draws, 2, fcst_fun)
      
      forecast_ret[, i] <- annual_fcst
      lower_ret[, i]    <- annual_lower
      upper_ret[, i]    <- annual_upper
      
      if (i %in% plot_idx) {
        time_full   <- c(tail(time_hist, 1), time_fore)
        m_full      <- c(tail(annual_hist, 1), annual_fcst)
        lower_full  <- c(tail(annual_hist, 1), annual_lower)
        upper_full  <- c(tail(annual_hist, 1), annual_upper)
        uncond_full <- c(tail(annual_hist, 1), annual_fcst_uncond)
        plot.ts(annual_hist, main = paste(colnames(data_Y)[i], "(annual)"),
                xlab = "Time", ylab = NULL, col = "black", lwd = 2,
                xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                ylim = range(upper_full, lower_full, annual_hist, uncond_full, na.rm = TRUE))
        grid(col = "gray", lty = "dotted") 
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full,      col = "blue", lwd = 2)
        lines(time_full, uncond_full, col = "red",  lwd = 1, lty = 1)
        legend("bottomleft",
               legend = c("Conditional forecast", "Unconditional forecast"),
               col    = c("blue", "red"),
               lwd    = 2,
               bty    = "o",
               bg     = "white")
      }
    } else {
      forecast_ret[, i] <- fcst_m
      lower_ret[, i]    <- fcst_lower
      upper_ret[, i]    <- fcst_upper
      
      if (i %in% plot_idx) {
        time_full   <- c(tail(time_hist, 1), time_fore)
        m_full      <- c(tail(smply, 1), fcst_m)
        lower_full  <- c(tail(smply, 1), fcst_lower)
        upper_full  <- c(tail(smply, 1), fcst_upper)
        uncond_full <- c(tail(smply, 1), uncond_m)
        plot.ts(smply, main = colnames(data_Y)[i], xlab = "Time", ylab = NULL,
                col = "black", lwd = 2,
                xlim = c(head(time_hist, 1), tail(time_fore, 1)),
                ylim = range(upper_full, lower_full, smply, uncond_full))
        grid(col = "gray", lty = "dotted") 
        polygon(c(time_full, rev(time_full)), c(upper_full, rev(lower_full)),
                col = rgb(0, 0, 1, 0.2), border = NA)
        lines(time_full, m_full,      col = "blue", lwd = 2)
        lines(time_full, uncond_full, col = "red",  lwd = 1, lty = 1)
        legend("bottomleft",
               legend = c("Conditional forecast", "Unconditional forecast"),
               col    = c("blue", "red"),
               lwd    = 2,
               bty    = "o",
               bg     = "white")
      }
    }
  }
  
  invisible(list(
    forecast   = forecast_ret,
    lower      = lower_ret,
    upper      = upper_ret
  ))
}