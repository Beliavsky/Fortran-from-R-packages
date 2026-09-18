#' Impulse Response Functions for a fitted steady-state BVAR model
#'
#' Computes and plots impulse response functions (IRFs) from a fitted
#' steady-state \code{bvar} object. Supports both orthogonalized (OIRF) and generalized
#' (GIRF) impulse responses, with optional conversion to annual growth rates.
#'
#' @param x A steady-state \code{bvar} object that has been passed through \code{\link{fit}}.
#' @param H Integer. The forecast horizon for the IRF. Default \code{20}.
#' @param response Integer vector. Indices of the response variables to plot. If
#'   \code{NULL} (default), all responses are plotted.
#' @param impulse Integer vector. Indices of the impulse variables to plot. If \code{NULL}
#'   (default), all impulses are plotted.
#' @param type Character. Whether to use \code{"median"} or \code{"mean"} as
#'   the point estimate. Default \code{"median"}.
#' @param method Character. The IRF method: \code{"OIRF"} for orthogonalized
#'   or \code{"GIRF"} for generalized impulse responses. Default \code{"OIRF"}.
#' @param ci Numeric. The credible interval width. Default \code{0.95}, i.e. 95%.
#' @param t Integer. Time index for the covariance matrix when using stochastic
#'   volatility models. If \code{NULL} (default), the last time \code{t} is used.
#' @param growth_rate_idx Integer vector. Indices of variables for which the impulse
#'   response is converted from a quarterly or monthly log first difference to an
#'   annual growth rate response, i.e. \eqn{100 (\ln x_{t} - \ln x_{t-f})}, where \eqn{f} is
#'   the frequency of the data (4 for quarterly, 12 for monthly).
#'   Only suitable for variables specified as \eqn{100 (\ln x_{t} - \ln x_{t-1})}, i.e. 
#'   \code{100*diff(log(x))}. Computed by summing up to \eqn{f} periods of the impulse response. Default is \code{NULL}.
#'
#' @return Invisibly returns a list with three arrays: the point estimate IRF, \code{lower}, and
#'   \code{upper} credible bounds, each of dimension \code{k x k x (H+1)}.
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
#' IRF(bvar_obj)
#' }
IRF <- function(x, H = 20, response = NULL, impulse = NULL,
                type = c("median", "mean"), method = c("OIRF", "GIRF"),
                ci = 0.95, t = NULL, growth_rate_idx = NULL) {
  
  method <- match.arg(method)
  type   <- match.arg(type)
  freq   <- frequency(x$data)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  compute_OIRF <- function(A, Sigma, e, N) {
    k            <- nrow(Sigma)
    OIRF_matrix  <- matrix(NA, k * k, N + 1)
    count        <- 1
    P            <- t(chol(Sigma))
    for (j in 1:k) {
      for (i in 1:k) {
        for (n in 0:N) {
          OIRF_matrix[count, n + 1] <- (A[[n + 1]] %*% P %*% e[[j]])[i]
        }
        count <- count + 1
      }
    }
    return(OIRF_matrix)
  }
  
  compute_GIRF <- function(A, Sigma, e, N) {
    k           <- nrow(Sigma)
    GIRF_matrix <- matrix(NA, k * k, N + 1)
    count       <- 1
    for (j in 1:k) {
      for (i in 1:k) {
        for (n in 0:N) {
          GIRF_matrix[count, n + 1] <- (Sigma[j, j]^(-1/2) * A[[n + 1]] %*% Sigma %*% e[[j]])[i]
        }
        count <- count + 1
      }
    }
    return(GIRF_matrix)
  }
  
  k         <- x$setup$k
  p         <- x$setup$p
  var_names <- if (is.null(colnames(x$data))) paste0("y_", 1:k) else colnames(x$data)
  
  stan_draws <- rstan::extract(x$fit$stan)
  N_draws    <- dim(stan_draws$beta)[1]
  
  irf_array <- array(NA, dim = c(N_draws, k, k, H + 1))
  
  e <- vector("list", k)
  for (nE in 1:k) {
    tmp      <- rep(0, k)
    tmp[nE]  <- 1
    e[[nE]]  <- tmp
  }
  
  for (draw in 1:N_draws) {
    Phi <- t(stan_draws$beta[draw, , ])
    if (length(dim(stan_draws$Sigma_u)) == 4) {
      t_use <- if (is.null(t)) dim(stan_draws$Sigma_u)[2] else t - x$setup$p
      Sigma  <- stan_draws$Sigma_u[draw, t_use, , ]
    } else {
      Sigma <- stan_draws$Sigma_u[draw, , ]
    }
    
    Psi  <- MTS::VARpsi(Phi, lag = H)$psi
    kk   <- ncol(Psi) / nrow(Psi)
    A    <- vector("list", length = kk)
    for (nA in 1:kk) {
      col_idx  <- ((nA - 1) * nrow(Psi) + 1):(nA * nrow(Psi))
      A[[nA]]  <- Psi[, col_idx]
    }
    
    tmp_IRF <- if (method == "OIRF") compute_OIRF(A, Sigma, e, H) else compute_GIRF(A, Sigma, e, H)
    
    for (j in 1:k) {
      for (i in 1:k) {
        irf_array[draw, i, j, ] <- tmp_IRF[(j - 1) * k + i, ]
      }
    }
  }
  
  alpha <- 1 - ci
  
  if (!is.null(growth_rate_idx)) {
    for (i in growth_rate_idx) {
      for (j in 1:k) {
        for (draw in 1:N_draws) {
          irf_path <- irf_array[draw, i, j, ]
          all_irf  <- c(rep(0, freq - 1), irf_path)
          tmp      <- numeric(H + 1)
          for (h in 1:(H + 1)) {
            tmp[h] <- sum(all_irf[h:(h + freq - 1)])
          }
          irf_array[draw, i, j, ] <- tmp
        }
      }
    }
  }
  
  m_irf     <- apply(irf_array, c(2, 3, 4), type)
  lower_irf <- apply(irf_array, c(2, 3, 4), quantile, probs = alpha / 2)
  upper_irf <- apply(irf_array, c(2, 3, 4), quantile, probs = 1 - alpha / 2)
  
  horizon    <- 0:H
  type_label <- if (type == "median") "Median" else "Mean"
  
  plot_single <- function(i, j) {
    max_abs    <- max(abs(c(lower_irf[i, j, ], upper_irf[i, j, ], m_irf[i, j, ])))
    ylim_range <- c(-max_abs, max_abs)
    ylab_label <- if (!is.null(growth_rate_idx) && i %in% growth_rate_idx) {
      paste0("Response: ", var_names[i], " (annual)")
    } else {
      paste0("Response: ", var_names[i])
    }
    
    plot(NA, xlim = c(0, H), ylim = ylim_range, type = "n",
         xlab = "Horizon", ylab = ylab_label, font.lab = 2)
    polygon(c(horizon, rev(horizon)),
            c(upper_irf[i, j, ], rev(lower_irf[i, j, ])),
            col = rgb(0, 0, 1, 0.2), border = NA)
    lines(horizon, m_irf[i, j, ], col = "blue", lwd = 2)
    abline(h = 0, col = "black", lty = 1)
    
    main_title <- if (isTRUE(x$priors$SV)) {
      paste0("Posterior ", type_label, " ", method, " (", round(ci * 100), "% CI)\nt=", t, "\nImpulse: ", var_names[j])
    } else {
      paste0("Posterior ", type_label, " ", method, " (", round(ci * 100), "% CI)\n\nImpulse: ", var_names[j])
    }
    title(main = main_title)
  }
  
  if (is.null(response) && is.null(impulse)) {
    par(mfrow = c(k, k))
    for (j in 1:k) for (i in 1:k) plot_single(i, j)
  } else if (is.null(response) && !is.null(impulse)) {
    par(mfrow = c(k, length(impulse)))
    for (j in impulse) for (i in 1:k) plot_single(i, j)
  } else if (!is.null(response) && is.null(impulse)) {
    par(mfrow = c(length(response), k))
    for (i in response) for (j in 1:k) plot_single(i, j)
  } else {
    par(mfrow = c(length(response), length(impulse)))
    for (i in response) for (j in impulse) plot_single(i, j)
  }
  
  if (type == "median") {
    invisible(list(median_irf = m_irf, lower = lower_irf, upper = upper_irf))
  } else {
    invisible(list(mean_irf = m_irf, lower = lower_irf, upper = upper_irf))
  }
}