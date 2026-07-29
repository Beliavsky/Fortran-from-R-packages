# ============================================================================
# Internal PSIS-LOO helpers
#
# The previous implementation passed a 3-D array [draws, time, sector] to
# loo::loo(), which interprets a 3-D array as [iterations, chains, observations]
# -> it silently treated the time axis as "chains" and only the S sectors as
# observations. These helpers reshape the per-(t,s) log-likelihood into a proper
# [draws x observations] matrix over the fitted window 2:T_lik, with a chain_id
# so relative_eff() can estimate the relative effective sample size correctly.
# ============================================================================

#' Extract a log-likelihood draws matrix with chain ids (column-named log_lik[t,s])
#' @keywords internal
#' @noRd
.extract_loglik_mat <- function(fit) {
  if (inherits(fit, "CmdStanMCMC")) {
    mat <- fit$draws("log_lik", format = "matrix")
    n_chains <- tryCatch(fit$num_chains(), error = function(e) 1L)
  } else {
    a <- rstan::extract(fit, pars = "log_lik", permuted = FALSE) # [iter, chain, var]
    n_chains <- dim(a)[2]
    vnames <- dimnames(a)[[3]]
    mat <- do.call(rbind, lapply(seq_len(n_chains), function(ci) {
      m <- a[, ci, , drop = FALSE]
      dim(m) <- dim(a)[c(1, 3)]
      colnames(m) <- vnames
      m
    }))
  }
  total <- nrow(mat)
  if (n_chains < 1L || total %% n_chains != 0) n_chains <- 1L
  list(mat = mat, chain_id = rep(seq_len(n_chains), each = total / n_chains))
}

#' Build a PSIS-LOO object from a log-lik matrix over the window 2:T_lik
#' @keywords internal
#' @noRd
.loo_from_loglik <- function(ll_mat, chain_id, T_lik, S) {
  t_idx <- 2:T_lik
  want <- as.vector(vapply(
    seq_len(S),
    function(s) sprintf("log_lik[%d,%d]", t_idx, s),
    character(length(t_idx))
  ))
  miss <- setdiff(want, colnames(ll_mat))
  if (length(miss) > 0L) {
    stop(sprintf("Missing %d expected log_lik columns (e.g. '%s').",
                 length(miss), miss[1]), call. = FALSE)
  }
  ll_used <- ll_mat[, want, drop = FALSE]
  r_eff <- tryCatch(
    loo::relative_eff(exp(ll_used), chain_id = chain_id),
    error = function(e) NULL
  )
  loo::loo(ll_used, r_eff = r_eff)
}

#' @keywords internal
#' @noRd
.compute_loo <- function(fit, T_lik, S) {
  ex <- .extract_loglik_mat(fit)
  .loo_from_loglik(ex$mat, ex$chain_id, T_lik, S)
}

#' Summarize Pareto-k diagnostics from a loo object
#' @keywords internal
#' @noRd
.summarize_pareto_k <- function(loo_res) {
  k <- tryCatch(loo_res$diagnostics$pareto_k, error = function(e) NULL)
  if (is.null(k)) return(NULL)
  list(
    n = length(k),
    max = max(k, na.rm = TRUE),
    prop_gt_0.7 = mean(k > 0.7, na.rm = TRUE),
    prop_gt_1   = mean(k > 1,   na.rm = TRUE)
  )
}


#' Count HMC divergences
#'
#' Extracts the number of divergent transitions from a fitted Stan model.
#'
#' @param fit Fitted Stan model object (CmdStanMCMC or stanfit)
#'
#' @return Integer. Number of divergent transitions (post-warmup).
#'
#' @examples
#' \donttest{
#' # Create a "mock" CmdStanMCMC object for demonstration
#' # (This simulates a model with 0 divergences)
#' mock_fit <- structure(list(
#'   sampler_diagnostics = function() {
#'     # Return a 3D array: [iterations, chains, variables]
#'     # Variable 1 is usually accept_stat__, let's say var 2 is divergent__
#'     ar <- array(0, dim = c(100, 4, 6)) 
#'     dimnames(ar)[[3]] <- c("accept_stat__", "divergent__", "energy__", 
#'                            "n_leapfrog__", "stepsize__", "treedepth__")
#'     return(ar)
#'   }
#' ), class = "CmdStanMCMC")
#' 
#' # Now the example can run without errors:
#' n_div <- count_divergences(mock_fit)
#' print(n_div)
#' }
#'
#' @export
count_divergences <- function(fit) {
  tryCatch({
    if (inherits(fit, "CmdStanMCMC")) {
      sd <- fit$sampler_diagnostics()
      vnames <- dimnames(sd)[[3]]
      if (!("divergent__" %in% vnames)) {
        return(NA_integer_)
      }
      i <- match("divergent__", vnames)
      return(sum(sd[, , i], na.rm = TRUE))
    } else {
      sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      return(sum(vapply(sp, function(m) sum(m[, "divergent__"]), numeric(1))))
    }
  }, error = function(e) {
    if (inherits(fit, "CmdStanMCMC")) {
      csvs <- tryCatch(fit$output_files(), error = function(e) NULL)
      if (is.null(csvs) || length(csvs) == 0) {
        return(NA_integer_)
      }
      sd2 <- cmdstanr::read_cmdstan_csv(csvs)$sampler_diagnostics
      if (length(dim(sd2)) == 3 && "divergent__" %in% dimnames(sd2)[[3]]) {
        return(sum(sd2[, , "divergent__"], na.rm = TRUE))
      }
    }
    NA_integer_
  })
}


#' Evaluate out-of-sample forecast metrics
#'
#' Computes RMSE and MAE for multiple forecast horizons.
#'
#' @param summ List. Posterior summary from \code{\link{extract_posterior_summary}}
#' @param Yz Numeric matrix. Standardized Y values (T x S)
#' @param Xz Numeric matrix. Standardized X values (T x S)
#' @param zTMG Numeric vector. Standardized TMG series
#' @param T_train Integer. End of training period
#' @param COM_ts Numeric matrix. COM values by time and sector (T x S)
#' @param K_ts Numeric matrix. Capital values by time and sector (T x S)
#' @param com_in_mean Logical. Whether COM is included in mean equation
#' @param horizons Integer vector. Forecast horizons to evaluate
#'
#' @return Named list with one element per horizon, each a list with
#'   \code{h}, \code{RMSE}, \code{MAE} and \code{n_obs} (number of pooled
#'   sector-by-origin errors; \code{0} / \code{NA} when the horizon exceeds the
#'   test window).
#'
#' @details
#' Two caveats matter for interpretation:
#' \enumerate{
#'   \item \strong{Conditional forecast.} The recursion uses the \emph{realized}
#'     future covariates \eqn{X_{t-1}} and \eqn{zTMG_t} at each step. It is
#'     therefore a forecast of \eqn{Y} \emph{conditional on} the future paths of
#'     \eqn{X} and TMG, not an unconditional forecast. If those covariates are
#'     not known ex ante, treat these numbers as conditional/nowcasting metrics.
#'   \item \strong{Plug-in medians.} The path is propagated with posterior
#'     medians of the parameters and ignores parameter uncertainty, the
#'     stochastic volatility and the Student-t innovations. It is a point
#'     (mean-equation) forecast, not the full posterior predictive distribution.
#'     For genuinely out-of-sample numbers, fit with \code{fit_window = "train"}.
#' }
#'
#' @examples
#' # 1. Generate dummy data for testing
#' T_obs <- 20
#' S <- 2
#' Yz <- matrix(rnorm(T_obs * S), nrow = T_obs, ncol = S)
#' Xz <- matrix(rnorm(T_obs * S), nrow = T_obs, ncol = S)
#' COM_ts <- matrix(abs(rnorm(T_obs * S)), nrow = T_obs, ncol = S)
#' K_ts <- matrix(abs(rnorm(T_obs * S)) + 1, nrow = T_obs, ncol = S)
#' zTMG <- rnorm(T_obs)
#' 
#' # 2. Create a dummy summary list (mimicking extract_posterior_summary)
#' summ <- list(
#'   theta_s = runif(S),
#'   kappa_s = runif(S),
#'   a3_s = runif(S),
#'   beta0_s = runif(S),
#'   beta1 = 0.5,
#'   gamma = 0.1
#' )
#' 
#' # 3. Run the function
#' metrics <- evaluate_oos(summ, Yz, Xz, zTMG, T_train = 15, 
#'                         COM_ts, K_ts, horizons = c(1, 2))
#' print(metrics)
#'
#' @export
evaluate_oos <- function(summ, Yz, Xz, zTMG, T_train,
                         COM_ts, K_ts, com_in_mean = FALSE,
                         horizons = c(1, 4, 8)) {
  
  Tn <- nrow(Yz)
  S <- ncol(Yz)
  
  com_wmean <- numeric(S)
  com_wsd <- numeric(S)
  
  for (s in seq_len(S)) {
    denom <- sum(K_ts[seq_len(T_train), s], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- 1
    w <- K_ts[seq_len(T_train), s] / denom
    com_wmean[s] <- sum(COM_ts[seq_len(T_train), s] * w, na.rm = TRUE)
    v <- sum(w * (COM_ts[seq_len(T_train), s] - com_wmean[s])^2, na.rm = TRUE)
    com_wsd[s] <- sqrt(max(v, 1e-16))
  }
  
  res <- lapply(horizons, function(hh) {
    last_origin <- Tn - hh + 1
    if (last_origin < (T_train + 1)) {
      # Horizon longer than the available test window: nothing to evaluate.
      # (Previously this produced a descending ':' sequence and out-of-range
      # indexing.)
      return(list(h = hh, RMSE = NA_real_, MAE = NA_real_, n_obs = 0L))
    }
    errs <- c()

    for (t in seq.int(T_train + 1, last_origin)) {
      if (t - 1 < 1) next
      
      y_pred <- Yz[t - 1, ]
      
      for (h in seq_len(hh)) {
        t_pred <- t - 1 + h
        if (t_pred > Tn) break
        
        ztmg <- if (t_pred <= Tn) zTMG[t_pred] else zTMG[Tn]
        
        for (s in seq_len(S)) {
          zlag <- y_pred[s] - summ$theta_s[s]
          drift <- summ$kappa_s[s] * (
            summ$theta_s[s] - y_pred[s] + summ$a3_s[s] * zlag^3
          )
          xterm <- (summ$beta0_s[s] + summ$beta1 * ztmg) *
            (if (t_pred <= Tn) Xz[t_pred - 1, s] else Xz[Tn, s])
          
          com_std_denom <- if (com_wsd[s] > 0) com_wsd[s] else 1
          com_std <- (COM_ts[min(t_pred - 1, Tn), s] - com_wmean[s]) / com_std_denom
          
          mu <- drift + xterm + (if (com_in_mean) summ$gamma * com_std else 0)
          y_pred[s] <- y_pred[s] + mu
        }
      }
      
      errs <- c(errs, Yz[t + hh - 1, ] - y_pred)
    }
    
    list(
      h = hh,
      RMSE = if (length(errs)) sqrt(mean(errs^2, na.rm = TRUE)) else NA_real_,
      MAE  = if (length(errs)) mean(abs(errs), na.rm = TRUE) else NA_real_,
      n_obs = length(errs)
    )
  })
  
  names(res) <- paste0("h", horizons)
  res
}


#' Validate OU model fit
#'
#' Prints diagnostic summaries and hypothesis tests for a fitted model.
#'
#' @param fit_res List returned by \code{\link{fit_ou_nonlinear_tmg}}
#' @param verbose Logical. Print detailed output. Default TRUE.
#'
#' @return Invisibly returns the diagnostics list
#'
#' @examples
#' \donttest{
#' # Create a dummy results list that mimics the output of fit_ou_nonlinear_tmg
#' dummy_results <- list(
#'   diagnostics = list(
#'     rhat = c(alpha = 1.01, beta = 1.00),
#'     ess = c(alpha = 400, beta = 350),
#'     loo = list(estimates = matrix(c(1, 0.1), ncol=2, 
#'                dimnames=list("elpd_loo", c("Estimate", "SE")))),
#'     oos = list(h1 = list(RMSE = 0.5))
#'   ),
#'   factor_ou = list(beta1 = 0.3),
#'   nonlinear = list(a3 = -0.5),
#'   sv = list(rho_s = 0.2)
#' )
#' 
#' # Run validation on the dummy object
#' validate_ou_fit(dummy_results)
#' }
#'
#' @export
validate_ou_fit <- function(fit_res, verbose = TRUE) {
  # ---- Hard structural checks (this is a real validator, not just a report) ----
  if (!is.list(fit_res)) {
    stop("`fit_res` must be a list returned by fit_ou_nonlinear_tmg().",
         call. = FALSE)
  }
  required <- c("factor_ou", "diagnostics")
  missing_top <- setdiff(required, names(fit_res))
  if (length(missing_top) > 0) {
    stop("`fit_res` is missing required components: ",
         paste(missing_top, collapse = ", "), call. = FALSE)
  }
  dg <- fit_res$diagnostics
  if (is.null(dg$rhat)) {
    warning("No R-hat found in diagnostics; the fit may be incomplete.",
            call. = FALSE)
  }

  if (verbose) {
    message("\n==== NONLINEAR OU MODEL VALIDATION ====")

    message("\n-- MCMC convergence --")
    message(sprintf("Max R-hat: %.4f   share(R-hat > 1.01): %.3f",
                    dg$rhat_max %||% NA_real_, dg$rhat_share %||% NA_real_))
    if (!is.null(dg$rhat_max) && is.finite(dg$rhat_max) && dg$rhat_max > 1.01) {
      message("  WARNING: R-hat > 1.01 -> chains may not have converged.")
    }
    message(sprintf("Divergences: %s", dg$divergences %||% NA))
    message("R-hat (head):"); print(utils::head(dg$rhat))
    message("ESS (head):");   print(utils::head(dg$ess))

    if (!is.null(dg$loo)) {
      message("\n-- PSIS-LOO --")
      print(dg$loo)
      pk <- dg$loo_pareto_k
      if (!is.null(pk)) {
        message(sprintf(
          "Pareto-k: max = %.2f | share > 0.7 = %.3f | share > 1 = %.3f",
          pk$max, pk$prop_gt_0.7, pk$prop_gt_1))
        if (pk$prop_gt_0.7 > 0.05) {
          message("  WARNING: many Pareto-k > 0.7 -> PSIS-LOO unreliable ",
                  "(expected with one latent SV state per observation; ",
                  "consider LFO-CV).")
        }
      }
    }

    message("\n-- OOS metrics --")
    print(dg$oos)

    beta1 <- fit_res$factor_ou$beta1
    message(sprintf(
      "\nH1: beta1 (TMG effect) median point estimate: %.4f (%s)",
      beta1, if (beta1 > 0) ">0" else "<=0"
    ))
    
    a3_med <- fit_res$nonlinear$a3
    message(sprintf(
      "H4: a3<0 (increasing restoring force): %.2f proportion of sectors",
      mean(a3_med < 0, na.rm = TRUE)
    ))
    
    rho <- fit_res$sv$rho_s
    message(sprintf("H6: median rho_s: %.3f", stats::median(rho, na.rm = TRUE)))
  }
  
  invisible(dg)
}


#' Compare models using PSIS-LOO
#'
#' Compares two fitted models using PSIS-LOO cross-validation.
#'
#' @param results_new List. Results from new model.
#' @param results_base List. Results from base model.
#'
#' @return List with components:
#'   \describe{
#'     \item{loo_table}{Comparison table from loo::loo_compare}
#'     \item{deltaELPD}{Numeric difference in ELPD}
#'   }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   # 1. Create mock 'loo' objects manually to avoid computation errors
#'   # Structure required by loo_compare: list with 'estimates' and class 'psis_loo'
#'   
#'   # Mock Model 1: ELPD = -100
#'   est1 <- matrix(c(-100, 5), nrow = 1, dimnames = list("elpd_loo", c("Estimate", "SE")))
#'   # Pointwise data (required for diff calculation). 10 observations.
#'   pw1 <- matrix(c(rep(-10, 10)), ncol = 1, dimnames = list(NULL, "elpd_loo"))
#'   
#'   loo_obj1 <- list(estimates = est1, pointwise = pw1)
#'   class(loo_obj1) <- c("psis_loo", "loo")
#'   
#'   # Mock Model 2: ELPD = -102 (worse)
#'   est2 <- matrix(c(-102, 5), nrow = 1, dimnames = list("elpd_loo", c("Estimate", "SE")))
#'   pw2 <- matrix(c(rep(-10.2, 10)), ncol = 1, dimnames = list(NULL, "elpd_loo"))
#'   
#'   loo_obj2 <- list(estimates = est2, pointwise = pw2)
#'   class(loo_obj2) <- c("psis_loo", "loo")
#'   
#'   # 2. Wrap in the structure expected by your package
#'   res_new <- list(diagnostics = list(loo = loo_obj1))
#'   res_base <- list(diagnostics = list(loo = loo_obj2))
#'   
#'   # 3. Compare (This will now run cleanly without warnings)
#'   cmp <- compare_models_loo(res_new, res_base)
#'   print(cmp)
#' }
#' }
#'
#' @export
compare_models_loo <- function(results_new, results_base) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required for model comparison.")
  }
  
  loo_new <- results_new$diagnostics$loo %||% results_new$loo
  loo_base <- results_base$diagnostics$loo %||% results_base$loo

  if (is.null(loo_new) || is.null(loo_base)) {
    stop("Both results must contain a 'loo' object (diagnostics$loo).",
         call. = FALSE)
  }
  if (!inherits(loo_new, "loo") || !inherits(loo_base, "loo")) {
    stop("The 'loo' components must be objects of class 'loo'.", call. = FALSE)
  }
  n_new <- nrow(loo_new$pointwise)
  n_base <- nrow(loo_base$pointwise)
  if (!is.null(n_new) && !is.null(n_base) && n_new != n_base) {
    stop(sprintf(paste0("LOO objects have different numbers of observations ",
                        "(%d vs %d); they are not comparable. Ensure both models ",
                        "use the same data and fit_window."), n_new, n_base),
         call. = FALSE)
  }

  cmp <- loo::loo_compare(loo_new, loo_base)

  # Unambiguous new-minus-base difference (positive => new model preferred),
  # computed from the estimates rather than relying on loo_compare row order.
  elpd_new <- loo_new$estimates["elpd_loo", "Estimate"]
  elpd_base <- loo_base$estimates["elpd_loo", "Estimate"]

  list(
    loo_table = cmp,
    deltaELPD = as.numeric(elpd_new - elpd_base)
  )
}
