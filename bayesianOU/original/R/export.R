#' Export model comparison to Excel
#'
#' Creates an Excel workbook with model comparison results,
#' parameter summaries, and fit information.
#'
#' @param results_new List. Results from new model.
#' @param results_base List. Results from base model.
#' @param path Character. Output file path. Default "model_comparison.xlsx".
#' @param verbose Logical. Print progress messages. Default FALSE.
#'
#' @return TRUE invisibly on success.
#'
#' @details
#' Creates three worksheets:
#' \itemize{
#'   \item loo_comparison: PSIS-LOO comparison table
#'   \item param_summary: Sector-specific parameter estimates
#'   \item fit_info: Model configuration and diagnostics
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("openxlsx", quietly = TRUE) &&
#'     requireNamespace("loo", quietly = TRUE)) {
#'   
#'   # 1. Create mock results objects
#'   # Mock Model New
#'   res_new <- list(
#'     factor_ou = list(
#'       kappa_s = c(0.5, 0.6), a3_s = c(-0.1, -0.2), beta0_s = c(1, 2),
#'       gamma = 0.05, model = "TestModel", beta1 = 0.3, nu = 4,
#'       factor_ou_info = list(T_train = 50, com_in_mean = TRUE)
#'     ),
#'     diagnostics = list(
#'       divergences = 0,
#'       # Mock LOO object
#'       loo = list(
#'          estimates = matrix(c(-100, 2), 1, 2, dimnames=list("elpd_loo", c("Estimate","SE"))),
#'          pointwise = matrix(rep(-2, 50), ncol=1)
#'       )
#'     )
#'   )
#'   class(res_new$diagnostics$loo) <- c("psis_loo", "loo")
#'   
#'   # Mock Model Base
#'   res_base <- list(
#'     diagnostics = list(
#'       loo = list(
#'          estimates = matrix(c(-110, 2), 1, 2, dimnames=list("elpd_loo", c("Estimate","SE"))),
#'          pointwise = matrix(rep(-2.2, 50), ncol=1)
#'       )
#'     )
#'   )
#'   class(res_base$diagnostics$loo) <- c("psis_loo", "loo")
#'
#'   # 2. Define a safe temporary path
#'   out_path <- file.path(tempdir(), "comparison.xlsx")
#'   
#'   # 3. Run export (This writes to tempdir, allowed by CRAN)
#'   try({
#'     export_model_comparison(res_new, res_base, path = out_path)
#'     # unlink(out_path) # Cleanup
#'   })
#' }
#' }
#'
#' @export
export_model_comparison <- function(results_new, results_base,
                                    path = file.path(tempdir(), "model_comparison.xlsx"),
                                    verbose = FALSE) {
  
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for Excel export.")
  }
  
  vmsg("Creating Excel workbook", verbose)
  wb <- openxlsx::createWorkbook()
  
  vmsg("Adding LOO comparison sheet", verbose)
  openxlsx::addWorksheet(wb, "loo_comparison")
  cmp <- tryCatch(
    compare_models_loo(results_new, results_base)$loo_table,
    error = function(e) {
      warning("Could not compare LOO: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(cmp)) {
    openxlsx::writeData(wb, "loo_comparison", as.data.frame(cmp))
  }
  
  vmsg("Adding parameter summary sheet", verbose)
  openxlsx::addWorksheet(wb, "param_summary")
  kappa <- results_new$factor_ou$kappa_s
  a3 <- results_new$factor_ou$a3_s
  beta0 <- results_new$factor_ou$beta0_s
  
  ps <- data.frame(
    sector = seq_along(kappa),
    kappa = kappa,
    a3 = a3,
    beta0 = beta0
  )
  
  gamma_global <- results_new$factor_ou$gamma
  if (!is.null(gamma_global) && is.finite(gamma_global)) {
    ps$gamma_global <- rep(gamma_global, length(kappa))
  }
  openxlsx::writeData(wb, "param_summary", ps)
  
  vmsg("Adding fit info sheet", verbose)
  openxlsx::addWorksheet(wb, "fit_info")
  info <- data.frame(
    model = as.character(results_new$factor_ou$model %||% NA),
    T_train = as.integer(
      results_new$factor_ou$factor_ou_info$T_train %||% NA
    ),
    com_in_mean = as.logical(
      results_new$factor_ou$factor_ou_info$com_in_mean %||% NA
    ),
    beta1 = as.numeric(results_new$factor_ou$beta1 %||% NA),
    nu = as.numeric(results_new$factor_ou$nu %||% NA),
    divergences = as.integer(results_new$diagnostics$divergences %||% NA)
  )
  openxlsx::writeData(wb, "fit_info", info)
  
  vmsg(sprintf("Saving workbook to %s", path), verbose)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  
  message(sprintf("Saved: %s", path))
  invisible(TRUE)
}


#' Export posterior draws to CSV
#'
#' Saves posterior draws for selected parameters to a CSV file.
#'
#' @param fit Fitted Stan model object
#' @param params Character vector. Parameters to export.
#' @param path Character. Output file path.
#' @param verbose Logical. Print progress messages. Default FALSE.
#'
#' @return Path to the created file, invisibly.
#'
#' @keywords internal
export_draws_csv <- function(fit, params, path, verbose = FALSE) {
  
  vmsg(sprintf("Extracting draws for: %s", paste(params, collapse = ", ")), verbose)
  
  if (inherits(fit, "CmdStanMCMC")) {
    draws_df <- fit$draws(params, format = "df")
  } else {
    draws_list <- rstan::extract(fit, pars = params)
    draws_df <- as.data.frame(do.call(cbind, draws_list))
  }
  
  vmsg(sprintf("Writing to %s", path), verbose)
  utils::write.csv(draws_df, path, row.names = FALSE)
  
  invisible(path)
}


#' Export model summary to text file
#'
#' Creates a plain text summary of model results.
#'
#' @param fit_res List returned by \code{\link{fit_ou_nonlinear_tmg}}
#' @param path Character. Output file path.
#' @param verbose Logical. Print progress messages. Default FALSE.
#'
#' @return Path to the created file, invisibly.
#'
#' @keywords internal
export_summary_txt <- function(fit_res, path, verbose = FALSE) {
  
  vmsg(sprintf("Creating summary at %s", path), verbose)
  
  lines <- c(
    "========================================",
    "BAYESIAN NONLINEAR OU MODEL SUMMARY",
    "========================================",
    "",
    sprintf("Model: %s", fit_res$factor_ou$model %||% "unknown"),
    sprintf("Training observations: %d", 
            fit_res$factor_ou$factor_ou_info$T_train %||% NA),
    "",
    "--- Key Parameters ---",
    sprintf("beta1 (TMG effect): %.4f", fit_res$factor_ou$beta1 %||% NA),
    sprintf("nu (df Student-t): %.2f", fit_res$factor_ou$nu %||% NA),
    sprintf("gamma (COM effect): %.4f", fit_res$factor_ou$gamma %||% NA),
    "",
    "--- Diagnostics ---",
    sprintf("Max R-hat: %.4f", fit_res$diagnostics$rhat_max %||% NA),
    sprintf("Divergences: %d", fit_res$diagnostics$divergences %||% NA),
    "",
    "--- Convergence (kappa_s) ---",
    sprintf("Mean kappa: %.4f", mean(fit_res$factor_ou$kappa_s, na.rm = TRUE)),
    sprintf("Range: [%.4f, %.4f]", 
            min(fit_res$factor_ou$kappa_s, na.rm = TRUE),
            max(fit_res$factor_ou$kappa_s, na.rm = TRUE)),
    ""
  )
  
  writeLines(lines, path)
  
  invisible(path)
}
