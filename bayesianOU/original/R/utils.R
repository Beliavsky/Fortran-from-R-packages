#' Null-coalescing operator
#'
#' Returns the left operand if not NULL, otherwise the right operand.
#'
#' @param x Left operand
#' @param y Right operand (default value)
#'
#' @return x if not NULL, otherwise y
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Z-score standardization using training period statistics
#'
#' Standardizes a matrix using mean and standard deviation computed
#' from the training period only.
#'
#' @param M Numeric matrix of dimensions T x S (time by sectors)
#' @param T_train Integer. Number of observations in training period
#' @param eps Numeric. Minimum standard deviation threshold to avoid
#'   division by zero. Default is 1e-8.
#'
#' @return A list with components:
#'   \describe{
#'     \item{Mz}{Standardized matrix of same dimensions as M}
#'     \item{mu}{Vector of training means for each column}
#'     \item{sd}{Vector of training standard deviations for each column}
#'   }
#'
#' @examples
#' M <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' result <- zscore_train(M, T_train = 14)
#' str(result)
#'
#' @export
zscore_train <- function(M, T_train, eps = 1e-8) {
  if (!is.matrix(M)) {
    M <- as.matrix(M)
  }

  if (!is.numeric(T_train) || length(T_train) != 1L ||
      !is.finite(T_train)) {
    stop("`T_train` must be a single finite number.", call. = FALSE)
  }
  T_train <- as.integer(T_train)
  if (T_train < 2L) {
    stop("`T_train` must be at least 2 to estimate a standard deviation.",
         call. = FALSE)
  }
  if (T_train > nrow(M)) {
    stop(sprintf("`T_train` (%d) exceeds the number of rows (%d).",
                 T_train, nrow(M)), call. = FALSE)
  }

  train_data <- M[seq_len(T_train), , drop = FALSE]
  mu <- colMeans(train_data, na.rm = TRUE)
  sdv <- apply(train_data, 2, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv < eps] <- 1
  
  Mz <- sweep(M, 2, mu, `-`)
  Mz <- sweep(Mz, 2, sdv, `/`)
  
  list(Mz = Mz, mu = mu, sd = sdv)
}


#' Align matrix columns to reference
#'
#' Reorders and filters columns of matrix B to match column names of A.
#'
#' @param A Reference matrix with target column names
#' @param B Matrix to align
#' @param verbose Logical. Print alignment messages. Default FALSE.
#'
#' @return Matrix B with columns reordered to match A
#'
#' @keywords internal
align_columns <- function(A, B, verbose = FALSE) {
  if (!is.null(colnames(A)) && !is.null(colnames(B))) {
    missing_in_B <- setdiff(colnames(A), colnames(B))
    extra_in_B <- setdiff(colnames(B), colnames(A))
    
    if (length(missing_in_B) > 0) {
      stop(
        "Missing columns in B: ",
        paste(missing_in_B, collapse = ", ")
      )
    }
    
    if (length(extra_in_B) > 0 && verbose) {
      message(
        "Extra columns in B will be ignored: ",
        paste(extra_in_B, collapse = ", ")
      )
    }
    
    B[, colnames(A), drop = FALSE]
  } else {
    B
  }
}


#' Check availability of Stan backend
#'
#' Verifies if cmdstanr or rstan is available for model fitting.
#'
#' @param verbose Logical. Print status messages. Default FALSE.
#'
#' @return Character string: "cmdstanr", "rstan", or "none"
#'
#' @keywords internal
check_stan_backend <- function(verbose = FALSE) {
  have_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE)
  
  if (have_cmdstan) {
    cmdstan_ok <- tryCatch(
      !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
      error = function(e) FALSE
    )
    if (cmdstan_ok) {
      if (verbose) message("Using cmdstanr backend")
      return("cmdstanr")
    }
  }
  
  if (requireNamespace("rstan", quietly = TRUE)) {
    if (verbose) message("Using rstan backend")
    return("rstan")
  }
  
  "none"
}


#' Verbose message helper
#'
#' Prints a message only if verbose is TRUE.
#'
#' @param msg Character. Message to print.
#' @param verbose Logical. Whether to print.
#'
#' @return NULL invisibly
#'
#' @keywords internal
vmsg <- function(msg, verbose) {
  if (isTRUE(verbose)) {
    message(msg)
  }
  invisible(NULL)
}


#' Compute common factor from matrix
#'
#' Extracts first principal component and standardizes using training period.
#'
#' @param Mz Standardized matrix (T x S)
#' @param T_train Number of training observations
#' @param use_train_loadings Logical. If TRUE (default), compute the factor
#'   loadings using only the training period and then project the full series
#'   onto them. This avoids look-ahead leakage. If FALSE, the loadings are
#'   computed from the full sample (SVD over all T), which lets future
#'   information influence the constructed factor; only use FALSE for purely
#'   in-sample, descriptive analyses.
#' @param verbose Logical. Print progress messages. Default FALSE.
#'
#' @return Numeric vector of factor scores (length T)
#'
#' @keywords internal
compute_common_factor <- function(Mz, T_train, use_train_loadings = TRUE,
                                  verbose = FALSE) {
  if (isTRUE(use_train_loadings)) {
    vmsg("Computing factor using training-period loadings only", verbose)
    M_train <- Mz[seq_len(T_train), , drop = FALSE]
    
    pc <- tryCatch(
      stats::prcomp(M_train, center = FALSE, scale. = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(pc)) {
      sv <- base::svd(M_train)
      load1 <- sv$v[, 1, drop = FALSE]
    } else {
      load1 <- pc$rotation[, 1, drop = FALSE]
    }
    
    Ft_all <- as.numeric(Mz %*% load1)
    muF <- mean(Ft_all[seq_len(T_train)])
    sdF <- stats::sd(Ft_all[seq_len(T_train)])
    if (!is.finite(sdF) || sdF < 1e-8) sdF <- 1
    
    Ft <- (Ft_all - muF) / sdF
  } else {
    vmsg("Computing factor using full-sample SVD", verbose)
    M_centered <- scale(Mz, center = TRUE, scale = FALSE)
    sv <- base::svd(M_centered)
    Ft <- sv$u[, 1]
    
    muF <- mean(Ft[seq_len(T_train)])
    sdF <- stats::sd(Ft[seq_len(T_train)])
    Ft <- (Ft - muF) / sdF
  }
  
  as.numeric(Ft)
}
