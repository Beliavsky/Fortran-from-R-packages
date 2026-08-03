#' Estimates a standard VAR(1)
#'
#' @param RiskFactors numeric matrix (\code{K x Td}). Time series of risk factors.
#' @param VARtype character. Permissible choices: "unconstrained" or "constrained".
#' @param Bcon_Mat matrix (\code{K x K + 1}). Constraints matrix (includes intercept). Entries containing NAs are treated as free parameters. Default is NULL.
#'
#' @return list. Contains:
#'   \itemize{
#'   \item intercept (K x 1)
#'   \item feedback matrix (K x K)
#'   \item variance-covariance matrix (K x K) of a VAR(1)
#' }
#'
#' @section General Notation:
#' \itemize{
#'   \item \code{Td}: model time series dimension
#'   \item \code{N}: number of country-specific spanned factors
#'   \item \code{K}: total number of risk factors
#' }
#'
#' @examples
#' data(RiskFacFull)
#' # Example 1: unconstrained case
#' VAR_para1 <- VAR(RiskFacFull, VARtype = "unconstrained")
#'
#' # Example 2: constrained case
#' K <- nrow(RiskFacFull)
#' Bcon_Mat <- matrix(0, nrow = K, ncol = K + 1)
#' Bcon_Mat[, 1:3] <- NaN
#' VAR_para2 <- VAR(RiskFacFull, VARtype = "constrained", Bcon_Mat)
#'
#' @export

VAR <- function(RiskFactors, VARtype, Bcon_Mat = NULL) {
  K <- nrow(RiskFactors)
  T_dim <- ncol(RiskFactors)
  LHS <- RiskFactors[, 2:T_dim]
  RHS <- RiskFactors[, 1:(T_dim - 1)]

  if (VARtype == "unconstrained") {
    RegVAR <- stats::lm(t(LHS) ~ t(RHS)) # VAR(1) under the P.
    K0Z <- t(t(RegVAR$coefficients[1, ]))
    K1Z <- t(RegVAR$coefficients[2:(K + 1), ])
    eZ <- RegVAR$residuals
    SSZ <- crossprod(eZ) / (T_dim - 1)
  } else { # i.e. if VARtype == 'constrained'
    intercept <- rep(1, times = T_dim - 1)
    RHS <- rbind(intercept, RHS)
    Coeff <- Est_RestOLS(LHS, RHS, Bcon_Mat)
    colnames(Coeff) <- rownames(RHS)
    rownames(Coeff) <- rownames(RHS)[-1]

    K0Z <- as.matrix(Coeff[, 1])
    K1Z <- Coeff[, 2:(K + 1)]
    eZ <- LHS - Coeff %*% RHS
    SSZ <- crossprod(t(eZ)) / (T_dim - 1)
  }

  return(list(K0Z = K0Z, K1Z = K1Z, SSZ = SSZ))
}

#########################################################################################################
#' Estimate a restricted OLS model
#'
#' @param LHS left hand side variables (M x T).
#' @param RHS right hand side variables (should include the intercept, if desired) (N x T).
#' @param Rmat matrix of constraints (M x N). Entries containing NAs are treated as free parameters.
#'
#' @return matrix of coefficient (M x N)
#' @keywords internal

Est_RestOLS <- function(LHS, RHS, Rmat) {
  T_dim <- ncol(LHS)

  # Identify constrained vs free parameters
  idx_FreePara <- is.nan(Rmat) # TRUE = free parameter
  Betas <- Rmat
  Betas[is.nan(Betas)] <- 0

  # Precompute common matrices
  XX <- tcrossprod(RHS) / T_dim
  YX <- LHS %*% t(RHS) / T_dim
  X_lg <- kronecker(XX, diag(nrow(LHS)))
  Y_lg <- as.vector(YX)

  # Solve only for free parameters
  Betas[idx_FreePara] <- solve(X_lg[idx_FreePara, idx_FreePara, drop = FALSE], Y_lg[idx_FreePara], tol = 1e-25)

  return(Betas)
}
