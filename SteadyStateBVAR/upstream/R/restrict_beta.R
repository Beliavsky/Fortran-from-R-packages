#' Restrict VAR coefficients to zero
#'
#' Applies zero restrictions to the VAR coefficient matrix \eqn{\beta} by setting the
#' corresponding prior variances in \code{Omega_beta} to a value near zero (prior means are zero by default in the Minnesota prior).
#' 
#' @details
#' The steady-state BVAR model takes the form
#'
#' \deqn{y_t = \Psi d_t + \Pi_1(y_{t-1}-\Psi d_{t-1})+\dots+\Pi_p(y_{t-p}-\Psi d_{t-p})+u_t}
#'
#' where \eqn{y_t} is an \eqn{k}-dimensional vector of endogenous variables at time \eqn{t}, and \eqn{d_t} is
#' a \eqn{q}-dimensional vector of deterministic (exogenous) variables at time \eqn{t}.
#' Here \eqn{\Pi_\ell} for \eqn{\ell=1,\dots,p} is a \eqn{(k \times k)} matrix of autoregressive parameters,
#' and \eqn{\Psi} is a \eqn{(k \times q)} matrix of steady-state parameters.
#' One can stack the (transposed) \eqn{\Pi_i} matrices in the \eqn{(kp \times k)} matrix \eqn{\beta}
#' \deqn{\beta=\begin{bmatrix}\Pi'_1 \\ \vdots  \\\Pi'_p\end{bmatrix}}
#' 
#' This function puts zero restrictions on the autoregressive parameters/coefficients, i.e the elements of \eqn{\beta}.
#'
#' @param x A steady-state \code{bvar} object that has been passed through \code{\link{priors}}.
#' @param restriction_matrix A numeric matrix of dimension \code{k*p x k} where
#'   entries of \code{0} indicate coefficients to be restricted to zero and
#'   entries of \code{1} indicate unrestricted coefficients.
#'
#' @return The steady-state \code{bvar} object with the restriction matrix stored in
#'   \code{setup} and the prior covariance matrix \code{Omega_beta} updated accordingly.
#' @export
#'
#' @examples
#' yt <- matrix(rnorm(50), 25, 2)
#' colnames(yt) <- c("y1", "y2")
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=2)
#'
#' bvar_obj <- priors(bvar_obj,
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2))
#'                    
#' p <- bvar_obj$setup$p
#' k <- bvar_obj$setup$k
#' 
#' restriction_matrix <- matrix(1,k*p,k)
#' 
#' #restrict beta so y2 does not granger-cause y1
#' 
#' restriction_matrix[2, 1] <- 0
#' restriction_matrix[4, 1] <- 0
#' 
#' bvar_obj <- restrict_beta(bvar_obj, restriction_matrix)
restrict_beta <- function(x, restriction_matrix) {
  
  k <- x$setup$k
  p <- x$setup$p
  
  if (!all(dim(restriction_matrix) == c(k * p, k))) {
    stop("restriction_matrix must have dimension (k*p x k)")
  }
  
  var_names <- colnames(x$data)
  if (is.null(var_names)) {
    var_names <- paste0("Var", 1:ncol(x$data))
  }
  
  var_rep   <- rep(var_names, p)
  lag_names <- paste0("l", rep(1:p, each = k))
  
  rownames(restriction_matrix) <- paste0(var_rep, ".", lag_names)
  colnames(restriction_matrix) <- var_names
  
  x$setup$restriction_matrix <- restriction_matrix
  zero_indices <- which(c(restriction_matrix) == 0)
  
  if (!is.null(x$priors$Omega_beta)) {
    diag(x$priors$Omega_beta)[zero_indices] <- 0.0000001
    cat("Restrictions applied using restriction matrix:\n\n",sep="")
    print(restriction_matrix)
    cat("\n1 indicates that the parameter is free\n0 indicates that the parameter is restricted to zero\n")
  } else {
    warning("Omega_beta not found in priors: restriction applied but Omega_beta not updated")
  }
  
  return(x)
}