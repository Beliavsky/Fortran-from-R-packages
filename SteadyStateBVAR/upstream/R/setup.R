#' Set up the steady-state BVAR model
#'
#' Prepares the components needed for prior specification and estimation. 
#' Also computes OLS estimates.
#'
#' @param x A steady-state \code{bvar} object created by \code{\link{bvar}}.
#' @param p Integer. The lag order of the VAR. Default \code{1}.
#' @param deterministic Character. The deterministic components to include in \eqn{d_t}.
#'   One of \code{"constant"} (default), \code{"constant_and_dummy"}, or
#'   \code{"constant_and_trend"}.
#' @param dummy Numeric vector of a dummy variable. Only
#'   used when \code{deterministic = "constant_and_dummy"}. Default \code{NULL}.
#'
#' @return The steady-state \code{bvar} object with a \code{setup} list containing the
#' required components for prior specification and estimation, and also the OLS estimates.
#' @export
#'
#' @examples
#' yt <- matrix(rnorm(50), 25, 2)
#' 
#' bvar_obj <- bvar(data = yt)
#' 
#' bvar_obj <- setup(bvar_obj)
setup <- function(x, p = 1, deterministic=c("constant", "constant_and_dummy", "constant_and_trend"), dummy=NULL) {
  
  if (!inherits(x, "bvar")) {
    stop("x must be a 'bvar' object")
  }
  
  if (is.null(x$data)) {
    stop("x must have data. Create bvar object with bvar(data = ...)")
  }
  
  if (!is.numeric(p) || p <= 0 || p != floor(p)) {
    stop("p must be a positive integer")
  }
  
  yt <- x$data
  
  if (p >= nrow(yt)) {
    stop("p must be less than the number of observations")
  }
  
  deterministic <- match.arg(deterministic)
  
  if (deterministic == "constant_and_dummy" && is.null(dummy)) {
    stop("dummy must be provided when deterministic = 'constant_and_dummy'")
  }
  
  if (!is.null(dummy) && length(dummy) != nrow(yt)) {
    stop("dummy must have length equal to number of observations")
  }
  
  if (!is.null(dummy) && !is.vector(dummy)) {
    stop("dummy must be a vector")
  }
  
  
  N = nrow(yt) - p
  k = ncol(yt)
  
  if (deterministic == "constant") {
    dt <- cbind(rep(1, nrow(yt)))
    q <- 1
  } else if (deterministic == "constant_and_dummy") {
    dt <- cbind(rep(1, nrow(yt)), dummy)
    q <- 2
  } else {
    trend <- 1:nrow(yt)
    dt <- cbind(rep(1, nrow(yt)), trend)
    q <- 2
  }
  
  Y <- yt[-c(1:p), ]
  W <- embed(yt, dimension = p+1)[, -(1:k)]
  X <- dt[-c(1:p), , drop=FALSE]
  Q <- embed(dt, dimension = p+1)[, -(1:q), drop=FALSE]
  
  Z <- cbind(W, X)
  beta_hat = solve(crossprod(Z), crossprod(Z, Y))
  U = Y - Z %*% beta_hat
  Sigma_u_OLS <- crossprod(U) / (N - k*p - q)
  
  if (q == 1) {
    C_hat <- beta_hat[(k*p+1):(k*p+q), ]
  } else {
    C_hat <- t(beta_hat[(k*p+1):(k*p+q), ])
  }
  
  A <- vector("list", p)
  for (i in 1:p) {
    rows_idx <- ((i - 1) * k + 1):(i * k)
    A[[i]] <- matrix(t(beta_hat[rows_idx, ]), nrow = k, ncol = k)
  }
  
  A_L <- diag(k)
  for (i in 1:p) {
    A_L <- A_L - A[[i]]
  }
  
  Psi_OLS <- solve(A_L, C_hat)
  beta_OLS = beta_hat[1:(k*p), ]
  
  n_free_params_A <- k*(k-1)/2
  
  Sigma_AR <- diag(0, k)
  
  for (i in 1:k){
    
    y <- yt[,i]
    
    Y_AR <- y[-c(1:p)]
    W_AR <- embed(y, dimension = p+1)[, -1]
    X_AR <- dt[-c(1:p), ,drop=FALSE]
    
    Z_AR <- cbind(W_AR,X_AR)
    beta_hat_AR = solve(crossprod(Z_AR,Z_AR),crossprod(Z_AR,Y_AR))
    U_AR = Y_AR-Z_AR%*%beta_hat_AR
    sigma2 <- crossprod(U_AR,U_AR)/(nrow(Z_AR)-ncol(Z_AR))
    Sigma_AR[i,i] <- sigma2
  }
  
  x$setup <- list(N=N,
                  k=k,
                  p=p,
                  Y=Y,
                  X=X,
                  W=W,
                  Q=Q,
                  q=q,
                  dummy=dummy,
                  beta_OLS=beta_OLS,
                  Sigma_u_OLS=Sigma_u_OLS,
                  Psi_OLS=Psi_OLS,
                  dt=dt,
                  D=X,
                  n_free_params_A=n_free_params_A,
                  Sigma_AR = Sigma_AR,
                  deterministic=deterministic)
  return(x)
}