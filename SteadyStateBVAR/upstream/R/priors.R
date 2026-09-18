#' Specify priors for the steady-state BVAR model
#'
#' This function prepares the priors. The Minnesota prior is used for the autoregressive parameters, and is determined by
#' the overall tightness, cross-equation tightness, and the lag decay rate.
#' For the steady-state parameters, a normal prior is used. Regarding the covariance matrix of the reduced-form innovations, for the homoscedastic model,
#' the user can choose between Jeffreys prior or an uninformative inverse-Wishart prior.
#' The user may also choose a stochastic volatility specification for the innovations, where the covariance matrix varies over time.
#'
#' @param x A steady-state \code{bvar} object that has been passed through \code{\link{setup}}.
#' @param lambda_1 Numeric. Overall tightness of the Minnesota prior.
#'   Default \code{0.2}.
#' @param lambda_2 Numeric. Cross-equation tightness of the Minnesota prior. Default \code{0.5}.
#' @param lambda_3 Numeric. Lag decay rate of the Minnesota prior. Default \code{1}.
#' @param first_own_lag_prior_mean Numeric vector of length \code{k}. Prior means for the first own lags
#'   of the variables. If \code{NULL} (default), a zero vector is used.
#' @param theta_Psi Numeric vector. Prior mean vector for \eqn{\text{vec}(\Psi)}, i.e. the steady-state parameters. If \code{NULL} (default),
#'   the OLS estimates are used.
#' @param Omega_Psi Numeric matrix. Prior covariance matrix for \eqn{\text{vec}(\Psi)}, i.e. the steady-state parameters. Assumed to be diagonal. If \code{NULL} (default),
#'    the identity matrix is used.
#' @param Jeffreys Logical. If \code{TRUE} (default), uses Jeffreys prior for the innovation covariance matrix.
#'   If \code{FALSE}, uses an uninformative inverse-Wishart prior. Only considered if \code{SV=FALSE}.
#' @param SV Logical. If \code{TRUE}, enables stochastic volatility specification.
#'   Default \code{FALSE}.
#' @param SV_type Character. Type of stochastic volatility model. Must be either \code{"RW"} or \code{"AR1"}.
#'   Required if \code{SV = TRUE}.
#' @param SV_priors List. User-supplied stochastic volatility priors. Required when \code{SV = TRUE}.
#'   The list must contain the following named elements depending on \code{SV_type}:
#'   \itemize{
#'     \item For \code{"RW"}: \code{theta_A}, \code{Omega_A}, \code{theta_log_lambda_1},
#'       \code{Omega_log_lambda_1}, \code{alpha_phi}, \code{beta_phi}.
#'     \item For \code{"AR1"}: \code{theta_A}, \code{Omega_A}, \code{theta_gamma_0},
#'       \code{Omega_gamma_0}, \code{theta_gamma_1}, \code{Omega_gamma_1},
#'       \code{theta_log_lambda_1}, \code{Omega_log_lambda_1}, \code{V_Phi}, \code{m_Phi}.
#'   }
#'
#' @return The steady-state \code{bvar} object with an appended \code{priors} list containing:
#'   \item{theta_beta}{Prior mean vector for \eqn{\text{vec}(\beta)}, i.e. the autoregressive parameters, constructed with the Minnesota prior}
#'   \item{Omega_beta}{Prior covariance matrix for \eqn{\text{vec}(\beta)}, i.e. the autoregressive parameters, constructed with the Minnesota prior}
#'   \item{theta_Psi}{Prior mean vector for \eqn{\text{vec}(\Psi)}, i.e. the steady-state parameters}
#'   \item{Omega_Psi}{Prior covariance matrix for \eqn{\text{vec}(\Psi)}, i.e. the steady-state parameters}
#'   \item{Jeffreys}{Indicator for Jeffreys prior usage}
#'   \item{Sigma_AR}{Residual variance estimates from univariate AR fits, which are used by the Minnesota prior}
#'   \item{m}{Inverse-Wishart prior degrees of freedom (if \code{Jeffreys = FALSE} and \code{SV = FALSE})}
#'   \item{V}{Inverse-Wishart prior scale matrix (if \code{Jeffreys = FALSE} and \code{SV = FALSE})}
#'   \item{SV}{Logical indicator for stochastic volatility specification}
#'   \item{SV_type}{Stochastic volatility specification type}
#'   \item{SV_priors}{User-supplied SV prior list (if \code{SV = TRUE})}
#'   
#' @details
#' 
#' The goal is to estimate the parameters \eqn{\Theta = \begin{bmatrix} \beta & \Psi & \Sigma_u \end{bmatrix}},
#' and as such priors are needed. Following Villani (2009), prior independence between \eqn{\beta, \Psi} and \eqn{\Sigma_u} is assumed.
#' For \eqn{\beta}, which contains the autoregressive parameters, the Minnesota prior is used
#'
#' \deqn{\mathrm{vec}(\beta) \sim \mathrm{N}_{kpk} (\theta_\beta,\Omega_\beta)}
#'
#' The prior means (the elements of \eqn{\theta_\beta}) are specified according to
#'
#' \deqn{
#' \begin{aligned}
#' \mathrm{E}\left(\Pi_{\ell}^{(i,j)}\right)&=
#' \begin{cases}\kappa & \text{if } \ell = 1 \ \text{and} \ i = j \\0 & \text{otherwise}\end{cases}\\
#' \kappa&=\begin{cases}\kappa^{levels} & \text{if } \text{variable} \ i \ \text{is in levels} \\
#' \kappa^{\Delta} & \text{if } \text{variable} \ i \ \text{is differenced}
#' \end{cases}\\
#' \end{aligned}
#' }
#'
#' Here, the autoregressive coefficient \eqn{\Pi_{\ell}^{(i,j)}} is element \eqn{\left(i,j\right)}
#' of \eqn{\Pi_{\ell}} for \eqn{\ell=1,\dots,p}. As such, the Minnesota prior sets all prior means
#' for the elements in \eqn{\beta} to \eqn{0}, except for the elements that relate to the first
#' own lags of the variables, which are set to \eqn{\kappa}. If variable \eqn{i} is in levels
#' (e.g. nominal interest rate), then \eqn{\kappa=\kappa^{levels}}, and typical choices for
#' \eqn{\kappa^{levels}} are \eqn{1} or \eqn{0.9}. Evaluating the equations at their prior means,
#' equation \eqn{i} becomes a random walk if \eqn{\kappa^{levels}=1} and a persistent stationary
#' AR(1) process if \eqn{\kappa^{levels}=0.9}. Since the steady state only exists if the process
#' is stationary, \eqn{0.9} is recommended for the steady-state BVAR. If variable \eqn{i} is differenced
#' (e.g. GDP growth), then \eqn{\kappa=\kappa^{\Delta}}, and the most common choice for
#' \eqn{\kappa^{\Delta}} is \eqn{0}, i.e. equation \eqn{i} becomes (when evaluating it at its prior means)
#' a random walk expressed in first differences (a white noise process). If a differenced variable still shows some degree
#' of persistence (can be examined with an ACF plot), a suitable value for \eqn{\kappa^{\Delta}} can
#' be (for example) \eqn{0.6} instead of \eqn{0}. Moving on to the prior variances, \eqn{\Omega_\beta} is a
#' diagonal matrix containing the prior variances for the elements in \eqn{\beta}. They are specified as
#' 
#' \deqn{\mathrm{Var}\left(\Pi_{\ell}^{(i,j)}\right)=
#' \begin{cases}\left(\frac{\lambda_1}{\ell^{\lambda_3}}\right)^2 & \text{if } i = j \\
#' \left(\frac{\lambda_1 \lambda_2\sigma_i}{\ell^{\lambda_3}\sigma_j}\right)^2& \text{if } i \neq j
#' \end{cases}}
#' 
#' Here \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\lambda_3} are scalar hyperparameters known as
#' the overall tightness, the cross-equation tightness and the lag decay rate. Furthermore,
#' \eqn{\sigma_i^2} is the \eqn{(i,i)}:th element of \eqn{\Sigma_u}, which is unknown and therefore
#' needs to be replaced with an estimate. In this package, it is replaced by the least squares residual variance
#' from a univariate autoregression for variable \eqn{i} with \eqn{p} lags
#' (including the constant and dummy/trend variable if applicable). Moving on to \eqn{\Psi}, the steady-state parameter matrix, the
#' (steady-state) prior is
#' 
#' \deqn{\mathrm{vec}(\Psi) \sim \mathrm{N}_{kq}(\theta_\Psi,\Omega_\Psi)}
#' 
#' This is the core of the steady-state BVAR model.
#' In \eqn{\theta_\Psi}, one specifies the prior beliefs about the location of the steady state,
#' and in \eqn{\Omega_\Psi}, which is assumed to be a diagonal matrix, one specifies the degree
#' of certainty in those prior beliefs. Too see how to specify steady-state priors in practice, please see the package vignettes.
#' The prior for \eqn{\Sigma_u} is either the usual non-informative Jeffreys prior
#' 
#' \deqn{p(\Sigma_u) \propto\left|\Sigma_u \right|^{-(k+1)/2}}
#' 
#' or a proper uninformative inverse-Wishart prior
#' 
#' \deqn{\Sigma_u \sim \mathrm{IW}(V,m)}
#' 
#' where \eqn{V} is the scale matrix and \eqn{m} is the number of degrees of freedom.
#' An uninformative prior is specified by setting
#' \eqn{V=(m-k-1)\hat{\Sigma}_u} where \eqn{\hat{\Sigma}_u} is the least squares estimate
#' from the VAR(\eqn{p}) (including the constant and dummy/trend variable if applicable), and \eqn{m=k+2}.
#' However, if we opt to use stochastic volatility, the innovation covariance matrix is time-varying \eqn{\Sigma_{u,t}},
#' instead of constant \eqn{\Sigma_{u}}. Therefore, stochastic volatility priors are needed.
#' See \link{bvar} for more details on the stochastic volatility specifications. Please note that
#' \eqn{\lambda} below (volatilities) has nothing to do with the \eqn{\lambda} from the Minnesota prior (hyperparameters).
#' Now, for the Random Walk (\code{"RW"}) stochastic volatility specification, the following priors are available
#' 
#' \deqn{\begin{aligned}a &\sim \mathrm{N}(\theta_A, \Omega_A) \\
#' \ln \lambda_{1} &\sim \mathrm{N}(\theta_{\ln \lambda_{1}}, \Omega_{\ln \lambda_{1}}) \\
#' \phi_i &\sim \mathrm{IG}(\alpha_{\phi_i},\beta_{\phi_i})\end{aligned}}
#' 
#' Here \eqn{a} is a \eqn{k(k-1)/2}-dimensional vector that collects the free parameters in \eqn{A} in row-major order,
#' and \eqn{\ln \lambda_1} is a \eqn{k}-dimensional vector containing the time \eqn{t=1} values (initial conditions) of \eqn{\ln \lambda_{t}}.
#' We assume that \eqn{\Omega_A} and \eqn{\Omega_{\ln \lambda_{1}}} are diagonal matrices.
#' Furthermore, \eqn{\phi_i} for \eqn{i=1,\dots,k} are the log volatility innovation variances. For the AR(1) (\code{"AR1"}) stochastic volatility specification, the following priors are available
#' 
#' \deqn{\begin{aligned}a &\sim \mathrm{N}(\theta_A, \Omega_A) \\
#' \gamma_{0} &\sim \mathrm{N}(\theta_{\gamma_0}, \Omega_{\gamma_0}) \\
#' \gamma_{1} &\sim \mathrm{N}(\theta_{\gamma_1}, \Omega_{\gamma_1}) \ I(|\gamma_{1,i}| < 1)\\
#' \ln \lambda_{1} &\sim \mathrm{N}(\theta_{\ln \lambda_{1}}, \Omega_{\ln \lambda_{1}}) \\
#' \Phi &\sim \mathrm{IW}(V_{\Phi},m_{\Phi})\end{aligned}}
#' 
#' Here \eqn{a} is again the \eqn{k(k-1)/2}-dimensional vector that collects the free parameters in \eqn{A} in row-major order,
#' and \eqn{\ln \lambda_1} is a \eqn{k}-dimensional vector containing the time \eqn{t=1} values (initial conditions) of \eqn{\ln \lambda_{t}}.
#' Furthermore, \eqn{\gamma_{0}} is a \eqn{k}-dimensional vector of log volatility intercepts, \eqn{\gamma_{1}} is a \eqn{k}-dimensional vector of log volatility
#' slopes, and \eqn{\Phi} is the \eqn{k \times k} log volatility innovation covariance matrix.
#' We assume that \eqn{\Omega_A}, \eqn{\Omega_{\gamma_0}}, \eqn{\Omega_{\gamma_1}}, and \eqn{\Omega_{\ln \lambda_{1}}} are diagonal matrices.
#' Note that the prior for \eqn{\gamma_1} is truncated normal.
#' 
#' For details on the homoscedastic steady-state BVAR model, see Villani (2009).
#' For details on the Random Walk stochastic volatility steady-state BVAR model, see Clark (2011).
#' See Carriero, Clark, and Marcellino (2024) for the AR(1) stochastic volatility
#' specification applied to a conventional BVAR.
#' 
#' To see examples of how to specify the priors, simply view the relevant vignette
#' 
#' \itemize{
#'   \item \code{vignette("Homoscedastic-steady-state-BVAR")}
#'   \item \code{vignette("RW-stochastic-volatility-steady-state-BVAR")}
#'   \item \code{vignette("AR1-stochastic-volatility-steady-state-BVAR")}
#' }
#' 
#' @references
#' Carriero, A., Clark, T. E., and Marcellino, M. (2024).
#' Capturing macro-economic tail risks with Bayesian vector autoregressions. 
#' \emph{Journal of Money, Credit and Banking}, 56(5), pp. 1099–1127.
#' 
#' Clark, T. E. (2011). Real-time density forecasts from Bayesian vector autoregressions
#' with stochastic volatility. \emph{Journal of Business & Economic Statistics}, 29(3), pp. 327–341.
#' 
#' Villani, M. (2009). Steady-state priors for vector autoregressions.
#' \emph{Journal of Applied Econometrics}, 24(4), pp. 630–650. 
#' 
#' @export
#'
#' @examples
#' #homoscedastic with Jeffreys prior
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
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
#' #RW stochastic volatility
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#' 
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#' 
#' SV_priors_RW <- list(
#' theta_A              =  rep(0, n_free_params_A),
#' Omega_A              =  diag(1000, n_free_params_A),
#' theta_log_lambda_1   =  rep(0, k),
#' Omega_log_lambda_1   =  diag(1000, k),
#' alpha_phi            =  rep(5, k),
#' beta_phi             = (rep(5, k) - 1) * rep(0.1, k)
#' )
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    SV = TRUE,
#'                    SV_type = "RW",
#'                    SV_priors = SV_priors_RW)
#'                    
#' #AR1 stochastic volatility
#' yt <- matrix(rnorm(50), 25, 2)
#'
#' bvar_obj <- bvar(data = yt)
#'
#' bvar_obj <- setup(bvar_obj, p=1)
#' 
#' k <- bvar_obj$setup$k
#' n_free_params_A <- bvar_obj$setup$n_free_params_A
#' 
#' SV_priors_AR1 <- list(
#' theta_A               =  rep(0, n_free_params_A),
#' Omega_A               =  diag(1000, n_free_params_A),
#' theta_gamma_0         =  rep(0.1, k),
#' Omega_gamma_0         =  diag(1000, k),
#' theta_gamma_1         =  rep(0.9, k),
#' Omega_gamma_1         =  diag(10, k),
#' theta_log_lambda_1    =  rep(0.1, k)/(1-rep(0.9, k)),
#' Omega_log_lambda_1    =  diag(1000, k),
#' V_Phi                 = (10 - k - 1) * diag(k),
#' m_Phi                 =  10
#' )
#'
#' bvar_obj <- priors(bvar_obj,
#'                    lambda_1 = 0.2,
#'                    lambda_2 = 0.5,
#'                    lambda_3 = 1,
#'                    first_own_lag_prior_mean = rep(1,2),
#'                    theta_Psi = rep(0, 2),
#'                    Omega_Psi = diag(0.1, 2, 2),
#'                    SV = TRUE,
#'                    SV_type = "AR1",
#'                    SV_priors = SV_priors_AR1)
#'
priors<- function(x,
                  lambda_1=0.2,
                  lambda_2=0.5,
                  lambda_3 = 1,
                  first_own_lag_prior_mean=NULL,
                  theta_Psi=NULL,
                  Omega_Psi=NULL,
                  Jeffreys=TRUE,
                  SV = FALSE,
                  SV_type = NULL,
                  SV_priors = NULL){
  
  if (!inherits(x, "bvar")) stop("must be a 'bvar' object")
  if (is.null(x$setup)) stop("must be passed through setup")
  
  setup <- x$setup
  
  if (!all(c("k","p","q","dt") %in% names(setup))) {
    stop("setup is incomplete")
  }
  
  k <- setup$k
  p <- setup$p
  q <- setup$q
  dt <- setup$dt
  yt <- x$data
  
  if (!is.null(theta_Psi) && !is.null(Omega_Psi)) {
    if (nrow(Omega_Psi) != length(theta_Psi)) {
      stop("dimensions must match")
    }
  }
  
  if (any(c(lambda_1, lambda_2, lambda_3) <= 0)) {
    stop("must be positive")
  }
  
  if (isTRUE(SV)) {
    if (is.null(SV_type)) stop("SV_type is needed")
    if (!SV_type %in% c("RW", "AR1")) stop("SV_type needs to be RW or AR1")
    if (is.null(SV_priors)) stop("SV_priors is required")
    
    n_free <- x$setup$n_free_params_A
    
    if (SV_type == "RW") {
      required_names <- c("theta_A", "Omega_A", "theta_log_lambda_1",
                          "Omega_log_lambda_1", "alpha_phi", "beta_phi")
      missing_names <- setdiff(required_names, names(SV_priors))
      if (length(missing_names) > 0)
        stop(paste("SV_priors is missing elements:", paste(missing_names, collapse = ", ")))
      
      if (length(SV_priors$theta_A) != n_free)
        stop(paste("theta_A must be a vector of length k*(k-1)/2 =", n_free))
      if (!all(dim(SV_priors$Omega_A) == n_free))
        stop(paste("Omega_A must be a", n_free, "x", n_free, "matrix"))
      if (length(SV_priors$theta_log_lambda_1) != k)
        stop(paste("theta_log_lambda_1 must be a vector of length k =", k))
      if (!all(dim(SV_priors$Omega_log_lambda_1) == k))
        stop(paste("Omega_log_lambda_1 must be a", k, "x", k, "matrix"))
      if (length(SV_priors$alpha_phi) != k)
        stop(paste("alpha_phi must be a vector of length k =", k))
      if (any(SV_priors$alpha_phi <= 0))
        stop("alpha_phi must be strictly positive")
      if (length(SV_priors$beta_phi) != k)
        stop(paste("beta_phi must be a vector of length k =", k))
      if (any(SV_priors$beta_phi <= 0))
        stop("beta_phi must be strictly positive")
      
    } else if (SV_type == "AR1") {
      required_names <- c("theta_A", "Omega_A", "theta_gamma_0", "Omega_gamma_0",
                          "theta_gamma_1", "Omega_gamma_1", "theta_log_lambda_1",
                          "Omega_log_lambda_1", "m_Phi", "V_Phi")
      missing_names <- setdiff(required_names, names(SV_priors))
      if (length(missing_names) > 0)
        stop(paste("SV_priors is missing elements:", paste(missing_names, collapse = ", ")))
      
      if (length(SV_priors$theta_A) != n_free)
        stop(paste("theta_A must be a vector of length k*(k-1)/2 =", n_free))
      if (!all(dim(SV_priors$Omega_A) == n_free))
        stop(paste("Omega_A must be a", n_free, "x", n_free, "matrix"))
      if (length(SV_priors$theta_gamma_0) != k)
        stop(paste("theta_gamma_0 must be a vector of length k =", k))
      if (!all(dim(SV_priors$Omega_gamma_0) == k))
        stop(paste("Omega_gamma_0 must be a", k, "x", k, "matrix"))
      if (length(SV_priors$theta_gamma_1) != k)
        stop(paste("theta_gamma_1 must be a vector of length k =", k))
      if (!all(dim(SV_priors$Omega_gamma_1) == k))
        stop(paste("Omega_gamma_1 must be a", k, "x", k, "matrix"))
      if (length(SV_priors$theta_log_lambda_1) != k)
        stop(paste("theta_log_lambda_1 must be a vector of length k =", k))
      if (!all(dim(SV_priors$Omega_log_lambda_1) == k))
        stop(paste("Omega_log_lambda_1 must be a", k, "x", k, "matrix"))
      if (!is.numeric(SV_priors$m_Phi) || length(SV_priors$m_Phi) != 1 || SV_priors$m_Phi < k)
        stop(paste("m_Phi must be a scalar integer >= k =", k))
      if (!all(dim(SV_priors$V_Phi) == k))
        stop(paste("V_Phi must be a", k, "x", k, "matrix"))
    }
  }
  
  priors <- list()
  
  Sigma_AR <- x$setup$Sigma_AR
  sigma <- sqrt(diag(Sigma_AR))
  
  Vv <- lapply(1:p, function(x) matrix(0, k, k))
  
  for (l in 1:p) {
    for (i in 1:k) {
      for (j in 1:k) {
        if (i == j) {
          Vv[[l]][i,j] <- (lambda_1/(l^lambda_3))^2
        } else {
          Vv[[l]][i,j] <- ((lambda_1*lambda_2*sigma[i])/(l^lambda_3*sigma[j]))^2
        }
      }
    }
  }
  
  V_mat <- do.call(cbind, Vv)
  Omega_beta <- diag(c(t(V_mat)))
  
  if (!is.null(first_own_lag_prior_mean)) {
    if (length(first_own_lag_prior_mean) != k) {
      stop("must have length k")
    }
  }
  
  if (is.null(first_own_lag_prior_mean)) first_own_lag_prior_mean <- rep(0,k)
  
  mat <- matrix(0, nrow = k*p, ncol = k)
  for (i in 1:k){
    mat[i,i] <- first_own_lag_prior_mean[i]
  }
  theta_beta = c(mat)
  
  if (is.null(theta_Psi)) theta_Psi <- c(x$setup$Psi_OLS)
  if (is.null(Omega_Psi)) Omega_Psi <- diag(1, k*q, k*q)
  
  if (isFALSE(Jeffreys) && isFALSE(SV)){
    m=k+2
    V = (m-k-1)*setup$Sigma_u_OLS
    priors$V <- V
    priors$m <- m
  }
  
  priors$theta_beta <- theta_beta
  priors$Omega_beta <- Omega_beta
  priors$theta_Psi <- theta_Psi
  priors$Omega_Psi <- Omega_Psi
  priors$Jeffreys <- Jeffreys
  priors$Sigma_AR <- Sigma_AR
  
  priors$SV <- SV
  priors$SV_type <- SV_type
  
  if (isTRUE(SV)) {
    priors$SV_priors <- SV_priors
    priors$Jeffreys <- FALSE
  }
  
  x$priors <- priors
  return(x)
}