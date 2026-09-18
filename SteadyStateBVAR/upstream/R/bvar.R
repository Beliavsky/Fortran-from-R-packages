#' Create a steady-state BVAR model object
#'
#' Initialises a steady-state \code{bvar} object. This is the starting point for all models in \code{SteadyStateBVAR}.
#' After creation, pass the object sequentially to \code{\link{setup}}, \code{\link{priors}},
#' and \code{\link{fit}} to build and estimate the model.
#' 
#' For a detailed theoretical introduction to the steady-state BVAR model, run the following
#' 
#' \code{vignette("SteadyStateBVAR-intro")}
#' 
#' @details
#' The steady-state BVAR model takes the form
#'
#' \deqn{y_t = \Psi d_t + \Pi_1(y_{t-1}-\Psi d_{t-1})+\dots+\Pi_p(y_{t-p}-\Psi d_{t-p})+u_t}
#'
#' where \eqn{y_t} is an \eqn{k}-dimensional vector of endogenous variables at time \eqn{t}, and \eqn{d_t} is
#' a \eqn{q}-dimensional vector of deterministic (exogenous) variables at time \eqn{t}.
#' Here \eqn{\Pi_\ell} for \eqn{\ell=1,\dots,p} is a \eqn{(k \times k)} matrix of autoregressive parameters,
#' and \eqn{\Psi} is a \eqn{(k \times q)} matrix of steady-state parameters. Now
#' 
#' \deqn{\mathrm{E}(y_t)=\mu_t=\Psi d_t}
#' 
#' is the unconditional mean, or the \strong{steady state} of the process.
#' Note that the current version of this package only allows for \eqn{d_t} to contain either a constant,
#' a constant and a dummy variable, or a constant and a time trend.
#' We may stack the (transposed) \eqn{\Pi_i} matrices in the \eqn{(kp \times k)} matrix \eqn{\beta}
#' \deqn{\beta=\begin{bmatrix}\Pi'_1 \\ \vdots  \\\Pi'_p\end{bmatrix}}
#' Then the model can be rewritten as a nonlinear regression (Karlsson, 2013)
#' \deqn{y_t' =d_t'\Psi' + \left[w_t'-q_t'(I_p \otimes \Psi') \right]\beta +u_t'}
#' where \eqn{w_t'=(y_{t-1}',\dots,y_{t-p}')} is a \eqn{kp}-dimensional vector of lagged endogenous variables
#' and \eqn{q_t'=(d_{t-1}',\dots,d_{t-p}')} is a \eqn{qp}-dimensional vector of lagged deterministic (exogenous) variables,
#' \eqn{I_p} is the \eqn{(p \times p)} identity matrix and \eqn{\otimes} denotes the Kronecker product.
#' The goal is to estimate the parameters \eqn{\Theta = \begin{bmatrix} \beta & \Psi & \Sigma_u \end{bmatrix}}, and as such priors are needed.
#' Please see \link{priors} for more details.
#' 
#' For the innovations to the model, in the case of the homoscedastic (the original) steady-state BVAR, they are \eqn{u_t \overset{\text{iid}}{\sim} \mathrm{N_k}(0,\Sigma_u)}.
#' However, with stochastic volatility, there is instead a time-varying covariance matrix \eqn{u_t \sim \mathrm{N_k}(0,\Sigma_{u,t})}.
#' The innovations then take the form
#' 
#' \deqn{\begin{aligned} u_t &= A^{-1} \Lambda^{0.5}_t \epsilon_t \\
#' \epsilon_t &\overset{\text{iid}}{\sim} \mathrm{N}(0, \mathrm{I}_k)\end{aligned}}
#' 
#' where \eqn{A} is a lower triangular matrix with ones on the diagonal that describes 
#' the contemporaneous interaction of the endogenous variables, and
#' 
#' \deqn{\Lambda_t = \mathrm{diag}(\lambda_{1,t},\dots,\lambda_{k,t})}
#' 
#' contains the time-varying volatilities.
#' For the \code{AR1} stochastic volatility specification, the log volatilities follow AR(1) processes
#' 
#' \deqn{\ln \lambda_{i,t} = \gamma_{0,i} + \gamma_{1,i} \ln \lambda_{i,t-1} + \nu_{i,t}, \ i=1,\dots,k}
#' 
#' where the log volatility AR(1) processes are restricted to the stationary region, i.e. \eqn{|\gamma_{1,i}|<1 \ \forall i}.
#' For the \code{RW} stochastic volatility specification, the log volatilities follow (driftless) Random Walk processes
#' 
#' \deqn{\gamma_{0,i}=0, \ \gamma_{1,i}=1 \ \forall i}
#' 
#' The innovations to the log volatilities follow in the AR(1) case
#' 
#' \deqn{\nu_{t} = (\nu_{1,t},\dots,\nu_{k,t})'\overset{\text{iid}}{\sim} \mathrm{N}(0, \Phi)}
#' 
#' where \eqn{\Phi} \emph{is not diagonal} and as such the innovations to the log
#' volatilities are allowed to be correlated across variables. For the Random Walk case, \eqn{\Phi} \emph{is diagonal}
#' with variances \eqn{\phi_i} for \eqn{i=1,\dots,k}.
#' 
#' Note that under both stochastic volatility specifications, the time-varying covariance matrix is
#' 
#' \deqn{\Sigma_{u,t} = A^{-1} \Lambda_t (A^{-1})'}
#' 
#' For details on the homoscedastic steady-state BVAR model, see Villani (2009).
#' For details on the Random Walk stochastic volatility steady-state BVAR model, see Clark (2011).
#' See Carriero, Clark, and Marcellino (2024) for the above-mentioned AR(1) stochastic volatility
#' specification applied to a conventional BVAR.
#' 
#' To see examples for each type of steady-state BVAR model, simply run one of the following
#' 
#' \itemize{
#'   \item \code{vignette("Homoscedastic-steady-state-BVAR")}
#'   \item \code{vignette("RW-stochastic-volatility-steady-state-BVAR")}
#'   \item \code{vignette("AR1-stochastic-volatility-steady-state-BVAR")}
#' }
#' 
#'
#' @param data A numeric matrix or time series of data where each column is a
#'   variable and each row is a time period.
#'
#' @return A steady-state \code{bvar} object.
#' 
#' @references
#' Carriero, A., Clark, T. E., and Marcellino, M. (2024).
#' Capturing macro-economic tail risks with Bayesian vector autoregressions. 
#' \emph{Journal of Money, Credit and Banking}, 56(5), pp. 1099–1127.
#' 
#' Clark, T. E. (2011). Real-time density forecasts from Bayesian vector autoregressions
#' with stochastic volatility. \emph{Journal of Business & Economic Statistics}, 29(3), pp. 327–341.
#' 
#' Karlsson, S. (2013). Forecasting with Bayesian vector autoregression.
#' In: Elliott, G. and Timmermann, A. (eds), \emph{Handbook of Economic Forecasting}.
#' Elsevier B.V., Vol. 2, Part B, pp. 791–897.
#' 
#' Villani, M. (2009). Steady-state priors for vector autoregressions.
#' \emph{Journal of Applied Econometrics}, 24(4), pp. 630–650. 
#' 
#'
#'
#' @export
#'
#' @examples
#' yt <- matrix(rnorm(50), 25, 2)
#' 
#' bvar_obj <- bvar(data = yt)
bvar <- function(data) {
  
  if (!is.matrix(data) && !is.ts(data)) {
    stop("data must be a matrix or time series object")
  }
  
  obj <- list(
    data    = data,
    setup   = NULL,
    priors  = NULL,
    fit     = NULL,
    predict = list()
  )
  
  class(obj) <- "bvar"
  return(obj)
}