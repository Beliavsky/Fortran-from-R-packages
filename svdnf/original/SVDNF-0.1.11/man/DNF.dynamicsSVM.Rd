\name{DNF.dynamicsSVM}
\alias{DNF.dynamicsSVM}
\alias{DNF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete Nonlinear Filtering Algorithm for Stochastic Volatility Models
}
\description{
The \code{DNF} function applies the discrete nonlinear filter (DNF) of Kitagawa (1987) as per the implementation of Bégin & Boudreault (2020) to obtain likelihood evaluations and filtering distribution estimates for a wide class of stochastic volatility models.
}
\usage{
\method{DNF}{dynamicsSVM}(dynamics, data, factors = NULL, N = 50, K = 20, R = 1, grids, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dynamics}{A dynamicsSVM object representing the model dynamics to be used by the DNF.}
\item{data}{A series of asset returns for which we want to run the DNF. This should be a vector or an \code{xts} object.}
\item{factors}{Series of values taken by \code{d} explanatory variables. This should be a matrix or an \code{xts} object with \code{d} rows and \code{T} columns.}
\item{N}{Number of nodes in the variance grid.}
\item{K}{Number of nodes in the jump size grid.}
\item{R}{Maximum number of jumps used in the numerical integration at each timestep.} 
\item{grids}{Grids to be used for numerical integration by the \code{DNF} function. The \code{DNF} function creates grids for built-in models. However, this arguments must be provided for custom models. It should contain a list of three sequences: \code{var_mid_points} (variance mid-point sequence), \code{j_nums} (sequence for the number of jumps), and \code{jump_mid_points} (jump mid-point sequence). If there are no variance jumps in the model, set \code{jump_mid_points} equal to zero. If there are no jumps in the model, both \code{j_nums} and \code{jump_mid_points} should be set to zero.}
  \item{...}{Further arguments passed to or from other methods.}

}
\value{
\item{log_likelihood}{Log-likelihood evaluation based on the DNF.}
\item{filter_grid}{Grid of dimensions \code{N} by \code{T+1} that stores each time-step's filtering distributions (we assume the filtering distribution is uniform at \eqn{t = 0}).}
\item{likelihoods}{Likelihood contribution at each time-step throughout the series.}
\item{grids}{List of grids used for numerical integration by the DNF.}
\item{dynamics}{The model dynamics used by the DNF.}
\item{data}{The series of asset returns to which the DNF was applied.}

}
\references{
Bégin, J.F., Boudreault, M. (2021) Likelihood evaluation of jump-diffusion models using deterministic nonlinear filters. \emph{Journal of Computational and Graphical Statistics}, 30(2), 452--466.

Kitagawa, G. (1987) Non-Gaussian state-space modeling of nonstationary time series.
\emph{Journal of the American Statistical Association}, 82(400), 1032--1041.
}

\examples{
set.seed(1)
# Generate 200 returns from the DuffiePanSingleton model
DuffiePanSingleton_mod <- dynamicsSVM(model = "DuffiePanSingleton") 
DuffiePanSingleton_sim <- modelSim(t = 200, dynamics = DuffiePanSingleton_mod) 

# Run DNF on the data
dnf_filter <- DNF(data = DuffiePanSingleton_sim$returns,
  dynamics = DuffiePanSingleton_mod) 

# Print log-likelihood evaluation.
logLik(dnf_filter)

# Using a custom model.
# Here, we define the DuffiePanSingleton model as a custom model
# to get the same log-likelihood found using the built-in option

# Daily observations
h <- 1/252

# Parameter values 
mu <- 0.038; kappa <- 3.689; theta <- 0.032
sigma <- 0.446; rho <- -0.745; omega <- 5.125
delta <- 0.03; alpha <- -0.014; rho_z <- -1.809; nu <- 0.004

# Jump compensator
alpha_bar <- exp(alpha + 0.5 * delta^2) / (1 - rho_z * nu) - 1

# Returns drift and diffusion
mu_y <- function(x, mu, alpha_bar, omega, h) {
  return(h * (mu - x / 2 - alpha_bar * omega))
}
mu_y_params <- list(mu, alpha_bar, omega, h)
sigma_y <- function(x, h) {
  return(sqrt(h * pmax(x, 0)))
}
sigma_y_params <- list(h)

# Volatility factor drift and diffusion
mu_x <- function(x, kappa, theta, h) {
  return(x + h * kappa * (theta - pmax(0, x)))
}
mu_x_params <- list(kappa, theta, h)

sigma_x <- function(x, sigma, h) {
  return(sigma * sqrt(h * pmax(x, 0)))
}
sigma_x_params <- list(sigma, h)

# Jump distribution for the DuffiePanSingleton Model
jump_density <- dpois
jump_dist <- rpois
jump_params <- c(h * omega)

# Create the custom model
custom_mod <- dynamicsSVM(model = 'Custom',
  mu_x = mu_x, mu_y = mu_y, sigma_x = sigma_x, sigma_y = sigma_y,
  mu_x_params = mu_x_params, mu_y_params = mu_y_params,
  sigma_x_params = sigma_x_params, sigma_y_params = sigma_y_params,
  jump_params = jump_params, jump_dist = jump_dist, jump_density = jump_density,
  nu = nu, rho_z = rho_z, rho = rho)
# Define the grid for DNF
N <- 50; R <- 1; K <- 20
var_mid_points <- seq(from = sqrt(0.0000001),
  to = sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), length = N)^2
  
j_nums <- seq(from = 0, to = R, by = 1)

jump_mid_points <- seq(from = 0.000001, to = (3 + log(K)) * sqrt(R) * nu, length = K)

grids <- list(var_mid_points = var_mid_points,
  j_nums = j_nums, jump_mid_points = jump_mid_points)

# Run the DNF function with the custom model
dnf_custom <- DNF(data = DuffiePanSingleton_sim$returns, grids = grids, 
  dynamics = custom_mod)
  
# Check if we get the same log-likelihoods
dnf_custom$log_likelihood; dnf_filter$log_likelihood 
}
