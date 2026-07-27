\name{dynamicsSVM}
\alias{dynamicsSVM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic Volatility Models Dynamics 
}
\description{
\code{dynamicsSVM} creates stochastic volatility model dynamics by either choosing from a set of built-in model dynamics or using custom drift and diffusion functions, as well as custom jump distributions. See Note for information about how to define custom functions.
}
\usage{
dynamicsSVM(mu = 0.038, kappa = 3.689, theta = 0.032, sigma = 0.446,
rho = -0.745, omega = 5.125, delta = 0.03, alpha = -0.014,
rho_z = -1.809, nu = 0.004, p = 0.01, phi = 0.965, h = 1/252, coefs = NULL,
model = "Heston", mu_x, mu_y, sigma_x, sigma_y,
jump_dist = rpois, jump_density = dpois, jump_params = 0,
mu_x_params, mu_y_params, sigma_x_params, sigma_y_params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{mu}{Annual expected rate of return.}
 \item{kappa}{Variance rate of mean reversion.}
 \item{theta}{Unconditional mean variance.}
 \item{sigma}{Volatility of the variance.}
 \item{rho}{Correlation between the return and the variance noise terms.}
 \item{omega}{Jump arrival intensity for models with Poisson jumps.}
 \item{delta}{Standard deviation of return jumps.}
 \item{alpha}{Mean of return jumps.}
 \item{rho_z}{Pseudo-correlation parameter between return and variance jumps.}
 \item{nu}{Mean for variance jumps.}
 \item{p}{Jump probability for models with Bernoulli jumps.}
 \item{phi}{Volatility persistence parameter.}
 \item{h}{Time interval between observations (e.g., \code{h = 1/252} for daily data).}
 \item{coefs}{Vector of regression coefficients for factor stochastic volatility models. This vector should be of the same lenght as the number of explanatory variables provided.}
 \item{model}{Model used by the discrete nonlinear filter. The options are \code{"Heston"},    \code{"Bates"}, \code{"DuffiePanSingleton"}, \code{"Taylor"}, \code{"TaylorWithLeverage"}, \cr  \code{"PittMalikDoucet"}, \code{"CAPM_SV"} and \code{"Custom"}. If \code{model = "Custom"}, users should pass the drift  functions (i.e., \code{mu_x} and \code{mu_y}), the diffusion functions (i.e., \code{sigma_x} and  \code{sigma_y}), and the jump distribution, (i.e., \code{jump_dist}) as well as their parameters to the \code{dynamicsSVM} function. See Examples.}
 \item{mu_x}{Function for variance drift (to be used with a custom model).}
 \item{mu_y}{Function for returns drift (to be used with a custom model).}
 \item{sigma_x}{Function for variance diffusion (to be used with a custom model).}
 \item{sigma_y}{Function for returns diffusion (to be used with a custom model).}
 \item{jump_dist}{Distribution used to generate return or volatility jumps at each  timestep (if both types of jumps are in the model, they are assumed to occur simulaneously).}
  \item{jump_density}{Probability mass function used to compute the probability of return or volatility jumps at each  timestep (if both types of jumps are in the model, they are assumed to occur simulaneously).}
 \item{jump_params}{List of parameters to be used as arguments in the \code{jump_dist} and \code{jump_density} function (parameters should be listed in the order that jump_dist uses them).}
 \item{mu_x_params}{List of parameters to be used as arguments in the \code{mu_x} function (parameters should be listed in the order that \code{mu_x} uses them).}
 \item{mu_y_params}{List of parameters to be used as arguments in the \code{mu_y} function (parameters should be listed in the order that \code{mu_y} uses them).}
 \item{sigma_x_params}{List of parameters to be used as arguments in the \code{sigma_x} function (parameters should be listed in the order that \code{sigma_x} uses them).}
 \item{sigma_y_params}{List of parameters to be used as arguments in the \code{sigma_y} function (parameters should be listed in the order that \code{sigma_y} uses them).}
}

\note{
Custom functions should have \code{x} (the volatility factor) as well as the function's other parameters as arguments.

If the custom function does not use any parameters, one should include an argument called dummy and its parameters as a \code{list(0)}. For example, for the Taylor model, 

\code{sigma_y_taylor <- function(x, dummy) {
  return(exp(x / 2))
}} \cr
\code{sigma_y_params <- list(0)}

It should also be noted that the custom function is a vector for \code{x}. This means that users should use vectorized version of functions. For example, \code{pmax(0,x)} instead of \code{max(0,x)} as code seen in the Example section below. 
}


\value{
Returns an object of type dynamicsSVM.
}
\examples{
# Create a dynamicsSVM object with model DuffiePanSingleton and default parameters
DuffiePanSingleton_mod <- dynamicsSVM(model = "DuffiePanSingleton") 

# Here, we define the same DuffiePanSingleton model 
# using the custom model option.

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
custom_DPS <- dynamicsSVM(model = 'Custom',
  mu_x = mu_x, mu_y = mu_y, sigma_x = sigma_x, sigma_y = sigma_y,
  mu_x_params = mu_x_params, mu_y_params = mu_y_params,
  sigma_x_params = sigma_x_params, sigma_y_params = sigma_y_params,
  jump_params = jump_params, jump_dist = jump_dist, jump_density = jump_density,
  nu = nu, rho_z = rho_z)
}