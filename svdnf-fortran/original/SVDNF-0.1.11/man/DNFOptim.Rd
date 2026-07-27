\name{DNFOptim.dynamicsSVM}
\alias{DNFOptim.dynamicsSVM}
\alias{DNFOptim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete Nonlinear Filter Maximum Likelihood Estimation Function
}
\description{
The \code{DNFOptim} function finds maximum likelihood estimates for stochastic volatility models parameters using the \code{DNF} function.
}
\usage{
\method{DNFOptim}{dynamicsSVM}(dynamics, data, par, factors, tol, N = 50, K = 20, R = 1,
  grids = 'Default', 
  rho = 0, delta = 0, alpha = 0, rho_z = 0, nu = 0, jump_params_list = "dummy",
  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dynamics}{A dynamicsSVM object representing the model dynamics to be used by the optimizer to find maximum likelihood parameter estimates.This should be a vector or an \code{xts} object.}
  \item{data}{A series of asset returns for which we want to find maximum likelihood estimates.}
 \item{par}{Initial values for the parameters to be optimized over. Information about how to pass the initial parameters is given in the 'Note' section.}
  \item{factors}{Series of values taken by \code{d} explanatory variables. This should be a matrix or an \code{xts} object with \code{d} rows and \code{T} columns.}
  \item{tol}{Tolerance hyperparameter for the optimization. The \code{optim} function will be rerun until the difference between the maximum values for the likelihood function is less than \code{tol}. After the first \code{optim} call, which uses the initial parameters, the next optimizers are run using the previous optimizer's maximum likelihood parameter estimates as initial values.}
  \item{N}{Number of nodes in the variance grid.}
  \item{K}{Number of nodes in the jump size grid.}
  \item{grids}{Grids to be used for numerical integration by the \code{DNF} function. The \code{DNF} function creates grids for built-in models. However, this arguments must be provided for custom models. It should contain a list of three sequences: \code{var_mid_points} (variance mid-point sequence), \code{j_nums} (sequence for the number of jumps), and \code{jump_mid_points} (jump mid-point sequence). If there are no variance jumps in the model, set \code{jump_mid_points} equal to zero. If there are no jumps in the model, both \code{j_nums} and \code{jump_mid_points} should be set to zero.}
  \item{R}{Maximum number of jumps used in the numerical integration at each timestep.} 
  \item{rho, delta, alpha, rho_z, nu}{See \code{help(dynamicsSVM)} for a description of each of these arguments individually. These arguments should be used only for custom models and can be fixed to a certain
value (e.g., \code{rho = -0.75}). If they are estimated, they should be set to \code{'var'} (e.g., to estimate \code{rho} set \code{rho = 'var'}) and include it in the vector par to be passed to the \code{optim} function.
See Note for more details on the order in which custom models should receive parameters.}
\item{jump_params_list}{List of the names of the arguments in the jump parameter distribution in the order that they are used by the \code{jump_dist} function. This is used by \code{DNFOptim} to check for parameters that occur both in the \code{jump_dist} function and as arguments in drift or diffusion functions.}
  \item{...}{Further arguments to be passed to the \code{optim} function. See Note.}
}
\note{
When passing the initial parameter vector \code{par} to the \code{optim} function (via \code{...}), the parameters should follow a specific order.

For the \code{PittMalikDoucet} model, the parameters should be in the following order: \code{phi}, \code{theta}, \code{sigma}, \code{rho}, \code{delta}, \code{alpha}, and \code{p}.

For the \code{DuffiePanSingleton} model, the parameters should be in the following order: \code{mu}, \code{alpha}, \code{delta}, \code{rho_z}, \code{nu}, \code{omega}, \code{kappa}, \code{theta}, \code{sigma}, and \code{rho}.

For the \code{CAPM_SV} model, the parameters should be in the following order: \code{c_0}, \code{c_1}, \code{phi}, \code{theta}, and \code{sigma}.

All other built-in models can be seen as being nested within these models (i.e., \code{Heston} and \code{Bates} models are nested in the \code{DuffiePanSingleton} model, while \code{Taylor} and \cr \code{TaylorWithLeverage} are nested in the \code{PittMalikDoucet} model). Their parameters should be passed in the same order as those in the more general models, minus the parameters not found in these nested models.

For example, the \code{Taylor} model contains neither jumps nor correlation between volatility and returns innovations. Thus, its three parameters are passed in the order: \code{phi}, \code{theta}, and \code{sigma}.

When \code{models = "Custom"}, parameters should be passed in the following order: \code{mu_y_params}, \code{sigma_y_params}, \code{mu_x_params}, \code{sigma_x_params}, \code{rho}, \code{delta}, \code{alpha}, \code{rho_z}, \code{nu}, and \code{jump_params}. If an argument is repeated (e.g., both \code{mu_y_params} and  \code{sigma_y_params} use the same parameter), write it only when it first appears in the custom model order.

}
\value{
\item{optim}{Returns a list obtained from \proglang{R}'s \code{optim} function. See \code{help(optim)} for details about the output.}
\item{SVDNF}{Returns a SVDNF object obtained from running the \code{DNF} function at the MLE parameter values. See \code{help(DNF)} for details about the output}
\item{rho, delta, alpha, rho_z, nu}{See \code{help(dynamicsSVM)} for a description of each of these arguments individually. If they are estimated, they are set to \code{'var'}. If the parameters were fixed during the estimation, this will return the value at which they were fixed.}
}
\references{
 \proglang{R} Core Team (2019). \proglang{R}: A language and environment for statistical computing. \proglang{R} Foundation
  for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
  }
\examples{
set.seed(1)

# Generating return data
Taylor_mod <- dynamicsSVM(model = "Taylor", phi = 0.9,
  theta = -7.36, sigma = 0.363)
Taylor_sim <- modelSim(t = 30, dynamics = Taylor_mod, init_vol = -7.36)
plot(Taylor_sim$volatility_factor, type = 'l')
plot(Taylor_sim$returns, type = 'l')

# Initial values and optimization bounds
init_par <- c(0.7, -5, 0.3)
lower <- c(0.01, -20, 0.1); upper <- c(0.99, 0, 1)

# Running DNFOptim to get MLEs
optim_test <- DNFOptim(data = Taylor_sim$returns,
  dynamics = Taylor_mod,
  par = init_par, lower = lower, upper = upper, method = "L-BFGS-B")
                
# Parameter estimates
summary(optim_test)

# Plot prediction and filtering distributions
plot(optim_test)

}
