\name{pars.dynamicsSVM}
\alias{pars}
\alias{pars.dynamicsSVM}
\title{Parameters Names and Order for Stochastic Volatility Models with Jumps}
\description{
The \code{pars} function returns the names of the parameters in a \code{dynamicsSVM} object in the order that they should be passed to the \code{DNFOptim} function.}
\usage{
\method{pars}{dynamicsSVM}(dynamics, rho = NULL,
delta = NULL, alpha = NULL, rho_z = NULL, nu = NULL,
jump_params_list = "dummy", ...)
}
\arguments{
  \item{dynamics}{A \code{dynamicsSVM} object representing the model dynamics to be used for which we want to get parameter names.}
  \item{rho, delta, alpha, rho_z, nu}{See \code{help(dynamicsSVM)} for a description of each of these arguments individually. These arguments should be used only for custom models and can be fixed to a certain
value (e.g., \code{rho = -0.75}). If they are to be estimated, they should be set to \code{'var'} (e.g., to estimate \code{rho} set \code{rho = 'var'}) and include it in the vector par to be passed to the \code{DNFOptim} function.
See Note for more details on the order in which custom models should receive parameters.}
\item{jump_params_list}{List of the names of the arguments in the jump parameter distribution in the order that they are used by the \code{jump_dist} function.}
  \item{...}{Other parameters to be passed through to function.}

}
\value{Returns a vector of strings with the names of the parameters in the given \code{dynamicsSVM} object. The parameters names are returned in the order the parameters should be passed to the \code{DNFOptim} function 
}
\examples{
mod <- dynamicsSVM(model = "Taylor")
pars(mod)
}
