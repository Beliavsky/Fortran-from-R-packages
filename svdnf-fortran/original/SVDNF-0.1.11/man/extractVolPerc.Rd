\name{extractVolPerc.SVDNF}
\alias{extractVolPerc.SVDNF}
\alias{extractVolPerc}
\alias{extractVolPerc.DNFOptim}
\title{
Extract Filtering and Prediction Distribution Percentiles}
\description{
Function to extract filtering and prediction distribution percentiles from \code{SVDNF} and \code{DNFOptim} objects.
}
\usage{
\method{extractVolPerc}{SVDNF}(x, p = 0.5, pred = F, ...)
}
\arguments{
  \item{x}{An \code{SVDNF} or \code{DNFOptim} object.}
  \item{p}{Distribution percentile to return.}
  \item{pred}{Return the prediction distribution percentile if \code{pred = TRUE},
  otherwise, return the filtering distribution percentile.}
  \item{...}{Other parameters to be passed through to function.}
}
\value{
Return a vector of the \code{p}th percentile filtering or prediction distribution volatility factor values depending on whether \code{pred = TRUE} or not. 
}
\examples{

set.seed(1)

# Define the built-in model using the dynamicsSVM function
Heston_mod <- dynamicsSVM(model = "Heston")

# Generate the data from the Duffie, Pan, and Singleton model
Heston_sim <- modelSim(t = 10, dynamics = Heston_mod)
Heston_dnf <- DNF(dynamics = Heston_mod, data = Heston_sim$returns)

extractVolPerc(Heston_dnf, p = 0.75)
}

