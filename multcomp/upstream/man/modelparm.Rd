\name{modelparm}
\alias{modelparm}
\title{ Generic Accessor Function for Model Parameters }
\description{
  Extract model parameters and their covariance matrix as well as 
  degrees of freedom (if available) from a fitted model.
}
\usage{
modelparm(model, coef., vcov., df, ...)
}
\arguments{
  \item{model}{a fitted model,
               for example an object returned by \code{\link{lm}},
               \code{\link{glm}}, \code{\link{aov}}, \code{\link[survival]{survreg}},
               \code{\link[fixest]{fixest}}, or \code{\link[lme4]{lmer}} etc. }
  \item{coef.}{ an accessor function for the model parameters. Alternatively, the vector of coefficients.}
  \item{vcov.}{ an accessor function for the covariance matrix of the
                model parameters. Alternatively, the covariance matrix directly.}
  \item{df}{ an optional specification of the degrees of freedom to be used in
             subsequent computations. }
  \item{\dots}{ additional arguments, currently ignored. }
}
\details{
  
  One can't expect \code{\link{coef}} and \code{\link{vcov}} methods 
  for arbitrary models to
  return a vector of \eqn{p} fixed effects model parameters (\code{coef})
  and corresponding \eqn{p \times p} covariance matrix (\code{vcov}).

  The \code{coef.} and \code{vcov.} arguments can be used to define
  modified \code{coef} or \code{vcov} methods for a specific model.
  Methods for \code{\link[lme4]{lmer}}, \code{\link[fixest]{fixest}}, 
  and \code{\link[survival]{survreg}} objects are available (internally).

  For objects inheriting from class \code{\link{lm}} the degrees of 
  freedom are determined from \code{model} and the corresponding
  multivariate t distribution is used by all methods to \code{\link{glht}}
  objects. By default, the asymptotic multivariate normal distribution
  is used in all other cases unless \code{df} is specified by the user.

}
\value{
  An object of class \code{modelparm} with elements
  \item{coef}{model parameters}
  \item{vcov}{covariance matrix of model parameters}
  \item{df}{degrees of freedom}
}
\keyword{misc}
