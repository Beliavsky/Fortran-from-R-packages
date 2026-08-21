\name{ProfLogLikelihood}
\Rdversion{2.1}
\alias{ProfLogLikelihood}
\alias{prof.gb2}
\alias{proflogl.gb2}
\alias{profscores.gb2}
\title{
Profile Log-likelihood of the GB2 Distribution
}
\description{
Expression of the parameters \code{shape2} \eqn{=p} and \code{shape3} \eqn{=q}  of the GB2 distribution as functions of \code{shape1} \eqn{=a} and \code{scale} \eqn{=b}, 
profile log-likelihood of the GB2 distribution, scores of the profile log-likelihood.
}
\usage{
prof.gb2(x, shape1, scale, w=rep(1, length(x)))
proflogl.gb2(x, shape1, scale, w=rep(1, length(x)))
profscores.gb2(x, shape1, scale, w=rep(1, length(x)))
}
\arguments{
 \item{x}{numeric; vector of data values.}
 \item{shape1}{numeric; positive parameter.}
 \item{scale}{numeric; positive parameter.}
 \item{w}{numeric; vector of weights. Must have the same length as \code{x}. By default \code{w} is a vector of 1.}

}
\details{
Using the full log-likelihood equations for the GB2 distribution, the parameters \eqn{p} and \eqn{q} can be estimated as functions of \eqn{a} and \eqn{b}.  These functions are plugged into the log-likelihood expression,
which becomes a function of \eqn{a} and \eqn{b} only. This is obtained by reparametrizing the GB2, i.e. we set \eqn{r=\frac{p}{p+q}} and \eqn{s=p+q}.  More details can be found in Graf (2009).
}
\value{
\code{prof} returns a vector containing the values of \eqn{r}, \eqn{s}, \eqn{p}, \eqn{q} as well as two other parameters used in the calculation of the profile log-likelihood and its first derivatives.
\code{proflogl.gb2} returns the value of the profile log-likelihood and \code{profscores.gb2} returns the vector of the first derivatives of the profile log-likelihhod with respect to \eqn{a} and \eqn{b}.
}
\references{
Graf, M. (2009)
 The Log-Likelihood of the Generalized Beta Distribution of the Second Kind.
 \emph{working paper}, SFSO.
}
\author{
Monique Graf and Desislava Nedyalkova
}
\keyword{distribution}

