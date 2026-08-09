\name{jafwd}
\alias{jafwd}
\title{jafwd}
\usage{
jafwd(pars, resfn = NULL, bdmsk = NULL, resbest = NULL, ndstep = 1e-07, ...)
}
\arguments{
\item{pars}{a named numeric vector of parameters to the model}

\item{resfn}{a function to compute a vector of residuals}

\item{bdmsk}{Vector defining bounds and masks. Default is \code{NULL}}

\item{resbest}{If supplied, a vector of the residuals at the parameters
\code{pars} to save re-evaluation.}

\item{ndstep}{A tolerance used to alter parameters to compute numerical
approximations to derivatives. Default \code{1e-7}.}

\item{\dots}{Extra information needed to compute the residuals}
}
\description{
approximate Jacobian via forward differences
}
\author{
J C Nash 2014-7-16   nashjc _at_ uottawa.ca
}
