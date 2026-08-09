\name{resgr}
\alias{resgr}
\title{resgr}
\description{
Computes the gradient of the sum of squares function for nonlinear least
squares where \code{resfn} and \code{jacfn} supply the residuals and Jacobian
}
\usage{
resgr(prm, resfn, jacfn, ...)
}
\arguments{
\item{prm}{a numeric vector of parameters to the model}

\item{resfn}{a function to compute a vector of residuals}

\item{jacfn}{a function to compute the Jacobian of the sum of squares. If
    the value is quoted, then the function is assumed to be a numerical
    approximation. Currently one of \code{"jafwd"}, \code{"jaback"}, 
    \code{"jacentral"}, or  \code{"jand"}.}

\item{...}{Extra information needed to compute the residuals}
}
\details{
Does NOT (yet) handle calling of code built into selfStart models. That
is, codes that in \code{nlxb} employ control \code{japprox="SSJac"}.
}
\value{
The main object returned is the numeric vector of residuals computed at \code{prm}
by means of the function \code{resfn}. 
There are \code{Jacobian} and \code{gradient} attributes giving the Jacobian
(matrix of 1st partial derivatives whose row i contains the partial derivative
of the i'th residual w.r.t. each of the parameters) and the gradient of the
sum of squared residuals w.r.t. each of the parameters. Moreover, the Jacobian
is repeated within the \code{gradient} attribute of the Jacobian. This somewhat
bizarre structure is present for compatibility with \code{nls()} and some other
legacy functions, as well as to simplify the call to \code{nlfb()}.
}
\author{
J C Nash 2014-7-16 revised 2022-11-22   nashjc _at_ uottawa.ca
}
