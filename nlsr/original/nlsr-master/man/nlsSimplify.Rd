\name{nlsSimplify}
\alias{nlsSimplify}
\title{nlsSimplify}
\usage{
nlsSimplify(expr, simpEnv = sysSimplifications, verbose = FALSE)
}
\arguments{
\item{expr}{An expression represented in a variety of ways. See Details.}

\item{simpEnv}{The environment in which simplifications are stored.}

\item{verbose}{If \code{TRUE}, then diagnostic output will be printed as simplifications are recognized.}
}
\description{
Try to simplify an expression.
}
\details{
The \code{expr} can be an expression vector or other language object.
If it
is a complex expression (e.g. \code{(a + b) + c}), then simplifications
will be applied recursively.  Simplifications are applied 
from the database created by \code{\link{newSimplification}},
and the new expression (or call...) is returned.

This function is mainly for internal use by \code{\link{nlsDeriv}}.
}
\examples{
nlsSimplify(quote(x + 0 + y*1), verbose = TRUE)
}
\keyword{internal}