\name{findSubexprs}
\alias{findSubexprs}
\title{findSubexprs}
\usage{
findSubexprs(expr, simplify = FALSE, tag = ".expr", verbose = FALSE, ...)
}
\arguments{
\item{expr}{An expression represented in a variety of ways. See Details.}

\item{simplify}{Whether to simplify.}

\item{tag}{to be attached to the returned object(s)}

\item{verbose}{If \code{TRUE}, then diagnostic output will be printed as derivatives
and simplifications are recognized.}

\item{...}{Additional arguments, passed to \code{\link{nlsSimplify}}}.
}
\description{
Find common sub-expressions in \code{expr}.
}
\details{
If \code{expr} contains repeated subexpressions, this function
returns an expression to evaluate them just once and store them
in temporary variables with names starting with the \code{tag}.

This function is used internally by \code{\link{codeDeriv}}.
}
\examples{
findSubexprs(expression(x^2 + x^2), tag = ".example")
}
\keyword{internal}
