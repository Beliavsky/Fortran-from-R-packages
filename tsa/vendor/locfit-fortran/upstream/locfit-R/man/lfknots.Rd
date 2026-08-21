\name{lfknots}
\alias{lfknots}
\title{
Extraction of fit-point information from a Locfit object.
}
\usage{
lfknots(x, tr, what = c("x", "coef", "h", "nlx"), delete.pv = TRUE)
}
\description{
  Extracts information, such as fitted values, influence functions
  from a \code{"locfit"} object.
}
\arguments{
\item{x}{
Fitted object from \code{\link{locfit}()}.
}
\item{tr}{
Back transformation. Default is the invers link function from the Locfit
object.
}
\item{what}{
What to return; default is \code{c("x","coef","h","nlx")}.
Allowed fields are \code{x} (fit points); \code{coef} (fitted values);
\code{f1} (local slope); \code{nlx} (length of the weight diagram);
\code{nlx1} (estimated derivative of \code{nlx}); \code{se} (standard errors);
\code{infl} (influence function); \code{infla} (slope of influence function);
\code{lik} (maximixed local log-likelihood and local degrees of freedom);
\code{h} (bandwidth) and \code{deg} (degree of fit).
}
\item{delete.pv}{If \code{T}, pseudo-vertices are deleted.}
}
\value{
A matrix with one row for each fit point. Columns correspond to
the specified \code{what} vector; some fields contribute multiple columns.
}
\keyword{smooth}
