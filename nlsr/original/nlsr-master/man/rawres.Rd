\name{rawres}
\alias{rawres}
\title{rawres}
\description{
Prepare and display raw residuals of \code{nlsr} computations
NOTE: we use model - data form i.e., rhs - lhs
}
\usage{
rawres(object = NULL, data = parent.frame(), ...)
}
\arguments{
\item{object}{an object of class \code{nlsr}}

\item{data}{a data frame with the data for which fits are wanted}

\item{...}{additional data needed to evaluate the modeling functions}
}
\value{
A vector of the raw residuals
}
\author{
J C Nash 2014-7-16 revised 2022-11-22   nashjc _at_ uottawa.ca
}
