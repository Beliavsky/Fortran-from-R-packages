\name{predict.nlsr}
\alias{predict.nlsr}
\title{predict.nlsr}
\usage{
\method{predict}{nlsr}(object = NULL, newdata = list(), ...)
}
\arguments{
\item{object}{an object of class \code{nlsr}}

\item{newdata}{a dataframe containing columns that match the original dataframe
used to estimate the nonlinear model in the \code{nlsr} object}

\item{...}{additional data needed to evaluate the modeling functions
Default FALSE}
}
\description{
prepare and display predictions from an \code{nlsr} model
}
\author{
J C Nash 2014-7-16   nashjc _at_ uottawa.ca
}
