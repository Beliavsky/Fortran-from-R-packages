\name{nlsrSS}
\alias{nlsrSS}
\title{nlsrSS - solve selfStart nonlinear least squares with nlsr package}
\usage{
nlsrSS(formula, data)
}
\arguments{
\item{formula}{a model formula incorporating a selfStart function in the right hand side}

\item{data}{a data frame with named columns that allow evaluation of the \code{formula}}
}
\value{
A solution object of class \code{nlsr}.

List of solution elements
  
\item{resid}{weighted residuals at the proposed solution}
\item{jacobian}{Jacobian matrix at the proposed solution}
\item{feval}{residual function evaluations used to reach solution from starting parameters}
\item{jeval}{Jacobian function (or approximation) evaluations used}
\item{coefficients}{a named vector of proposed solution parameters}
\item{ssquares}{weighted sum of squared residuals (often the deviance)}
\item{lower}{lower bounds on parameters}
\item{upper}{upper bounds on parameters}
\item{maskidx}{vector if indices of fixed (masked) parameters}
\item{weights}{specified weights on observations}
\item{formula}{the modeling formula}
\item{resfn}{the residual function (unweighted) based on the formula}

}
\description{
This function uses the \code{getInitial()} function to estimate starting parameters for
a Gauss-Newton iteration, then calls \code{nlsr::nlxb()} appropriately to find a solution
to the required nonlinear least squares problem. 
}
\author{J C Nash 2022-9-14   nashjc _at_ uottawa.ca}
