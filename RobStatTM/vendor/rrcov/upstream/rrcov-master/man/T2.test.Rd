\name{T2.test}
\alias{T2.test}
\alias{T2.test.default}
\alias{T2.test.formula}
\title{Robust Hotelling T2 test}
\description{
  Performs one and two sample Hotelling T2 tests as well as robust one-sample Hotelling T2 test
}
\usage{
T2.test(x, \dots)

\method{T2.test}{default}(x, y = NULL, mu = 0, conf.level = 0.95, method=c("c", "mcd"), \dots)

\method{T2.test}{formula}(formula, data, subset, na.action, \dots)
}
\arguments{
  \item{x}{a (non-empty) numeric data frame or matrix.}
  \item{y}{an optional (non-empty) numeric data frame or matrix.}
  \item{mu}{an optional (non-empty) numeric vector of data values (or
    a single number which will be repeated p times) indicating the 
    true value of the mean (or difference in means if you are 
    performing a two sample test).}
  \item{conf.level}{confidence level of the interval}
  \item{method}{the method to be used - 'c' for sample mean and covariance matrix and
    'mcd' for minimum covariance determinant estimator. A two-sample MCD based T2-test is not 
    yet implemented.}
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs}
    is a numeric data frame or matrix giving the observations and \code{rhs} a factor
    with two levels giving the corresponding groups.}
  \item{data}{an optional matrix or data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used (currently not used)}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")} (currently only "na.rm" used)}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  The formula interface is only applicable for the two-sample tests.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the T2-statistic.}
  \item{parameter}{the degrees of freedom for the T2-statistic.}
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the mean vector appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{the estimated mean vector or vectors depending on
    whether it was a one-sample test or a two-sample test.}
  \item{null.value}{the specified hypothesized value of the mean or mean
    difference depending on whether it was a one-sample test or a
    two-sample test.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string indicating what type of T2-test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{ 
 Willems G., Pison G., Rousseeuw P. and Van Aelst S.  (2002), 
 A robust hotelling test, 
 \emph{Metrika}, \bold{55}, 125--138.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}

%\note{}

\seealso{ \code{\link{CovMcd}}}

\examples{

## One-sample classical test
data(delivery)
delivery.x <- delivery[,1:2]
T2.test(delivery.x)

## One-sample robust test
data(delivery)
delivery.x <- delivery[,1:2]
T2.test(delivery.x, method="mcd")

## Two-sample classical test
data(hemophilia)
grp <-as.factor(hemophilia[,3])
x <- hemophilia[which(grp==levels(grp)[1]),1:2]
y <- hemophilia[which(grp==levels(grp)[2]),1:2]
T2.test(x,y)

## or using the formula interface
T2.test(as.matrix(hemophilia[,-3])~hemophilia[,3])


\dontrun{
## Two-sample robust test
T2.test(x,y, method="mcd")    ## error - not yet implemented
}}
\keyword{htest}
\keyword{multivariate}
