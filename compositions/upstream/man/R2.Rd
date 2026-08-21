\name{R2}
\alias{R2}
\alias{R2.lm}
\alias{R2.default}
\alias{Rsquare}
\title{R square}
\description{
  The R2 measure of determination for linear models
}
\usage{
R2(object,...)
\S3method{R2}{lm}(object,...,adjust=TRUE,ref=0)
\S3method{R2}{default}(object,...,ref=0)
          }
\arguments{
  \item{object}{ a statistical model}
  \item{\dots}{further not yet used parameters}
  \item{adjust}{Logical, whether the estimate of R2 should be adjusted
    for the degrees of freedom of the model.}
  \item{ref}{A reference model for computation of a relative \eqn{R^2}.}
}
\value{
  The R2 measure of determination. 
}

\details{
  The \eqn{R^2}{R^2} measure of determination is defined as:
  \deqn{R^2=1-\frac{var(residuals)}{var(data)}}{1-var(residuals)/var(data)}

  and provides the portion of variance explained by the model. It is a
  number between 0 and 1, where 1 means the model perfectly explains the
  data and 0 means that the model has no better explanation of the data
  than a constant mean. In case of multivariate models metric variances
  are used. 
  
  If a reference model is given by \code{ref}, the variance of the
  residuals of that models rather than the variance of the data is
  used. The value of such a relative \eqn{R^2}{R^2} estimates how much
  of the residual variance is explained.
  
  If \code{adjust=TRUE} the unbiased estiamators for the variances are
  used, to avoid the automatisme that a more parameters automatically
  lead to a higher \eqn{R^2}{R^2}. 
  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{lm}}, \code{\link{mvar}}, \code{\link{AIC}} 
}
\examples{
data(Orange)
R2(lm(circumference~age,data=Orange))
R2(lm(log(circumference)~log(age),data=Orange))
}
\keyword{multivariate}
