\name{AutocorTest}
\alias{AutocorTest}
\title{ Box-Ljung autocorrelation test }
\description{
  Ljung-Box test for autocorrelation
}
\usage{
AutocorTest(x, lag = ceiling(log(length(x))),
            type = c("Ljung-Box", "Box-Pierce", "rank"),
            df = lag )
}
\arguments{
  \item{x}{
    a numeric vector or a univariate time series
  }
  \item{lag}{
    the statistic will be based on 'lag' autocorrelation coefficients.
    Tsay (p. 27-28) says, 'Simulation studies suggest that the choice of
    [lag = log(length(x))] provides better power performance.  This
    general rule needs modification in analysis of seasonal time series
    for which autocorrelations with lags at multiples of the seasonality
    are more important.'
  }
  \item{type}{
    which Box.test 'type' should be used?  Partial matching is used.

    The 'rank' alternative computes 'Ljung-Box' on rank(x);  see Burns
    (2002) and references therein.

    NOTE:  The default 'Ljung-Box' type generally seems to be more
    accurate and popular than the earlier 'Box-Pierce', which is however
    the default for 'Box.test'.
  }
  \item{df}{
    a positive number giving the degrees of freedom for the reference
    chi-squre distribution used to compute the p-value for the
    statistic.

    This makes it easy to call AutocorTest with the residuals from a fit
    and have the p-value computed with reference to a chi-square with
    degrees of freedom different from "lag".  See the discussion
    of degrees of freedom for 'Box.test in   \code{\link{ARIMA}}.
  }
}
\details{
  This is provided for compatibility with 'autocorTest' in the S-Plus
  script in Tsay (p. 30).  It is a wrapper for the R function
  \code{\link{Box.test}}.
}
\value{
  a list of class 'htest' containing the following components:

  \item{statistic}{
    a number giving the value of the test statistic.
  }
  \item{paramter}{
    a number giving the degrees of freedom of the approximate
    chi-squared distribution of the test statistic used to compute the
    p.value.
  }
  \item{p.value}{
    the p-value of the test.
  }
  \item{method}{
    a character string indicating which type of test was performed.
    If(df != lag), this character string ends with
    paste("(lag = ", lag, ")", sep="").
  }
}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley)

  Patrick Burns (2002) 'Robustness of the Ljung-Box Test and its Rank
  Equivalent',
  \url{https://www.burns-stat.com/pages/Working/ljungbox.pdf}, accessed
  2007.12.29.
}
\seealso{
  \code{\link{Box.test}}
  \code{\link{ARIMA}}
}
\examples{
data(m.ibm2697)
AutocorTest(m.ibm2697, 5)

AT4 <- AutocorTest(m.ibm2697, 5, df=4)
str(AT4) # $method = "Box-Ljung test (lag = 5)"
}
\keyword{ts}
