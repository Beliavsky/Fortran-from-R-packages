\name{testRmACD}
\alias{testRmACD}

\title{
LM test of no Remaining ACD (Meitz and Terasvirta, 2006)
}
\description{
Tests if there is any remaining ACD structure in the residuals (currently only works for \code{model = "ACD"}).
}
\usage{
testRmACD(fitModel, pStar = 2, robust = TRUE)
}
\arguments{
  \item{fitModel}{
a fitted ACD model, i.e. an object of class "acdFit".
}
  \item{pStar}{
the number of alpha parameters in the alternative hypothesis. See \eqn{p*} under 'Details'.
}
  \item{robust}{
if \code{TRUE} the LM statistic will be calculated using the "robust" version, making its asymptotic behavior unaffected by possible misspecification of the error term distribution (Meitz and Terasvirta, 2006).
}
}
\details{

For the model
\deqn{x_i = \mu_i \phi_i \epsilon_i,}
\deqn{\mu_i = \omega + \sum_{j=1}^{p} \alpha_j x_{i-j} + \sum_{j=1}^{q} \beta_j \mu_{i-j},}
\deqn{\phi_i = 1 + \sum_{j=1}^{p*} \frac{x_{i-j}}{\mu_{i-j}},}

the function tests the null hypothesis

\deqn{H_0: \phi_i = 1.}

}
\value{
a list of:
\item{chi2}{
the value of the LM statistic.
}
\item{pv}{
the pvalue of the test statistic.
}
}
\references{
Meitz, M. and Terasvirta, T. (2006).
\emph{Evaluating models of autoregressive conditional duration}. Journal of Business and Economic Statistics 24: 104-124. 
}
\author{
Markus Belfrage
}


\seealso{
\code{\link{testTVACD}}, \code{\link{testSTACD}}.
}
\examples{
data(adjDurData)
fitModel3000obs <- acdFit(adjDurData[1:3000,])
testRmACD(fitModel3000obs, pStar = 2, robust = TRUE)
}
