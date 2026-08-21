\name{RobustWeights}
\Rdversion{2.1}
\alias{RobustWeights}
\alias{robwts}
\title{
Robustification of the sampling weights
}
\description{
Calculation of a Huber-type correction factor by which the vector of weights is multiplied.
}
\usage{
robwts(x, w=rep(1,length(x)), c=0.01, alpha=0.001)
}
\arguments{
  \item{x}{numeric; vector of data values.}
  \item{w}{numeric; vector of weights. Must have the same length as \code{x}. By default \code{w} is a vector of 1.}
  \item{c}{numeric; a constant which can take different values, e.g. 0.01,0.02. By default \eqn{c=0.1}.}
  \item{alpha}{numeric; a probability in the interval \eqn{(0,1)}. By default \eqn{alpha=0.001}.}
}
\details{
If \eqn{x} denotes the observed value and \eqn{x_{\alpha}} the \eqn{\alpha}-th qiantile of the Fisk distribution, then we define our scale as:
\deqn{d = \displaystyle \frac{x_{1-\alpha}}{b} - \frac{x_{\alpha}}{b}}. Next, the correction factor is calculated as follows:
\deqn{corr = \max\left\{c, \min\left(1,\displaystyle \frac{d}{|b/x-1|},\frac{d}{|x/b-1|}\right)\right\}}
}
\value{\code{robwts} returns a list of two elements: the vector of correction factors by which the weights are multiplied and the vector of corrected (robustified) weights. 
}
\references{
 Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
}
\author{
Monique Graf
}
