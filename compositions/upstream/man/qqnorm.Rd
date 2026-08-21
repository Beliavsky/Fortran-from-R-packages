\name{qqnorm}
\alias{qqnorm.acomp}
\alias{qqnorm.rcomp}
\alias{qqnorm.rplus}
\alias{qqnorm.aplus}
\alias{vp.qqnorm}
\title{Normal quantile plots for compositions and amounts}
\description{
  The plots allow to check the normal distribution of multiple
  univaritate marginals by normal quantile-quantile plots. 
  For the different interpretations of amount data a different type of
  normality is assumed and checked. When an alpha-level is given the
  marginal displayed in each panel is checked for normality. 
}
\usage{
\method{qqnorm}{acomp}(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL)
\method{qqnorm}{rcomp}(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL)
\method{qqnorm}{aplus}(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL)
\method{qqnorm}{rplus}(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL)
vp.qqnorm(x,y,...,alpha=NULL)
          }
\arguments{
  \item{y}{a dataset}
  \item{fak}{a factor to split the dataset, not yet implemented in aplus
  and rplus}
  \item{panel}{the panel function to be used or a list of multiple panel
    functions}
  \item{alpha}{the alpha level of a test for normality to be performed
    for each of the displayed marginals. The levels are adjusted for
    multiple testing with a Bonferroni-correction (i.e. dividing each of
    the alpha-level by the number of test performed)}
  \item{\dots}{further graphical parameters}
  \item{x}{used by pairs only. Internal use}
}

\details{
  \code{qqnorm.rplus} and \code{qqnorm.rcomp} display qqnorm plots of
  individual amounts (on the diagonal), of pairwise differences of amounts 
  (above the diagonal) and of pairwise sums of amounts (below the
  diagonal).
  \cr
  \code{qqnorm.aplus}  displays qqnorm-plots of
  individual log-amounts (on the diagonal), of pairwise log-ratios of
  amounts (above the diagonal) and of pairwise sums of log amount (below the
  diagonal).
  \cr
  \code{qqnorm.acomp} displays qqnorm-plots of pairwise log-ratios of
  amounts in all of diagonal panels. Nothing is displayed on the
  diagonal.
  \cr
  In all cases a joint normality of the original data in the selected
  framework would imply normality in all displayed marginal
  distributions (although the reciprocal is in general not true!).
  \cr
  The marginal normality can be checked in each of the plots using a
  \code{\link[stats]{shapiro.test}}, by specifying an alpha level. The
  alpha level is corrected for multiple testing. Plots displaying a
  marginal distribution significantly deviating from a normal
  distribution at that alpha level are marked by a red exclamation mark.
  \cr
  \code{vp.qqnorm} is internally used as a panel function to make high dimensional
  plots. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}}, \code{\link{boxplot.acomp}},
  \code{\link{rnorm.acomp}}, \code{\link{rnorm.rcomp}}, 
    \code{\link{rnorm.aplus}}, \code{\link{rnorm.rplus}}
}

\examples{
data(SimulatedAmounts)
qqnorm(acomp(sa.lognormals),alpha=0.05)
qqnorm(rcomp(sa.lognormals),alpha=0.05)
qqnorm(aplus(sa.lognormals),alpha=0.05)
qqnorm(rplus(sa.lognormals),alpha=0.05)
}
\keyword{multivariate}
