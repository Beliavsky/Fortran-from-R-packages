\name{apca}
\alias{apca}
\title{
Asymptotic PCA
}
\description{
  Asymptotic Principal Components Analysis for a fixed number of factors
}
\usage{
apca(x,nf)
}
\arguments{
  \item{x}{
    a numeric matrix or other object for which 'as.matrix' will produce
    a numeric matrix.
  }
  \item{nf}{
    number of factors desired.
  }
}
\details{
  NOTE:  This is a preliminary version of this function, and it may be
  modified in the future.
}
\value{
  A list with four components:

  \item{eig}{eigenvalues}
  \item{factors}{estimated factor scores}
  \item{loadings}{estimated factor loadings}
  \item{rsq}{
    R-squared from the regression of each variable on the factor space
  }
}
\author{
  Ruey Tsay
}
\references{
  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley,
  sec. 9.6, pp. 436-440)
}
\seealso{
  \code{\link[stats]{princomp}}
}
\examples{
# Consider the monthly simple returns of 40 stocks on NYSE and NASDAQ
# from 2001 to 2003 with 36 observations.
data(m.apca0103)
dim(m.apca0103)
M.apca0103 <- with(m.apca0103,
    array(return, dim = c(36, 40),
          dimnames = list(as.character(date[1:36]),
                          paste("Co", CompanyID[seq(1, 1440, 36)], sep = ""))))

# The traditional PCA is not applicable to estimate the factor model
# because of the singularity of the covariance matrix.  The asymptotic
# PCA provides an approach to estimate factor model based on asymptotic
# properties. For the simple example considered, the sample size is
# $T$ = 36 and the dimension is $k$ = 40. If the number of factor is
# assumed to be 1, the APCA gives a summary of the factor loadings as
# below:
#
apca40 <- apca(M.apca0103, 1)
#
# (min, 1st Quartile, median, mean, 3rd quartile, max) =
# (0.069, 0.432, 0.629, 0.688, 1.071, 1.612).
#
# Note that the sign of any loading vector is not uniquely determined
# in the same way as the sign of an eigenvector is not uniquely
# determined. The output also contains the summary statistics of the
# R-squares of individual returns, i.e. the R-squares measuring the
# total variation of individual return explained by the factors. For
# the simple case considered, the summary of R-squares is (min, 1st
# Quartile, median, mean, 3rd quartile, max) =
# (0.090,0.287,0.487,0.456,0.574,0.831).
}
\keyword{ts}
