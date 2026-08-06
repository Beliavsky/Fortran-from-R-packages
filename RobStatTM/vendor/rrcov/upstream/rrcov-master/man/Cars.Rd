\name{Cars}
\alias{Cars}
\docType{data}
\title{
    Consumer reports car data: dimensions
}
\description{
A data frame containing 11 variables with different dimensions of 111 cars
}
\usage{data(Cars)}
\format{
  A data frame with 111 observations on the following 11 variables.
  \describe{
    \item{\code{length}}{a numeric vector}
    \item{\code{wheelbase}}{a numeric vector}
    \item{\code{width}}{a numeric vector}
    \item{\code{height}}{a numeric vector}
    \item{\code{front.hd}}{a numeric vector}
    \item{\code{rear.hd}}{a numeric vector}
    \item{\code{front.leg}}{a numeric vector}
    \item{\code{rear.seating}}{a numeric vector}
    \item{\code{front.shoulder}}{a numeric vector}
    \item{\code{rear.shoulder}}{a numeric vector}
    \item{\code{luggage}}{a numeric vector}
  }
}

%\details{}
\source{
    Consumer reports. (April 1990). http://backissues.com/issue/Consumer-Reports-April-1990, pp. 235--288.
}
    
\references{
    Chambers, J. M. and Hastie, T. J. (1992). Statistical models in S. Cole, Pacific Grove, CA:
    Wadsworth and Brooks, pp. 46--47.
    
    M. Hubert, P. J. Rousseeuw, K. Vanden Branden (2005), ROBPCA: A new approach to robust 
    principal components analysis, \emph{Technometrics}, 
    \bold{47}, 64--79.    
}

\examples{
    data(Cars)

## Plot a pairwise scaterplot matrix
    pairs(Cars[,1:6])

    mcd <- CovMcd(Cars[,1:6])    
    plot(mcd, which="pairs")
    
## Start with robust PCA
    pca <- PcaHubert(Cars, k=ncol(Cars), kmax=ncol(Cars))
    pca

## Compare with the classical PCA
    prcomp(Cars)

## or  
    PcaClassic(Cars, k=ncol(Cars), kmax=ncol(Cars))
    
## If you want to print the scores too, use
    print(pca, print.x=TRUE)

## Using the formula interface
    PcaHubert(~., data=Cars, k=ncol(Cars), kmax=ncol(Cars))

## To plot the results:

    plot(pca)                    # distance plot
    pca2 <- PcaHubert(Cars, k=4)  
    plot(pca2)                   # PCA diagnostic plot (or outlier map)
    
## Use the standard plots available for prcomp and princomp
    screeplot(pca)    # it is interesting with all variables    
    biplot(pca)       # for biplot we need more than one PCs
    
## Restore the covraiance matrix     
    py <- PcaHubert(Cars, k=ncol(Cars), kmax=ncol(Cars))
    cov.1 <- py@loadings \%*\% diag(py@eigenvalues) \%*\% t(py@loadings)
    cov.1      

}
\keyword{datasets}
