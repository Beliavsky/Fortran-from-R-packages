\name{octane}
\alias{octane}
\docType{data}
\title{  Octane data }
\description{
    The octane data contains near infrared absorbance spectra (NIR) of 
    \code{n=39} gasoline samples over \code{p=226} wavelengths ranging 
    from 1102 nm to 1552 nm with measurements every two nanometers. For
    each of the 39 production gasoline samples the octane number \code{y} was measured.
    Six of the samples (25, 26, and 36-39) contain added alcohol.
}
\usage{data(octane)}
\format{
  A data frame with 39 observations and 227 columns, the wavelengts are 
    by column and the first variable is the dependent one (\code{y}).
}

%\details{}

\source{
    K.H. Esbensen, S. Schoenkopf and T. Midtgaard \emph{Multivariate Analysis in Practice}, Trondheim, Norway: Camo, 1994.
}

\references{
    M. Hubert, P. J. Rousseeuw, K. Vanden Branden (2005), ROBPCA: a new 
    approach to robust principal components analysis, \emph{Technometrics}, \bold{47}, 64--79.   
    
    P. J. Rousseeuw,  M. Debruyne, S. Engelen and M. Hubert (2006), 
    Robustness and Outlier Detection in Chemometrics, 
    \emph{Critical Reviews in Analytical Chemistry}, \bold{36}(3--4), 221--242.
}
\examples{
data(octane)

octane <- octane[, -1]    # remove the dependent variable y

pca=PcaHubert(octane, k=10)
screeplot(pca, type="lines")

pca2 <- PcaHubert(octane, k=2)
plot(pca2, id.n.sd=6)

pca7 <- PcaHubert(octane, k=7)
plot(pca7, id.n.sd=6)

}
\keyword{datasets}
