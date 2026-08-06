\name{hemophilia}
\alias{hemophilia}
\docType{data}
\title{ Hemophilia Data
 }
\description{
The hemophilia data set contains two measured variables on 
75 women, belonging to two groups: n1=30 of them are 
non-carriers (normal group) and n2=45 are known hemophilia 
A carriers (obligatory carriers).
}
\usage{data(hemophilia)}
\format{
  A data frame with 75 observations on the following 3 variables.
  \describe{
    \item{\code{AHFactivity}}{AHF activity}
    \item{\code{AHFantigen}}{AHF antigen}
    \item{\code{gr}}{group - normal or obligatory carrier}
  }
}
\details{
Originally analized in the context of discriminant 
analysis by Habemma and Hermans (1974). The objective 
is to find a procedure for detecting potential hemophilia A 
carriers on the basis of two measured variables: X1=log10(AHV activity) and
X2=log10(AHV-like antigen). The first group of n1=30 women consists
of known non-carriers (normal group) and the second group of n2=45
women is selected from known hemophilia A carriers (obligatory
carriers). This data set was also analyzed by Johnson and Wichern
(1998) as well as, in the context of robust Linear Discriminant Analysis 
by Hawkins and McLachlan (1997) and Hubert and Van Driessen (2004).
}
\source{
Habemma, J.D.F, Hermans, J. and van den Broek, K. (1974) 
Stepwise Discriminant Analysis Program Using Density 
Estimation in \emph{Proceedings in Computational statistics, COMPSTAT`1974}
(Physica Verlag, Heidelberg, 1974, pp 101--110).
}

\references{
Johnson, R.A. and Wichern, D. W. \emph{Applied Multivariate
Statistical Analysis} (Prentice Hall, International
Editions, 2002, fifth edition)

Hawkins, D. M. and McLachlan, G.J. (1997) 
High-Breakdown Linear Discriminant Analysis
\emph{J. Amer. Statist. Assoc.} \bold{92} 136--143. 

Hubert, M., Van Driessen, K. (2004) Fast and robust discriminant analysis,
\emph{Computational Statistics and Data Analysis}, \bold{45}
301--320.

}

\examples{
data(hemophilia)
plot(AHFantigen~AHFactivity, data=hemophilia, col=as.numeric(as.factor(gr))+1)
##
## Compute robust location and covariance matrix and 
## plot the tolerance ellipses
(mcd <- CovMcd(hemophilia[,1:2]))
col <- ifelse(hemophilia$gr == "carrier", 2, 3) ## define clours for the groups
plot(mcd, which="tolEllipsePlot", class=TRUE, col=col)

%x <- hemophilia[,-3]
%gr <- hemophilia[,3]
%Linda(x, gr, method="ogk")
}
\keyword{datasets}
