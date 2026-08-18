\name{draw.correlated.binary}
\alias{draw.correlated.binary}
\title{Generation of Correlated Binary Data}
\description{This function implements pseudo-random number generation for a multivariate Bernoulli distribution (correlated binary data).
}
\usage{
draw.correlated.binary(no.row,d,prop.vec,corr.mat)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{prop.vec}{Vector of means.}
  \item{corr.mat}{Correlation matrix.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\references{
Park, C. G., Park, T., & Shin D. W. (1996). A simple method for generating correlated binary variates. \emph{The American Statistician}, \bold{50(4)}, 306-310.
}
\seealso{\code{\link{loc.min}}}
\examples{
cmat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
propvec=c(0.3,0.5,0.7)
\donttest{
mydata=draw.correlated.binary(no.row=1e5,d=3,prop.vec=propvec,corr.mat=cmat)
apply(mydata,2,mean)-propvec
cor(mydata)-cmat
}
}

