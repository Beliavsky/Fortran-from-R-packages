\name{firstValueRow}
\alias{firstValueRow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returning the first >0 value in each row of a matrix
}
\description{
Technical function implemented in the pricing functions of the package. It returns the first >0 value in each row of a matrix and assign zero to all subsequent values.
}
\usage{
firstValueRow(x)
}

\arguments{
  \item{x}{
A matrix.
}
}

\value{
A matrix.
}
\author{
Mikhail A. Beketov
}
\seealso{
Functions: \code{\link{AmerPutLSM}},  
\code{\link{AmerPutLSM_AV}}, 
\code{\link{AmerPutLSM_CV}}, 
\code{\link{AsianAmerPutLSM}}, 
\code{\link{QuantoAmerPutLSM}}, and 
\code{\link{QuantoAmerPutLSM_AV}}.
}

\examples{
mat<-matrix(c(0,0,2,0,4,0,3,0,1,9,8,7), ncol=4)
mat
firstValueRow(mat)
}

\keyword{ Monte Carlo }
\keyword{ Option pricing }
\keyword{ American put }
