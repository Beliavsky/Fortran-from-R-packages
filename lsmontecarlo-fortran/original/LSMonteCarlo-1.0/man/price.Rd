\name{price}
\alias{price}

\title{
Extracting price from the pricing functions outputs
}
\description{
The function is nothing else, but the \code{object$price} action, with the \code{object} returned by the pricing functions in the package.
}
\usage{
price(x)
}

\arguments{
  \item{x}{
Object returned by the pricing functions in the package
}
}

\value{
The function returns the price as a single number (class "numeric").
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
put<-AmerPutLSM()
price(put)
put$price
}

\keyword{ Monte Carlo }
\keyword{ Option pricing }
\keyword{ American put }
