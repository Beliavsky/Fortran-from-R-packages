\name{plotEliteDistrib}
\alias{plotEliteDistrib}

\title{
Plot the distribution of elite members
}
\description{
A simple function to create distribution plots of the elite members after the optimization procedure. The distribution is graphically represented using a kernel density estimation. Aditionally, this function also indicates the best and average members.
}
\usage{
plotEliteDistrib(elite)
}

\arguments{
  \item{elite}{
  A matrix containing parameters of the elite members.
}
}
\value{
A graphic representation of the elite members and also of the best and average members. 
}
\author{
Alberto Krone-Martins
}

\seealso{
\code{\link{ceimOpt}}
}
\examples{
# Solve a simple 2D problem and show the distribution of the parameters
po <- ceimOpt(OptimFunction=function(x){(x[1]+1)^2+(x[2]+2)^2}, maxIter=100, 
      epsilon=0.1, nParams=2)
plotEliteDistrib(po$EliteMembers)
rm(po)
}

\keyword{ misc }
