\name{readMps}
\alias{readMps}

\title{Read MPS Files}

\description{
  This function reads MPS files - the standard format for Linear Programming problems.
}

\usage{
   readMps( file, solve=FALSE, maximum=FALSE )
}

\arguments{
  \item{file}{a character string naming the file to read.}
  \item{solve}{logical. Should the problem be solved after reading it
               from the file (using \code{\link{solveLP}})?}
  \item{maximum}{logical. Should we maximize or minimize (the default)?}
}

\details{
   Equality constraints and 'greater than'-bounds are not implemented yet.
}

\value{
  \code{readMps} returns a list containing following objects:

  \item{name}{the name of the Linear Programming problem.}
  \item{cvec}{vector \eqn{c}.}
  \item{bvec}{vector \eqn{b}.}
  \item{Amat}{matrix \eqn{A}.}
  \item{res}{if \code{solve} is TRUE, it contains the results of the solving
              process (an object of class \code{\link{solveLP}}).}
}

\author{
   Arne Henningsen
}

\seealso{
   \code{\link{solveLP}}, \code{\link{writeMps}}
}


\examples{

## example of Steinhauser, Langbehn and Peters (1992)
## Production activities
cvec <- c(1800, 600, 600)  # gross margins
names(cvec) <- c("Cows","Bulls","Pigs")

## Constraints (quasi-fix factors)
bvec <- c(40, 90, 2500)  # endowment
names(bvec) <- c("Land","Stable","Labor")

## Needs of Production activities
Amat <- rbind( c(  0.7,   0.35,   0 ),
               c(  1.5,   1,      3 ),
               c( 50,    12.5,   20 ) )

## Write to MPS file
writeMps( "steinh.mps", cvec, bvec, Amat, "Steinhauser" )

## delete all LP objects
rm( cvec, bvec, Amat )

## Read LP data from MPS file and solve it.
lp <- readMps( "steinh.mps", TRUE, TRUE )

## Print the results
lp$res

## remove the MPS file
file.remove( "steinh.mps" )
}

\keyword{ optimize }
