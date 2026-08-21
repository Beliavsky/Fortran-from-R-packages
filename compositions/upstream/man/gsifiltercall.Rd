\name{gsifiltercall}
\alias{gsi.filtercall}
\title{Calling from a function with the own parameters}
\description{
  gsi.filtercall gives those arguments of the arguments of the calling
  function to called function that are made for this function. It is possible
  to modify the parameters on the flight.
}
\usage{
gsi.filtercall(fkt,...,overwrite=NULL,
               transmit=character(0),default=list(...),
               prefix=if(is.character(fkt)) fkt else NULL)
}
\arguments{
  \item{fkt}{a function or a name of a function to be called}
  \item{\dots}{default arguments, which are used if not overwritten by own
  arguments}
  \item{overwrite}{a named list overwriting arguments provided in any other way}
  \item{transmit}{a character vector containing all named arguments that
  should be given to \code{fkt}, although they are not formal parameters of
  \code{fkt}}
  \item{default}{the list of default parameters to be used if not otherwise
  specified}	
  \item{prefix}{the user can specify arguments for \code{fkt}
  explicitly by prefixing them with \code{prefix}.. Typically the prefix is the
  name of the function}
}
\value{
  The value of a call to \code{fkt} with those parameters modified by 
  \code{gsi.filtercall} itself
}
\details{
 gsi.filtercall is a kind of R technology help. It allows to pass parameters in a
 \dots{} fashion with less effort, more robustness and more control, since one
 does not need to copy all explicitly-given parameters. It simply manages the
 \dots{}-parameters accepted by \code{fkt}, allowing for giving default
 values, or overwriting specific parameters without need to mention them in
 the formal arguments. This also helps to clarify the documentation of \code{fkt}.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package. However I would like to propose a function like
  this to become standard technology}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi}}
}
\examples{
#	fxy <- function(x,y) x+y
#	fxz <- function(x,z) x*z
#	fxyz <- function(...) {
#	   list(gsi.filtercall("fxy"),gsi.filtercall("fxz"))
#        }	 
#	fxyz(x=1,y=2,z=7)
#	fxyz(x=1,fxz.x=2,y=10,z=55)
}
\keyword{internal}
