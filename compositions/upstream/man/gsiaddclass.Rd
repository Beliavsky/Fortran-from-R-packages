\name{gsi.addclass}
\alias{gsi.addclass}
\title{Internal function: give a derived subclass to an object}
\description{
  This function just extends the class of an object by the given class.
}
\usage{
gsi.addclass( x , cls)
}
\arguments{
  \item{x}{the object}
  \item{cls}{the new additional class}
}
\value{
  The object x with additional class attached.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\examples{
x=rplus(1:10)
x
#gsi.addclass(x,"goofy")
}
\keyword{internal}
