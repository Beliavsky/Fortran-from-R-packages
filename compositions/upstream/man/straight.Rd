\name{straight}
\alias{straight}
\alias{straight.rmult}
\alias{straight.acomp}
\alias{straight.rcomp}
\alias{straight.aplus}
\alias{straight.rplus}
\title{Draws straight lines.}
\description{
  The function draws lines in a given direction \code{d} through points \code{x}.
}
\usage{
          straight(x,...)
          \method{straight}{acomp}(x,d,...,steps=30,aspanel=FALSE)
          \method{straight}{rcomp}(x,d,...,steps=30,aspanel=FALSE)
          \method{straight}{aplus}(x,d,...,steps=30,aspanel=FALSE)
          \method{straight}{rplus}(x,d,...,steps=30,aspanel=FALSE)
          \method{straight}{rmult}(x,d,...,steps=30,aspanel=FALSE)
          }
\arguments{
  \item{x}{dataset of points of the given type to draw the line through}
  \item{d}{dataset of directions of the line}
  \item{\dots}{further graphical parameters}
  \item{steps}{the number of discretisation points to draw the segments, 
    since the representation might not visually be a straight line}
  \item{aspanel}{Logical, indicates use as slave to do acutal drawing only.}
}
\details{
The functions add lines to the graphics generated with the corresponding
plot functions.
\cr
Adding to multipaneled plots redraws the plot completely, and is only
possible when the plot has been created with the plotting routines from
this library.
\cr
Lines end when they leave the space (e.g. the simplex), which sometimes
leads to the impression of premature end (specially in \code{\link{rcomp}} 
geometry). 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}},\code{\link{lines.acomp}} 
}

\examples{
data(SimulatedAmounts)

plot(acomp(sa.lognormals))
straight(mean(acomp(sa.lognormals)),
         princomp(acomp(sa.lognormals))$Loadings[1,],
         col="red")
straight(mean(rcomp(sa.lognormals)),
         princomp(rcomp(sa.lognormals))$loadings[,1],
         col="blue")

plot(aplus(sa.lognormals[,1:2]))
straight(mean(aplus(sa.lognormals[,1:2])),
         princomp(aplus(sa.lognormals[,1:2]))$Loadings[1,],
         col="red")
straight(mean(rplus(sa.lognormals[,1:2])),
         princomp(rplus(sa.lognormals[,1:2]))$loadings[,1],
         col="blue")

plot(rplus(sa.lognormals[,1:2]))
straight(mean(aplus(sa.lognormals[,1:2])),
         princomp(aplus(sa.lognormals[,1:2]))$Loadings[1,],
         col="red")
straight(mean(rplus(sa.lognormals[,1:2])),
         princomp(rplus(sa.lognormals[,1:2]))$loadings[,1],
         col="blue")

}

\keyword{multivariate}
