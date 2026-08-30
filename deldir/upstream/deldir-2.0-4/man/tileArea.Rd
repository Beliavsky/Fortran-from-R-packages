\name{tileArea}
\alias{tileArea}
\title{
    Area of a Dirichlet tile.
}
\description{
    Calculates the area of a Dirichlet tile, applying a discrete
    version of Stoke's theorem.
}
\usage{
tileArea(x, y, rw)
}
\arguments{
  \item{x}{
  The \code{x}-coordinates of the vertices of the tile, in
  \bold{anticlockwise} direction.   The last coordinate should
  \bold{not} repeat the first.
}
  \item{y}{
  The \code{y}-coordinates of the vertices of the tile, in
  \bold{anticlockwise} direction.   The last coordinate should
  \bold{not} repeat the first.
}
  \item{rw}{
  A vector of length 4 specifying the rectangular window in
  which the relevant tessellation was constructed.  See
  \code{\link{deldir}()} for more detail.  Actually this
  can be any rectangle containing the tile in question.
}
}
\details{
  The heavy lifting is done by the Fortran subroutine \code{stoke()}
  which is called by the \code{.Fortran()} function.
}
\value{
  A positive scalar.
}
\author{\rolf}

\seealso{
  \code{\link{deldir}()}
  \code{\link{tilePerim}()}
}
\examples{
set.seed(42)
x <- runif(20)
y <- runif(20)
z <- deldir(x,y,rw=c(0,1,0,1))
w <- tile.list(z)
with(w[[1]],tileArea(x,y,rw=z$rw))
sapply(w,function(x,rw){tileArea(x$x,x$y,attr(w,"rw"))})
x <- c(0.613102,0.429294,0.386023,0.271880,0.387249,0.455900,0.486101)
y <- c(0.531978,0.609665,0.597780,0.421738,0.270596,0.262953,0.271532)
# The vertices of the Dirichlet tile for point 6.
tileArea(x,y,rw=c(0,1,0,1))
tileArea(x,y,rw=c(-1,2,-3,4)) # Same as above.
}

\keyword{ utilities }
