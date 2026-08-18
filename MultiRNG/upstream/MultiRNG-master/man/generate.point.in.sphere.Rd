\name{generate.point.in.sphere}
\alias{generate.point.in.sphere}
\title{Point Generation for a Sphere}
\description{This function generates s points on a d-dimensional sphere.
}
\usage{
generate.point.in.sphere(no.row,d)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
}
\value{A \eqn{no.row \times d} matrix of coordinates of points in sphere.}
\references{
Marsaglia, G. (1972). Choosing a point from the surface of a sphere. \emph{Annals of Mathematical Statistics}, \bold{43}, 645-646.
}
\examples{
generate.point.in.sphere(no.row=1e5,d=3)
}
