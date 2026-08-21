\name{kingTetrahedron}
\alias{kingTetrahedron}
\title{Ploting composition into rotable tetrahedron}
\description{
  Plots acomp/rcomp objects into tetrahedron exported in kinemage format.
  }
\usage{
   kingTetrahedron(X, parts=1:4, file="tmptetrahedron.kin",
   clu=NULL,vec=NULL, king=TRUE, scale=0.2, col=1, 
   title="Compositional Tetrahedron")
          }
\arguments{
  \item{X}{a compositional acomp or rcomp object of 4 or more parts}
  \item{parts}{a numeric or character vector specifying the 4 parts to
    be used.}
  \item{file}{
       file.kin for 3D display  with the KiNG (Kinemage, Next Generation) interactive
       system for three-dimensional vector graphics.
      }
  \item{clu}{partition determining the colors of points}
  \item{vec}{vector of values determining points sizes}
  \item{king}{FALSE for Mage; TRUE for King (described below)}
  \item{scale}{relative size of points}
  \item{col}{color of points if clu=NULL}
  \item{title}{The title of the plot}
  }
  \value{
    The function is called for its side effect of generating 
  a file  for 3D display  with the KiNG (Kinemage, Next Generation) interactive
       system for three-dimensional vector graphics.
       Works only with KiNG viewer. More information (and links to the actual viewers) 
       can be found at \url{https://en.wikipedia.org/wiki/Kinemage}
      }
\details{
  The  routine transforms a 4 parts mixture m quadrays into 3-dimensional
  XYZ coordinates and writes them as file.kin. For this transformation we
  apply K. Urner: Quadrays and XYZ at
  \url{http://www.grunch.net/synergetics/quadxyz.html}. 
  The kin file we display as 3-D animation with KiNG viewer.  
  A kinemage is a dynamic, 3-D illustration. 
  The best way to take advantage of that is by rotating it and twisting it around 
  with the mouse click near the center of the graphics window and slowly draging right or left,
  up or down. Furthermore by clicking on points with the mouse (left button again), 
  the label associated with each point will appear in the bottom left of the graphics area
  and also the distance from this point to the last will be displayed.
  With the right button drag we can zoom in and out of the picture.
  This animation supports coloring and different sizing of points.
  
  We can display the kin file as 3-D animation also with MAGE viewer a previous version of KiNG,
  %also a free software available at http://kinemage.biochem.duke.edu.
  more information (and links to the software) can be found at 
  \url{https://en.wikipedia.org/wiki/Kinemage}.
  For this one has to put king=FALSE as a parameter.
}
\note{
  This routine and the documentation is based on mix.Quad2net from the
  MixeR-package of Vladimir Batagelj and Matevz Bren, and has been
  contributed by Matevz Bren to this package. Modifications
  have been applied to make function compatible with the philosophy and
  objects of the compositions package.
}
\seealso{
  \code{\link{plot.acomp}}
}
\author{Vladimir Batagelj and Matevz Bren, with slight modifications of K.Gerald van den Boogaart}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

  %{http://vlado.fmf.uni-lj.si/pub/MixeR/}
  
  \url{http://www.grunch.net/synergetics/quadxyz.html}\cr


}
\examples{
\dontrun{
data(SimulatedAmounts)
dat <- acomp(sa.groups5)
hc <- hclust(dist(dat), method = "complete")  # data are clustered 
kingTetrahedron(dat,parts=1:4, file="myfirst.kin", clu=cutree(hc,7),
scale=0.2)
# the 3-D plot is written into Glac1.kin file to be displayed with KiNG viewer. 
# The seven clusters partition is notated with different colors of points.
}
}
\keyword{dynamic}
