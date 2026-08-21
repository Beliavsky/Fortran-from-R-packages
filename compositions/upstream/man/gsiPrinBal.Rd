\name{gsi.PrinBal}
\alias{gsi.PrinBal}
%\alias{gsi.ilrBase}
\title{The canonical basis in the clr plane used for ilr and ipt transforms.}
\description{
  Compute the basis of a clr-plane, to use with isometric log-ratio or planar 
  transform of a (dataset of) compositions.
}
\usage{
          gsi.PrinBal( x, method="PBhclust" )
          }
\arguments{
  \item{x}{a data set (required)}
  \item{method}{method to build the principal balances; currently, one of
           "PBhclust", "PBmaxvar" or "PBangprox"}
}
\value{
  All methods give a matrix containing by columns the basis elements for the
  canonical basis of the clr-plane used for the ilr and ipt transform. Principal
  Balances are balance bases approximating Principal Component in certain ways,
  depending on the specific method chosen.
}
\details{
  Method "PBhclust" generates a Ward cluster analysis of the parts (not of
  the observations!) using as distance between parts the variation matrix. 

  Method "PBmaxvar" splits all parts in two groups after a PCA (those with
  positive weight against those with negative weight in the first PC),
  and then checks for each element whether moving it increases the variance of
  the resulting balance. The largest variance balance is chosen. That variance 
  is removed, and each of the groups of parts is analysed recursively.

  Method "angprox" calculates all possible balances (exhaustive combinations of 
  "ncol(x)" elements split in three groups: +1, 0 and -1), and computes their angles
  with the first PC. The nearest balance is chosen. Its variance is substracted, and
  the algorithm is repeated recursively ensuring orthogonality of the next 
  splittings with respect to the balances already chosen. 
  This method creates a temporary file ".APtable" that can be removed
  after execution (or copypasted in your logfiles).
  
  Note that more methods can appear in the future, and that, except "PBhclust",
  all methods can be VERY slow. Extensive reprogramming of these routines can 
  be expected in the future.
}
\seealso{
  \code{\link{ilrBase}},\code{\link{ilr}},\code{\link{clr}}
}
\references{
   Pawlowsky-Glahn, V., J.J. Egozcue and R. Tolosana-Delgado (2011) Principal 
   Balances. In: Egozcue, J.J., R. Tolosana-Delgado and M.I. Ortego (eds) 
   Proceedings of CoDaWork'2011, the 4th International Workshop on Compositional
   Data Analysis. CIMNE, Barcelona (E)

   %\url{http://congress.cimne.com/codawork11/frontal/ProgIntro.asp}
 }
\examples{
data("Hydrochem")
x = acomp(Hydrochem[,c(6:10)])
(v1 = gsi.PrinBal(x, "PBhclust"))
(v2 = gsi.PrinBal(x, "PBmaxvar"))
(v3 = gsi.PrinBal(x, "PBangprox"))
}
\author{J.J. Egozcue (programming "PBmaxvar") and R. Tolosana-Delgado (programming other methods and wrapper)}
\keyword{internal}
