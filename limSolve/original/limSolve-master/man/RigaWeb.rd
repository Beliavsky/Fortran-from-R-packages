\name{RigaWeb}
\docType{data}
\alias{RigaWeb}
\title{
  An underdetermined linear inverse problem: the Gulf of Riga *spring*
  planktonic food web
}
\description{
  Input matrices and vectors for estimating the flows in the planktonic food
  web of the Gulf of Riga.

  (as in Donali et al. (1999)).

  The original input file can be found in the package subdirectory
  \code{/inst/docs/RigaSpring.input}

  There are 7 functional compartments:
  P1,P2,B,N,Z,D,OC
  (two phytoplankton groups, Bacteria, Nanozooplankton, Zooplankton,
  Detritus and DOC).

  and 2 externals: CO2 and SED (sedimentation)

  These are connected with 26 flows:
  P1->CO2, P2->CO2, Z->CO2, N->CO2, B->CO2, CO2->P1, CO2->P2, P1->Z, P1->N,
  P1->DOC, P1->SED, P2->DOC, P2->Z, P2->D, P2->SED, N->DOC, N->Z, Z->DOC,
  Z->D, Z->SED, D->Z, D->DOC, D->SED, B->N, B->SED, DOC->B

  The lsei model contains:
  \itemize{
    \item 14 equalities (Ax=B): the 7 mass balances (one for each compartment)
      and 7 measurement equations
    \item 26 unknowns (x), the flow values
    \item 45 inequalities (Gx>h). The first 19 inequalities impose bounds
      on some combinations of flows.
      The last 26 inequalities impose that the flows have to be positive.
  }
  As there are more unknowns (26) than equations (14), there exist an
  infinite amount of solutions (it is an underdetermined problem).
}
\usage{RigaWeb}
\format{
  A list with the matrices and vectors that constitute the mass balance problem:
  \code{A}, \code{B}, \code{G} and \code{H}.

  The columnames of \code{A} and \code{G} are the names of the unknown
  reaction rates;
  The first 14 rownames of \code{A} give the names of the components
  (these rows consitute the mass balance equations).

}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{

E <- RigaWeb$A
F <- RigaWeb$B
G <- RigaWeb$G
H <- RigaWeb$H

# 1. parsimonious (simplest) solution
pars <- lsei(E = E, F = F, G = G, H = H)$X

# 2.ranges of all unknowns, including the central value
xr   <- xranges(E = E, F = F, G = G, H = H, central = TRUE)

# the central point is a valid solution:
X <- xr[,"central"]
max(abs(E\%*\%X - F))
min(G\%*\%X - H)

\dontrun{   # this does not work on windows i386!
# 3. Sample solution space; the central value is a good starting point
#   for algorithms cda and rda - but these need many iterations   
xs  <- xsample(E = E, F = F, G = G, H = H,
               iter = 10000, out = 1000, type = "rda", x0 = X)$X
# better convergence using 50000 iterations, but this takes a while
xs  <- xsample(E = E, F = F, G = G, H = H,
               iter = 50000, out = 1000, type = "rda", x0 = X)$X

pairs(xs, pch = ".", cex = 2, gap = 0, upper.panel = NULL)


# using mirror algorithm takes less iterations,
# but an iteration takes more time ; it is better to start in a corner...
# (i.e. no need to use X as starting value)
xs  <- xsample(E = E, F = F, G = G, H = H,
               iter = 1500, output = 500, type = "mirror")$X
pairs(xs, pch = ".", cex = 2, gap = 0, upper.panel = NULL, 
      yaxt = "n", xaxt = "n")

# Print results:
data.frame(pars = pars, xr[ ,1:2], Mean = colMeans(xs), sd = apply(xs, 2, sd))
}
}
\references{
  Donali, E., Olli, K., Heiskanen, A.S., Andersen, T., 1999. Carbon flow
  patterns in the planktonic food web of the Gulf of Riga, the Baltic Sea:
  a reconstruction by the inverse method. Journal of Marine Systems 23,
  251..268.
}

\keyword{datasets}

