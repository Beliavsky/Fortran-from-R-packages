\name{Minkdiet}
\docType{data}
\alias{Minkdiet}
\title{
  An underdetermined linear inverse problem: estimating diet composition of
  Southeast Alaskan Mink.
}
\description{
  Input data for assessing the diet composition of mink in southeast Alaska,
  using C and N isotope ratios (d13C and d15N).

  The data consist of
  \enumerate{
    \item the input matrix \code{Prey}, which contains the C (1st row) and N
      (2nd row) isotopic values of the prey items (columns), corrected for
      fractionation.
    \item the input vector \code{Mink}, with the C and N isotopic value of
      the predator, mink
  }
  There are seven prey items as food sources:
  \itemize{
    \item fish
    \item mussels
    \item crabs
    \item shrimp
    \item rodents
    \item amphipods
    \item ducks
  }

  The d13C and d15N for each of these prey items, and for mink (the predator)
  was assessed. The isotopic values of the preys were corrected for
  fractionation.

  The problem is to find the diet composition of mink, e.g. the fraction of
  each of these food items in the diet.

  Mathematically this is by solving an lsei (least squares with equalities
  and inequalities) problem: \eqn{Ex=f} subject to \eqn{Gx>h}.

  The equalities \eqn{Ex=f}:
   \deqn{d13CMink = p1*d13Cfish+p2*d13Cmussels + .... + p7*d13Cducks}
   \deqn{d15NMink = p1*d15Nfish+p2*d15Nmussels + .... + p7*d15Nducks}
   \deqn{1 = p1+p2+p3+p4+p5+p6+p7}

  and inequalities \eqn{Gx>h}:
  \deqn{pi >= 0}
 
  are solved for p1,p2,...p7.
  

  The first two equations calculate the isotopic ratio of the consumer (Mink)
  as a weighted average of the ratio of the food sources

  Equation 3 assures that the sum of all fraction equals 1.

  As there are 7 unknowns and only 3 equations, the model is UNDERdetermined,
  i.e. there exist an infinite amount of solutions.

  This model can be solved by various techniques:
  \enumerate{
    \item least distance programming will select the "simplest" solution.
      See \code{\link{ldei}}.
    \item the remaining uncertainty ranges of the fractions can be estimated
      using linear programming. See \code{\link{xranges}}
    \item the statistical distribution of the fractions can be estimated using
      an MCMC algorithm which takes a sample of the solution space.
      See \code{\link{xsample}}
  }
}

\usage{Minkdiet}

\format{
  a list with matrix \code{Prey} and vector \code{Mink}.
  \itemize{
    \item \code{Prey} contains the isotopic composition (13C and 15N) of the
      7 possible food items of Mink
    \item \code{Mink} contains the isotopic composition (13C and 15N) of Mink
  }

  columnnames of \code{Prey} are the food items, rownames of \code{Prey}
  (=names of Mink) are the names of the isotopic elements.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
# 1. visualisation of the data
plot(t(Minkdiet$Prey), xlim = c(-25, -13), xlab = "d13C", ylab = "d15N",
     main = "Minkdiet", sub = "Ben-David et al. (1979)")

text(t(Minkdiet$Prey)-0.1, colnames(Minkdiet$Prey))   

points(t(Minkdiet$Mink), pch = 16, cex = 2)
text(t(Minkdiet$Mink)-0.15, "MINK", cex = 1.2)   
legend("bottomright", pt.cex = c(1, 2), pch = c(1, 16),
       c("food", "predator"))

# 2. Generate the food web model input matrices
# the equalities: 
E <- rbind(Minkdiet$Prey, rep(1, 7))
F <- c(Minkdiet$Mink, 1)

# the inequalities (all pi>0)
G <- diag(7)
H <- rep(0, 7)

# 3. Select the parsimonious (simplest) solution
parsimonious <- ldei(E, F, G = G, H = H)

# 4. show results
data.frame(food = colnames(Minkdiet$Prey),
           fraction = parsimonious$X)

dotchart(x = as.vector(parsimonious$X), labels = colnames(Minkdiet$A),
         main = "Estimated diet composition of Mink",
         sub = "using ldei and xranges", pch = 16)

# 5. Ranges of diet composition
iso   <- xranges(E, F, ispos = TRUE)
segments(iso[,1], 1:ncol(E), iso[,2], 1:ncol(E))
legend  ("topright", pch = c(16, NA), lty = c(NA, 1),
          legend = c("parsimonious", "range"))

pairs (xsample(E = E, F = F, G = diag(7), H = rep(0, 7), iter = 1000)$X,
       main = "Minkdiet 1000 solutions, using xsample")
}

\references{
  Ben-David M, Hanley TA, Klein DR, Schell DM (1997) Seasonal changes
  in diets of coastal and riverine mink: the role of spawning Pacific salmon.
  Canadian Journal of Zoology 75:803-811.
}

\seealso{
  \code{\link{ldei}} to solve for the parsimonious solution

  \code{\link{xranges}} to solve for the uncertainty ranges

  \code{\link{xsample}} to sample the solution space
}
\keyword{datasets}

