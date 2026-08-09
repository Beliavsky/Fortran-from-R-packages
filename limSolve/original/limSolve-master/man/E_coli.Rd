\name{E_coli}
\docType{data}
\alias{E_coli}
\title{
  An underdetermined linear inverse problem: the Escherichia Coli
  Core Metabolism Model.
}

\description{
  Input matrices and vectors for performing Flux Balance Analysis of the
  E.coli metabolism

  (as from \code{http://gcrg.ucsd.edu/Downloads/Flux_Balance_Analysis}).

  The original input file can be found in the package subdirectory
  \code{/inst/docs/E_coli.input}

  There are 53 substances:

  GLC, G6P, F6P, FDP, T3P2, T3P1, 13PDG, 3PG, 2PG, PEP, PYR, ACCOA, CIT, ICIT,
  AKG, SUCCOA, SUCC, FUM, MAL, OA, ACTP, ETH, AC, LAC, FOR, D6PGL, D6PGC, RL5P,
  X5P, R5P, S7P, E4P, RIB, GLX, NAD, NADH, NADP, NADPH, HEXT, Q, FAD, FADH,
  AMP, ADP, ATP, GL3P, CO2, PI, PPI, O2, COA, GL, QH2

  and 13 externals:

  Biomass, GLCxt, GLxt, RIBxt, ACxt, LACxt, FORxt, ETHxt, SUCCxt, PYRxt,
  PIxt, O2xt, CO2xt


  There are 70 unknown reactions (named by the gene encoding for it):

  GLK1,  PGI1,  PFKA,  FBP, FBA, TPIA, GAPA, PGK, GPMA, ENO, PPSA, PYKA,
  ACEE, ZWF, PGL, GND, RPIA, RPE, TKTA1, TKTA2, TALA, GLTA, ACNA, ICDA,
  SUCA, SUCC1, SDHA1, FRDA, FUMA, MDH, DLD1, ADHE2, PFLA, PTA, ACKA, ACS,
  PCKA, PPC, MAEB, SFCA, ACEA, ACEB, PPA, GLPK, GPSA1, RBSK, NUOA, FDOH,
  GLPD, CYOA, SDHA2, PNT1A, PNT2A, ATPA, GLCUP, GLCPTS, GLUP, RIBUP, ACUP,
  LACUP, FORUP, ETHUP, SUCCUP, PYRUP, PIUP, O2TX, CO2TX, ATPM, ADK, Growth

  The \code{lsei} model contains:
  \itemize{
    \item 54 equalities (Ax=B): the 53 mass balances (one for each substance)
      and one equation that sets the ATP drain flux for constant maintenance
      requirements to a fixed value (5.87)
   \item 70 unknowns (x), the reaction rates
   \item 62 inequalities (Gx>h). The first 28 inequalities impose bounds
     on some reactions.
     The last 34 inequalities impose that the reaction rates have to be
     positive (for unidirectional reactions only).
   \item 1 function that has to be maximised, the biomass production (growth).
  }
  As there are more unknowns (70) than equations (54), there exist an
  infinite amount of solutions (it is an underdetermined problem).
}

\usage{E_coli}

\format{
  A list with the matrices and vectors that constitute the mass balance problem:
  \code{A}, \code{B}, \code{G} and \code{H} and

  \code{Maximise}, with the function to maximise.

  The columnames of \code{A} and \code{G} are the names of the unknown
  reaction rates;
  The first 53 rownames of \code{A} give the names of the components
  (these rows consitute the mass balance equations).
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
\dontrun{
# 1. parsimonious (simplest) solution
pars <- lsei(E = E_coli$A, F = E_coli$B, 
             G = E_coli$G, H = E_coli$H)$X

# 2. the optimal solution - solved with linear programming
#    some unknowns can be negative

LP <- linp(E = E_coli$A, F = E_coli$B,G = E_coli$G, H = E_coli$H,
           Cost = -E_coli$Maximise, ispos = FALSE)
(Optimal <- LP$X)

# 3.ranges of all unknowns, including the central value and all solutions
xr   <- xranges(E = E_coli$A, F = E_coli$B, G = E_coli$G, H = E_coli$H,
                central = TRUE, full = TRUE)

# the central point is a valid solution:
X <- xr[ ,"central"]

max(abs(E_coli$A\%*\%X - E_coli$B))
min(E_coli$G\%*\%X - E_coli$H)
}
# 4. Sample solution space; the central value is a good starting point
#   for algorithms cda and rda - but these need many iterations
\dontrun{
xs <- xsample(E = E_coli$A, F = E_coli$B, G = E_coli$G,H = E_coli$H,
              iter = 50000, out = 5000, type = "rda", x0 = X)$X
pairs(xs[ ,10:20], pch = ".", cex = 2, main = "sampling, using rda")
}

# using mirror algorithm takes less iterations,
# but an iteration takes more time ; it is better to start in a corner...
# (i.e. no need to use X as starting value)
\dontrun{
xs <- xsample(E = E_coli$A, F = E_coli$B, G = E_coli$G, H = E_coli$H,
              iter = 2000, out = 500, jmp = 50, type = "mirror",
              x0 = X)$X
pairs(xs[ ,10:20], pch = ".", cex = 2, main = "sampling, using mirror")

# Print results:
data.frame(pars = pars, Optimal = Optimal, xr[ ,1:2],
           Mean = colMeans(xs), sd = apply(xs, 2, sd))

# Plot results
par(mfrow = c(1, 2))
nr <- length(Optimal)/2

ii <- 1:nr
dotchart(Optimal[ii], xlim = range(xr), pch = 16)
segments(xr[ii,1], 1:nr, xr[ii,2], 1:nr)

ii <- (nr+1):length(Optimal)
dotchart(Optimal[ii], xlim = range(xr), pch = 16)
segments(xr[ii,1], 1:nr, xr[ii,2], 1:nr)
mtext(side = 3, cex = 1.5, outer = TRUE, line = -1.5,
      "E coli Core Metabolism, optimal solution and ranges")
}
}

\references{
originated from the urlhttp://gcrg.ucsd.edu/Downloads/Flux_Balance_Analysis

Edwards,J.S., Covert, M., and Palsson, B., (2002)
Metabolic Modeling of Microbes: the Flux Balance Approach,
Environmental Microbiology, 4(3): pp. 133-140.
}
\keyword{datasets}

