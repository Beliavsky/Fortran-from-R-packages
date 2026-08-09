\name{Chemtax}
\docType{data}
\alias{Chemtax}
\title{
  An overdetermined linear inverse problem: estimating algal composition
  based on pigment biomarkers.
}
\description{
  Input files for assessing the algal composition of a field sample,
  based on the pigment composition of the algal groups (measured in the
  laboratory) and the pigment composition of the field sample.

  In the example there are 8 types of \bold{algae}:
  \itemize{
    \item Prasinophytes
    \item Dinoflagellates
    \item Cryptophytes
    \item Haptophytes type 3 (Hapto3s)
    \item Haptophytes type 4 (Hapto4s)
    \item Chlorophytes
    \item Cynechococcus
    \item Diatoms
  }

  and 12 \bold{pigments}:

  \code{Perid} = Peridinin,
  \code{19But} = 19-butanoyloxyfucoxanthin,
  \code{Fucox} = fucoxanthin,

  \code{19Hex} = 19-hexanoyloxyfucoxanthin,
  \code{Neo} = neoxanthin,
  \code{Pras} = prasinoxanthin,

  \code{Viol} = violaxanthin,
  \code{Allox} = alloxanthin,
  \code{Lutein} = lutein,

  \code{Zeax} = zeaxanthin,
  \code{Chlb} = chlorophyll b,
  \code{Chla} = chlorophyll a


  The input data consist of:
  \enumerate{
    \item the pigment composition of the various algal groups, for instance
      determined from cultures (the Southern Ocean example -table 4- in
      Mackey et al., 1996)
    \item the pigment composition of a field sample.
  }
  Based on these data, the algal composition of the field sample is estimated,
  under the assumption that the pigment composition of the field sample is a
  weighted avertage of the pigment composition of algae present in the sample,
  where weighting is proportional to their biomass.

  As there are more measurements (12 pigments) than unknowns (8 algae),
  the resulting linear system is overdetermined.

  It is thus solved in a least squares sense (using function \code{lsei}):

  \deqn{\min(||Ax-b||^2)} subject to
  \deqn{Ex=f}
  \deqn{Gx>=h}

  If there are 2 algae \code{A,B}, and 3 pigments \code{1,2,3} then the 3
  approximate equalities (\eqn{A*x=B}) would be:
  \deqn{f_{1,S} = p_{A}*f_{A,1}+p_{B}*f_{B,1}}
  \deqn{f_{2,S} = p_{A}*f_{A,2}+p_{B}*f_{B,3}}
  \deqn{f_{3,S} = p_{A}*f_{A,3}+p_{B}*f_{B,3}}

  where \eqn{p_{A}} and \eqn{p_{b}} are the (unknown) proportions of algae
  A and B in the field sample (S), and
  \eqn{f_{A,1}} is the relative amount of pigment 1 in alga A, etc...

  The equality ensures that the sum of fractions equals 1:
  \deqn{1 = p_{A}+p_{b}}


  and the inequalities ensure that fractions are positive numbers
  \deqn{p_{A} > 0}
  \deqn{p_{B} > 0}
  
  It should be noted that in the actual \emph{Chemtax} programme the problem
  is solved in a more complex way. In Chemtoz, the A-coefficients are also
  allowed to vary, while here they are taken as they are (constant).
  Chemtax then finds the "best fit" by fitting both the fractions, and the
  non-zero coefficients in A.
}

\usage{Chemtax}

\format{
  A list with the input ratio matrix (\code{Ratio}) and a vector with the
  field data (\code{Field})
  \itemize{
    \item The input ratio matrix \code{Ratio} contains the pigment compositions
      (columns) for each algal group (rows);
      the compositions are scaled relative to \code{Chla} (last column).
    \item The vector with the \code{Field} data contains the pigment composition
      of a sample in the field, also scaled relative to \code{Chla};
      the pigments are similarly ordened as for the input ratio matrix.
  }
  The rownames of matrix \code{Ratio} are the algal group names, columnames
  of \code{Ratio} (=names of \code{Field}) are the pigments
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>.
}

\examples{
# 1. Graphical representation of the chemtax example input data
palette(rainbow(12, s = 0.6, v = 0.75))

mp     <- apply(Chemtax$Ratio, MARGIN = 2, max)
pstars <- rbind(t(t(Chemtax$Ratio)/mp) ,
                  sample = Chemtax$Field/max(Chemtax$Field))
stars(pstars, len = 0.9, key.loc = c(7.2, 1.7),scale=FALSE,ncol=4,
      main = "CHEMTAX pigment composition", draw.segments = TRUE,
      flip.labels=FALSE)

# 2. Estimating the algal composition of the field sample
Nx     <-nrow(Chemtax$Ratio)

# equations that have to be met exactly Ex=f: 
# sum of all fraction must be equal to 1.
EE <- rep(1, Nx)
FF <- 1

# inequalities, Gx>=h:
# all fractions must be positive numbers
GG <- diag(nrow = Nx)
HH <- rep(0, Nx)

# equations that must be reproduced as close as possible, Ax ~ b
# = the field data; the input ratio matrix and field data are rescaled
AA     <- Chemtax$Ratio/rowSums(Chemtax$Ratio)
BB     <- Chemtax$Field/sum(Chemtax$Field)

# 1. Solve with lsei method
X <- lsei(t(AA), BB, EE, FF, GG, HH)$X
(Sample <- data.frame(Algae = rownames(Chemtax$Ratio),
                      fraction = X))

# plot results
barplot(X, names = rownames(Chemtax$Ratio), col = heat.colors(8),
        cex.names = 0.8, main = "Chemtax example solved with lsei")

# 2. Bayesian sampling; 
# The standard deviation on the field data is assumed to be 0.01
# jump length not too large or NO solutions are found!
xs <- xsample(t(AA), BB, EE, FF, GG, HH, sdB = 0.01, jmp = 0.025)$X
pairs(xs, main= "Chemtax, Bayesian sample")
}
\references{
  Mackey MD, Mackey DJ, Higgins HW, Wright SW, 1996.
  CHEMTAX - A program for estimating class abundances from chemical markers:
  Application to HPLC measurements of phytoplankton. Marine Ecology-Progress
  Series 144 (1-3): 265-283.

  Van den Meersche, K., Soetaert, K., Middelburg, J., 2008.
  A Bayesian compositional estimator for microbial taxonomy based on biomarkers.
  Limnology and Oceanography Methods, 6, 190-199.

  R-package \code{BCE}
}
\seealso{
  \code{\link{lsei}}, the function to solve for the algal composition of the
  field sample.
}

\keyword{datasets}

