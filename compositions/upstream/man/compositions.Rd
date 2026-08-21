\name{compositions-package}
\alias{compositions-package}
\alias{compositions}
\docType{package}
\title{
\packageTitle{compositions}
}
\description{
"compositions" is a package for the analysis of compositional and
multivariate positive data (generally called "amounts"), based on several alternative approaches.
}
\details{
  The DESCRIPTION file:
\packageDESCRIPTION{compositions}
\packageIndices{compositions}
   To get detailed "getting started" introduction use
   \code{help.start()} or \code{help.start(browser="myfavouritebrowser")}
   Go to "Packages" then "compositions" and then "overview"
   and then launch the file "UsingCompositions.pdf" from there. Please
   also check the web-site: http://www.stat.boogaart.de/compositions/ for
   improved material and our new book expected to appear spring 2009.  
   \cr
   The package is devoted to the analysis of multiple amounts. Amounts
   have typically non-negative values, and often sum up to 100\% or one. These
   constraints lead to spurious effects on the covariance structure, 
   as pointed out by Chayes (1960). The problem is treated rigorously 
   in the monography by Aitchison (1986), 
   who characterizes compositions as vectors having a relative scale,
   and identifies its sample space with the D-part simplex. 
   However still (i.e. 2005) most statistical packages do not
   provided any support for this scale.
   \cr
   The grounding idea of the package exploits the class concept: 
   the analyst gives the data a compositional or amount class, and
   then all further analysis are (should be) automatically  done
   in a consistent way, e.g. \code{x <- acomp(X); plot(x)} 
   should plot the data as a composition (in a ternary diagram) 
   directly without any further interaction of the user. 
   \cr
   The package provides four different approaches to analyse
   amounts. These approaches are associated to four R-classes,
   representing four different geometries of the sampling space of
   amounts. These geometries depend on two questions: whether the total sum
   of the amounts is a relevant information, and which is the meaningful
   measure of difference of the data.
   \cr
   
   \code{\link{rplus}} : (Real Plus) The total amount matters, and amounts should be
   compared on an absolute basis. i.e. the difference between 1g and
   2g is the same as the difference between 1kg and 1001g, one gram.
   \cr
   \code{\link{aplus}} : (Aitchison Plus) The total amount matters,
   but amounts should be compared relatively, i.e. the difference
   between 1mg and 2mg is the same as that of 1g and 2g: the double.
   \cr
   \code{\link{acomp}} : (Aitchison composition) the total amount is constant 
   (or an artifact of the sampling/measurement procedure), and the meaningful 
   difference is a relative one. This class follows
   the original proposals of Aitchison.
   \cr
   \code{\link{rcomp}} : (Real composition) the sum
   is a constant, and the difference in amount from 0\% to 1\% and from
   10\% to 11\% is regarded as equal. This class represents the
   raw/naive treatment of compositions as elements of the real simplex based
   on an absolute geometry. This treatment is implicitly used 
   in most amalgamation problems. However the whole approach suffers
   from the drawbacks and problems discussed in Chayes (1960) and Aitchison
   (1986).
   \cr
   The aim of the package is to provide all the functionality to do a
   consistent analysis in all of these approaches and to make the
   results obtained with different geometries as easy to compare as possible. 
 }
 \author{
\packageAuthor{compositions}

Maintainer: \packageMaintainer{compositions}
}
 \references{

Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

  Billheimer, D., P. Guttorp, W.F. and Fagan (2001) Statistical interpretation of species composition,
  \emph{Journal of the American Statistical Association}, \bold{96} (456), 1205-1214\cr

  Chayes, F. (1960). On correlation between variables of constant
  sum. Journal of Geophysical Research 65~(12), 4185--4193.

  Pawlowsky-Glahn, V. and J.J. Egozcue (2001) Geometric approach to
  statistical analysis on the simplex. \emph{SERRA} \bold{15}(5), 384-398\cr

  Pawlowsky-Glahn, V. (2003) Statistical modelling on coordinates. In: 
  Thi\'o -Henestrosa, S. and Mart\'in-Fern\'a ndez, J.A. (Eds.)
  \emph{Proceedings of the 1st International Workshop on Compositional Data Analysis},
  Universitat de Girona, ISBN 84-8458-111-X, \url{https://ima.udg.edu/Activitats/CoDaWork03/}\cr

  Mateu-Figueras, G. and Barcel\'o-Vidal, C. (Eds.)
  \emph{Proceedings of the 2nd International Workshop on Compositional Data Analysis},
  Universitat de Girona, ISBN 84-8458-222-1, \url{https://ima.udg.edu/Activitats/CoDaWork05/}\cr

    van den Boogaart, K.G. and R. Tolosana-Delgado (2008) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers &
    Geosciences}, 34 (4), pages 320-338, \doi{10.1016/j.cageo.2006.11.017}.
 }
\examples{
library(compositions)      # load library
data(SimulatedAmounts)     # load data sa.lognormals
x <- acomp(sa.lognormals)  # Declare the dataset to be compositional
                           # and use relative geometry
plot(x)                    # plot.acomp : ternary diagram
ellipses(mean(x),var(x),r=2,col="red")  # Simplex 2sigma predictive region
pr <- princomp(x)
straight(mean(x),pr$Loadings) 

x <- rcomp(sa.lognormals)  # Declare the dataset to be compositional
                           # and use absolute geometry
plot(x)                    # plot.acomp : ternary diagram
ellipses(mean(x),var(x),r=2,col="red")  # Real 2sigma predictive region
pr <- princomp(x)          
straight(mean(x),pr$Loadings) 
}
\seealso{
  \link{compositions-package}, \link{missingsInCompositions},
  \link{robustnessInCompositions}, \link{outliersInCompositions},
}
\note{The package compositions has grown a lot in the last year:
  missings, robust estimations, outlier detection and classification,
  codadendrogram. This makes everything much more complex especially
  from the side of programm testing. Thus we would like to urge our
  users to report all errors and problems of the lastest version (please
  check first) to
  support@boogaart.de.

}
\keyword{package}
