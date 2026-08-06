%  file sn/man/wines.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{wines}
\alias{wines}
\docType{data}
\encoding{UTF-8}
\title{Piedmont wines data}

\description{Data refer to chemical properties of 178 specimens of three 
  types of wine produced in the Piedmont region of Italy. }
   
\usage{data(wines)}
 
\format{
  A data frame with 178 observations on the following 28 variables.
  \tabular{ll}{
    \code{wine}\tab wine name (categorical,  
         levels: \code{Barbera}, \code{Barolo}, \code{Grignolino})\cr
    \code{alcohol}\tab alcohol percentage (numeric)\cr
    \code{sugar}\tab sugar-free extract (numeric)\cr
    \code{acidity}\tab fixed acidity (numeric)\cr
    \code{tartaric}\tab tartaric acid (numeric)\cr 
    \code{malic}\tab malic acid (numeric)\cr
    \code{uronic}\tab uronic acids (numeric)\cr
    \code{pH}\tab pH (numeric)\cr
    \code{ash}\tab ash (numeric)\cr
    \code{alcal_ash}\tab alcalinity of ash (numeric)\cr
    \code{potassium}\tab potassium (numeric)\cr
    \code{calcium}\tab calcium (numeric)\cr
    \code{magnesium}\tab magnesium (numeric)\cr
    \code{phosphate}\tab phosphate (numeric)\cr
    \code{cloride}\tab chloride (numeric)\cr
    \code{phenols}\tab total phenols (numeric)\cr
    \code{flavanoids}\tab flavanoids (numeric)\cr
    \code{nonflavanoids}\tab nonflavanoid phenols (numeric)\cr
    \code{proanthocyanins}\tab proanthocyanins (numeric)\cr
    \code{colour}\tab colour intensity (numeric)\cr
    \code{hue}\tab hue (numeric)\cr
    \code{OD_dw}\tab \eqn{OD_{280}/OD_{315}}{OD₂₈₀/OD₃₁₅} of diluted wines
          (numeric)\cr
    \code{OD_fl}\tab \eqn{OD_{280}/OD_{315}}{OD₂₈₀/OD₃₁₅} of flavanoids 
          (numeric)\cr          
    \code{glycerol}\tab glycerol (numeric)\cr
    \code{butanediol}\tab 2,3-butanediol (numeric)\cr
    \code{nitrogen}\tab total nitrogen (numeric)\cr
    \code{proline}\tab proline (numeric)\cr
    \code{methanol}\tab methanol (numeric)\cr
  }
}

\details{
The data represent 27 chemical measurements on each of 178 wine specimens 
belonging to three types of wine produced in the Piedmont region of Italy. 
The data have been presented and examined by Forina \emph{et al.} (1986) and 
were freely accessible from the \acronym{PARVUS} web-site until it was active.
These data or, more often, a subset of them are now available from various
places, including some \R packages. The present dataset includes 
all variables available on the \acronym{PARVUS} repository, which are
the variables listed by Forina \emph{et al.} (1986) with the exception 
of \sQuote{Sulphate}. Moreover, it reveals the undocumented fact that
the original dataset appears to include also the vintage year;
see the final portion of the \sQuote{Examples} below.}

\source{
Forina, M., Lanteri, S. Armanino, C., Casolino, C., Casale, M. and Oliveri, P. 
\acronym{V-PARVUS 2008}: an extendible package of programs for esplorative 
data analysis, classification and regression analysis. 
Dip. Chimica e Tecnologie Farmaceutiche ed Alimentari, 
Università di Genova, Italia. Web-site (not accessible as of 2014): 
\samp{http://www.parvus.unige.it}
}

\references{ 
 Forina M., Armanino C., Castino M. and Ubigli M. (1986).  
 Multivariate data analysis as a discriminating method of the origin of wines. 
 \emph{Vitis} \bold{25}, 189--201.
}

\examples{
data(wines)
pairs(wines[,c(2,3,16:18)], col=as.numeric(wines$wine))
#
code <- substr(rownames(wines), 1, 3)
table(wines$wine, code)
#
year <- as.numeric(substr(rownames(wines), 6, 7))
table(wines$wine, year)
# coincides with Table 1(a) of Forina et al. (1986)
}

\keyword{datasets}
