\name{OsloTransect}
\alias{OsloTransect}
\docType{data}
\title{ Oslo Transect Data }
\description{
The oslo Transect data set contains 360 samples of different 
plant species collected along a 120 km transect running through the city of Oslo,
Norway.
}
\usage{data(OsloTransect)}
\format{
  A data frame with 360 observations on the following 38 variables.
  \describe{
    \item{\code{X.ID}}{a numeric vector, unique ID of the sample}
    \item{\code{X.MAT}}{a factor with levels \code{BBA} \code{BIL} \code{BWO} \code{FER} \code{MOS} \code{ROG} \code{SNE} \code{STW} \code{TWI}}
    \item{\code{XCOO}}{a numeric vector, X coordinate}
    \item{\code{YCOO}}{a numeric vector, Y coordinate}
    \item{\code{XCOO_km}}{a numeric vector}
    \item{\code{YCOO_km}}{a numeric vector}
    \item{\code{X.FOREST}}{a factor with levels \code{BIRSPR} \code{MIXDEC} \code{PINE} \code{SPRBIR} \code{SPRPIN} \code{SPRUCE}}
    \item{\code{DAY}}{a numeric vector}
    \item{\code{X.WEATHER}}{a factor with levels \code{CLOUD} \code{MOIST} \code{NICE} \code{RAIN}}
    \item{\code{ALT}}{a numeric vector}
    \item{\code{X.ASP}}{a factor with levels \code{} \code{E} \code{FLAT} \code{N} \code{NE} \code{NW} \code{S} \code{SE} \code{SW} \code{W}}
    \item{\code{X.GRVEG}}{a factor with levels \code{BLGR} \code{BLLY} \code{BLMOLI} \code{BLUE} \code{BLUGRA} \code{GRAS} \code{GRBLU} \code{GRFE} \code{GRMO} \code{LYLI} \code{MIX} \code{MOGR} \code{MOSS}}
    \item{\code{X.FLITHO}}{a factor with levels \code{CAMSED} \code{GNEID_O} \code{GNEIS_O} \code{GNEIS_R} \code{MAGM} \code{MICSH}}
    \item{\code{Ag_ppb}}{a numeric vector}
    \item{\code{As_ash}}{a numeric vector}
    \item{\code{B}}{a numeric vector}
    \item{\code{Ba}}{a numeric vector}
    \item{\code{Ca}}{a numeric vector}
    \item{\code{Cd}}{a numeric vector}
    \item{\code{Co}}{a numeric vector}
    \item{\code{Cr}}{a numeric vector}
    \item{\code{Cu}}{a numeric vector}
    \item{\code{Fe}}{a numeric vector}
    \item{\code{Hg_ppb}}{a numeric vector}
    \item{\code{K}}{a numeric vector}
    \item{\code{La}}{a numeric vector}
    \item{\code{LOI}}{a numeric vector}
    \item{\code{Mg}}{a numeric vector}
    \item{\code{Mn}}{a numeric vector}
    \item{\code{Mo}}{a numeric vector}
    \item{\code{Ni}}{a numeric vector}
    \item{\code{P}}{a numeric vector}
    \item{\code{Pb}}{a numeric vector}
    \item{\code{S}}{a numeric vector}
    \item{\code{Sb}}{a numeric vector}
    \item{\code{Sr}}{a numeric vector}
    \item{\code{Ti}}{a numeric vector}
    \item{\code{Zn}}{a numeric vector}
  }
}
\details{
Samples of different plant species were collected along a 120 km 
transect running through the city of Oslo,
Norway (forty samples each of leaves, needles,roots or 
barks of several plant species), and the concentrations of 
25 chemical elements for the sample materials are reported. 
The factors that influenced the observed element 
concentrations in the sample materials were investigated.
This data set was used in Todorov and Filzmoser (2007) for 
illustration of the robust statistics for one-way MANOVA implemented in the function 
\code{\link{Wilks.test}}.
}
\source{
REIMANN,C., ARNOLDUSSEN,A., BOYD,R., FINNE,T.E., 
NORDGULEN,Oe., VOLDEN,T. and  ENGLMAIER,P. (2006) 
The Influence of a city on element contents of a terrestrial moss 
(Hylocomium splendens), \emph{The Science of the Total Environment} \bold{369} 419--432.

REIMANN,C., ARNOLDUSSEN,A., BOYD,R., FINNE,T.E., KOLLER,F., NORDGULEN,Oe., and ENGLMAIER,P. (2007) 
Element contents in leaves of four plant species (birch, mountain ash, fern and spruce) 
along anthropogenic and geogenic concentration gradients, 
\emph{The Science of the Total Environment} \bold{377} 416--433.

REIMANN,C., ARNOLDUSSEN,A., FINNE,T.E., KOLLER,F., NORDGULEN,Oe.,  and ENGLMAIER,P., (2007) 
Element contents in birch leaves, bark and wood under different 
anthropogenic and geogenic conditions,
\emph{Applied Geochemistry}, \bold{22} 1549--1566.
}
\references{
Todorov V. and Filzmoser P. (2007) Robust statistic for the one-way MANOVA,
  \emph{submetted to the Journal of Environmetrics}.
}
\examples{
data(OsloTransect)
str(OsloTransect)

##
##  Log-transform the numerical part of the data, 
##  choose the desired groups and variables and 
##  perform the classical Wilks' Lambda test
##
OsloTransect[,14:38] <- log(OsloTransect[,14:38])
grp <- OsloTransect$X.FLITHO
ind <- which(grp =="CAMSED" | grp == "GNEIS_O" |
    grp == "GNEIS_R" | grp=="MAGM")
(cwl <- Wilks.test(X.FLITHO~K+P+Zn+Cu,data=OsloTransect[ind,]))

##
## Perform now the robust MCD based Wilks' Lambda test. 
##  Use the already computed multiplication factor 'xd' and 
##  degrees of freedom 'xq' for the approximate distribution.
##

xd <- -0.003708238
xq <- 11.79073
(mcdwl <- Wilks.test(X.FLITHO~K+P+Zn+Cu,data=OsloTransect[ind,], 
    method="mcd", xd=xd, xq=xq))
}
\keyword{datasets}
