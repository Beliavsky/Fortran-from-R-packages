\name{jura}
\docType{data}
\alias{juraset}
\alias{jura259}
\title{The jura dataset}
\description{
  A geochemical dataset from the Swiss Jura.
}
\usage{
data(juraset)
data(jura259)
}
\format{A 359x11 or 259x11 dataframe }
\source{
  AI-Geostats
}
\details{

The JURA data set
provided by J.-P. Dubois, IATE-Paedologie, Ecole Polytechnique Federale de Lausanne,
1015 Lausanne, Switzerland.
Spatial coordinates and values of categorial and continuous attributes
at the 359 sampled sites. The 100 test locartions are denoted with a star.
Rock Types: 1: Argovian, 2: Kimmeridgian, 3: Sequanian, 4: Portlandian, 5: Quaternary.
Land uses: 1: Forest, 2: Pasture, 3: Meadow , 4: Tillage 

\tabular{rll}{
 \tab X\tab X location coordinate\cr
 \tab Y\tab Y location coordinate\cr
 \tab Rock\tab  Categorical: rocktype,\cr
 \tab Land\tab  Categorical: land usage\cr
 \tab Cd\tab  element amount,\cr
 \tab Cu\tab  element amount,\cr
 \tab Pb\tab  element amount,\cr
 \tab Co\tab  element amount,\cr
 \tab Cr\tab  element amount,\cr
 \tab Ni\tab  element amount,\cr
}
All 3-part compositions sum to one.
}
\examples{
\dontrun{
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
lrv <- logratioVariogram(comp,X,maxdist=1,nbins=10)
plot(lrv)
}

}
\references{
  
Atteia, O., Dubois, J.-P., Webster, R., 1994,
Geostatistical analysis of soil contamination in the Swiss Jura:
Environmental Pollution 86, 315-327

Webster, R., Atteia, O., Dubois, J.-P., 1994,
Coregionalization of trace metals in the soil in the Swiss Jura:
European Journal of Soil Science 45, 205-218
 }
\keyword{datasets}
