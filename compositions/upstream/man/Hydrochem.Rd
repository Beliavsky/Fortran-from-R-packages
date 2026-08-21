\name{Hydrochem}
\docType{data}
\alias{Hydrochem}
\title{Hydrochemical composition data set of Llobregat river basin water (NE Spain)}
\description{
  Contains a hydrochemical amount/compositional data set obtained from
  several rivers in the Llobregat river basin, in northeastern Spain.
}
\usage{
data(Hydrochem)
}
\format{Data matrix with 485 cases and 19 variables.}
\source{
  The dataset is also accessible in Otero et al. (2005), and are here
  included under the  GNU Public Library Licence Version 2 or newer.
}
\author{Raimon Tolosana-Delgado}
\details{
  This hydrochemical data set contains measurements of 14 components, H,
  Na, K, Ca, Mg, Sr, Ba, \eqn{\mathrm{NH}_4}{NH4},
  Cl, \eqn{\mathrm{HCO}_3}{HCO3}, \eqn{\mathrm{NO}_3}{NO3},
  \eqn{\mathrm{SO}_4}{SO4}, \eqn{\mathrm{PO}_4}{PO4},
  TOC. From them, hydrogen was derived by inverting the relationship
  between its stable form in water,
  \eqn{\mathrm{H}_3\mathrm{O}^+}{H3O+}, and pH. Details can be
  found in Otero et al. (2005). Each of these parameters is measured
  approximately once each month during 2 years in 31 stations, placed
  along the rivers and main tributaries of the Llobregat river, one of the
  medium rivers in northeastern Spain.

  The Llobregat river drains an area of 4948.2 \eqn{km^2}{km2}, and it
  is 156.6 km long, with two main tributaries,
  Cardener and Anoia. The headwaters of Llobregat and Cardener are in a
  rather unpolluted area of the
  Eastern Pyrenees. Mid-waters these rivers flow through a densely
  populated and industrialized area, where potash mining activity occurs
  and there are large salt mine tailings stored with no water
  proofing. There, the main land use is agriculture and stockbreeding. The
  lower course flows through one of the most densely populated areas of
  the Mediterranean region (around the city of Barcelona) and the river
  receives large inputs from industry and urban origin, while intensive
  agriculture activity is again present in the Llobregat delta. Anoia is
  quite different. Its headwaters are in an agricultural area, downwaters
  it flows through an industrialized zone (paper mills, tannery and
  textile industries), and near the confluence with Llobregat the main
  land use is agriculture again, mainly vineyards, with a decrease in
  industry and urban contribution. Given this variety in geological
  background and human activities, the sample has been splitted in four
  groups (higher Llobregat course, Cardener, Anoia and lower Llobregat
  course), which in turn are splitted into main river and tributaries
  (Otero et al, 2005). Information on these groupings, the sampling
  locations and sampling time is included in 5 complementary variables.
}

\references{
  Otero, N.; R. Tolosana-Delgado, A. Soler, V. Pawlowsky-Glahn
  and A. Canals (2005).  Relative vs. absolute statistical analysis of
  compositions: A comparative study of surface waters of a Mediterranean
  river. Water Research, 39(7): 1404-1414. \doi{10.1016/j.watres.2005.01.012}.

  Tolosana-Delgado, R.; Otero, N.; Pawlowsky-Glahn, V.; Soler, A. (2005). 
  Latent Compositional Factors in The Llobregat River Basin (Spain) 
  Hydrogeochemistry. Mathematical Geology 37(7): 681-702.
}
\examples{
data(Hydrochem)
cHydrochem=Hydrochem[, 6:19]
biplot(princomp(rplus(cHydrochem)))
biplot(princomp(rcomp(cHydrochem)))

biplot(princomp(aplus(cHydrochem)))
biplot(princomp(acomp(cHydrochem)))
}
\keyword{datasets}
