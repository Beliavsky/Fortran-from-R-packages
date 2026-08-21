\name{Aar}
\docType{data}
\alias{Aar}
\title{Composition of glaciar sediments from the Aar massif (Switzerland)}
\description{
 Geochemical composition of glaciar sediments from the Aar massif region (Switzerland), major oxides and trace elements. 
 }
\usage{
data(Aar)
}
\source{
von Eynatten H.; Tolosana-Delgado, R.; Karius, V (2012) Sediment generation in modern glacial settings: Grain-size and source-rock control on sediment composition. Sedimentary Geology 280 (1): 80-92 \doi{10.1016/j.sedgeo.2012.03.008}
}

\details{
 Composition of recent sediments of several morraines and streams from glaciers around the Aar massif, including both major oxides and trace elements. The major oxides are expressed in weight percent (total sum reported in column \code{SumOxides}), from Silica (\code{SiO2}, column 3) to total Iron 3 Oxide (\code{Fe2O3t}, column12, incorporating FeO recasted to Fe2O3). The trace elements are reported in parts per million (ppm, mg/Kg) between columns 14 (\code{Ba}) and 29 (\code{Nd}). Partial sum of the trace elements (in ppm) and of all traces and major oxides (in \%) are also reported.
 
 Apart of the compositional information, two covariables are included: Sample and GS. The variable \code{Sample} reports the ID of the sample material. This material was sieved in 11 grain size fractions, and each fraction was analysed separately after drying. The grain size fraction of each subsample is reported in variable \code{GS}, representing the upper limit of the size fraction reported in \eqn{\phi}{phi} scale, e.g. the binary log transformation of the average diameter \eqn{\bar{d}}{dmean}
 \deqn{ \phi = -\log_2 (\bar{d})}{phi=log_2(dmean)}

 The Aar is a granitic-granodioritic-gneissic massif of the Alps, in Switzerland, comprised of several intrusions with different compositions within the range of granitoid lithologies. Details of the region, mineralogy, procedures and study questions behind the data can be found in von Eynatten at al (2012) and references thereon.
}

\references{
   von Eynatten H.; Tolosana-Delgado, R.; Karius, V (2012) Sediment generation in modern glacial settings: Grain-size and source-rock control on sediment composition. Sedimentary Geology 280 (1): 80-92 \doi{10.1016/j.sedgeo.2012.03.008}
 }
\note{Courtesy of H. von Eynatten}
\keyword{datasets}
