\name{fish}
\alias{fish}
\docType{data}
\title{ Fish Catch Data Set}
\description{
  The Fish Catch data set contains measurements on 159 fish caught in the lake
    Laengelmavesi, Finland.
}
\usage{data(fish)}

\format{
  A data frame with 159 observations on the following 7 variables.

\describe{
    \item{\code{Weight}}{Weight of the fish (in grams)}
    \item{\code{Length1}}{Length from the nose to the beginning of the tail (in cm)}
    \item{\code{Length2}}{Length from the nose to the notch of the tail (in cm)}
    \item{\code{Length3}}{Length from the nose to the end of the tail (in cm)}
    \item{\code{Height}}{Maximal height as \% of Length3}
    \item{\code{Width}}{Maximal width as \% of Length3}
    \item{\code{Species}}{Species}
  }

}

\details{
    The Fish Catch data set contains measurements on 159 fish caught in the lake
    Laengelmavesi, Finland. For the 159 fishes of 7 species the weight,
    length, height, and width were measured. Three different length
    measurements are recorded: from the nose of the fish to the
    beginning of its tail, from the nose to the notch of its tail and
    from the nose to the end of its tail. The height and width are
    calculated as percentages of the third length variable. This results
    in 6 observed variables, Weight, Length1, Length2, Length3, Height, Width.
    Observation 14 has a missing value in variable Weight, therefore
    this observation is usually excluded from the analysis. The last
    variable, Species, represents the grouping structure: the 7
    species are 1=Bream, 2=Whitewish, 3=Roach, 4=Parkki, 5=Smelt,
    6=Pike, 7=Perch. This data set was also analyzed in the context 
    of robust Linear Discriminant Analysis by Todorov (2007), 
    Todorov and Pires (2007).
}
\source{
    Journal of Statistical Education, Fish Catch Data Set,
    [https://jse.amstat.org/datasets/fishcatch.dat.txt] accessed
    November, 2023.
}
\references{
Todorov, V. (2007 Robust selection of variables in linear
  discriminant analysis, \emph{Statistical Methods and Applications}, \bold{15},
  395--407, doi:10.1007/s10260-006-0032-6.

Todorov, V. and Pires, A.M. (2007) Comparative performance of several robust
  linear discriminant analysis methods, \emph{REVSTAT Statistical Journal}, \bold{5},
  63--83.
}
\examples{
    data(fish)

    # remove observation #14 containing missing value
    fish <- fish[-14,]

    # The height and width are calculated as percentages 
    #   of the third length variable
    fish[,5] <- fish[,5]*fish[,4]/100
    fish[,6] <- fish[,6]*fish[,4]/100
 
    # plot a matrix of scatterplots
    pairs(fish[1:6],
          main="Fish Catch Data",
          pch=21,
          bg=c("red", "green3", "blue", "yellow", "magenta", "violet", 
          "turquoise")[unclass(fish$Species)])

}
\keyword{datasets}
