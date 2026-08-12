\name{Kinship82_Consensus}
\alias{Kinship82_Consensus}
\title{Gordon-Vichi Kinship82 Consensus Partition Data}
\description{
  The soft (\dQuote{fuzzy}) consensus partitions for the Rosenberg-Kim
  kinship terms partition data given in \bibcitet{Gordon+Vichi:2001}.
}
\usage{data("Kinship82_Consensus")}
\format{
  A named cluster ensemble of three soft partitions of the 15 kinship
  terms into three classes.
}
\details{
  The elements of the ensemble are named \code{"MF1"}, \code{"MF2"}, and
  \code{"JMF"}, and correspond to the consensus partitions obtained by
  applying models 1, 2, and 3 in \bibcitet{Gordon+Vichi:2001} to the
  kinship   terms partition data in \bibcitet{Rosenberg:1982}, which are
  available as data set \code{\link{Kinship82}}.
}
\source{
  Table 6 in \bibcitet{Gordon+Vichi:2001}.
}
\references{
  \bibshow{Gordon+Vichi:2001, Rosenberg:1982}
}
\examples{
## Load the consensus partitions.
data("Kinship82_Consensus")
## Fuzziness using the Partition Coefficient.
cl_fuzziness(Kinship82_Consensus)
## (Corresponds to 1 - F in the source.)
## Dissimilarities:
cl_dissimilarity(Kinship82_Consensus)
cl_dissimilarity(Kinship82_Consensus, method = "comem")
}
\keyword{datasets}
