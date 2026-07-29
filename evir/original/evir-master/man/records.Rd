\name{records}
\alias{records}
\title{Calculate Record Development}
\description{
  Creates a data frame showing the development of records in 
  a dataset and calculating the expected behaviour for iid data.
}
\usage{
records(data, do.plot = TRUE, conf.level = 0.95, \dots)
}
\arguments{
  \item{data}{data vector}
  \item{do.plot}{whether a plot of record development should be
    created}
  \item{conf.level}{confidence level for record development plot}
  \item{\dots}{graphics parameters}
}
\value{
 A data frame. 
}
\details{
  Records are counted and the observations at which they occur 
  recorded. This is compared with the expected behaviour for iid 
  data. 
}

\examples{
\dontrun{data(danish)}
\dontrun{records(danish)} 
# Record fire insurance losses in Denmark 
}
\keyword{hplot}
