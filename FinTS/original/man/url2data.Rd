\name{url2data}
\alias{url2data}
\title{ Create local copies of files read from urls }
\description{
  Call 'download.file' with each element of a vector of character
  strings assumed to be URLs, create a local copy for each, and return a
  character matrix summarizing what was done.  
}
\usage{
url2data(url.) 
}
\arguments{
  \item{url.}{
    a vector of character strings assumed to be URLs, whose names are
    assumed to be the names to be used for local copies of the URLs.  
  }
}
\details{
  1.  \code{fili <- names(urls.)}

  2.  \code{dati <- fili} without its extension, i.e., the part
  following the last '.'

  3.  \code{for(i in 1:length(url.)try(download.file(url.[i], fili[i]));}
  \code{found[i] = TRUE} if something found and FALSE if not.  

  4.  Return a character matrix with 4 columns:  data = dati, file =
  fili, url = url., and found.  
}
\examples{
# See ~R\library\FinTS\scripts\TsayFiles.R  
}
\keyword{file}
