\name{TsayFiles}
\docType{data}
\alias{FinTS.url}
\alias{TsayFiles}
\title{
  List of the names of files downloaded from the "Analysis of Financial
  Data" web site.
}
\description{
  A list organized by chapter and text vs. exercises of the files
  downloaded from the web site associated with Tsay (2005) Analysis of
  Financial Time Series, 2nd ed. (Springer) and stored in
  "~library/FinTS/scripts/TsayFiles".  These facilitate the process of
  creating and updating the 'FinTS' package (and documenting the
  creation process).
}
\usage{
TsayFiles
}
\format{
  \describe{
    \item{TsayFiles}{
      A list with names 'ch01', 'ch02', ..., 'ch12' for components
      describing the files associated with the corresponding chapter.

      Each chapter component is a list with 'text' and 'exercises'
      components, where 'text' and 'exercises' are each a character
      array giving names for 'data', 'file', 'url', and 'found' for the
      data referenced in the text or exercises of that chapter:

      \describe{
	\item{data}{
	  'file' without the extension, e.g, 'd-ibmvwewsp6203' for daily
	  simple returns of IBM, VW, EW, SP (7/3/62-12/31/03) = 'file'
	  without the extension.
	}
	\item{file}{
	  short file name = 'data' plus the extension, e.g.,
	  'd-ibmvwewsp6203'.txt' for daily simple returns of IBM, VW,
	  EW, SP (7/3/62-12/31/03)
	}
	\item{url}{
	  universal resource locator for the data, e.g.,

	  "http://faculty.chicagogsb.edu/ruey.tsay/teaching/fts2/d-ibmvwewsp6203.txt"
	}
	\item{found}{
	  'TRUE' if the data were found, 'FALSE' if the attempt to
	  access the url failed.
	}
      }

      NOTES:

      (1) 13 files are referenced twice, and 2 are referenced three
      times on the web page.  This redundancy is retained in
      'TsayFiles'.

      (2) A few files (most noticably some with with '.dat' extension)
      are referenced in the HTML code without an apparent visible link.
      These 'invisible files' are retained in 'TsayFiles'.
    }
    \item{FinTS.url}{
      A character string giving the universal resource locator (URL)
      associated with the Tsay (2005) book:

      % was: "http://faculty.chicagogsb.edu/ruey.tsay/teaching/fts2"
      \code{FinTS.url <- "https://faculty.chicagobooth.edu/ruey-s-tsay/teaching"}
    }
  }
}
\source{
  Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Springer)
}

%\seealso{\code{\link{eval.penalty}}}
\examples{
data(TsayFiles)
TsayFiles$ch01$exercises
}
\keyword{datasets}
