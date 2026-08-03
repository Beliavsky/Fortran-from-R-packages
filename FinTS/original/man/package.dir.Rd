\name{package.dir}
\alias{package.dir}
\title{ Directory of a package  }
\description{
  Display a partial or complete directory of a package.  By default,
  suppress common package contents to focus on 'demo', 'doc', 'scripts',
  and similar subdirectories whose contents might contain examples that
  could make it easier to learn capabilities of the package.
}
\usage{
package.dir(package = 'base', lib.loc = NULL,
      exclude = c('chtml', 'data', 'help', 'html', 'latex', 'libs',
                  'man', 'Meta', 'po', 'R', 'R-ex', 'src'),
      include = NULL, pattern = NULL, recursive = FALSE)
}
\arguments{
  \item{package}{
    character string naming a locally installed package.  If 'package'
    is not locally installed, it is an error.
  }
  \item{lib.loc}{
    a character vector with path names of R libraries, or 'NULL'.  The
    default value of 'NULL' corresponds to all libraries currently
    known.  If the default is used, the loaded packages are searched
    before the libraries.
  }
  \item{exclude}{
    either NULL or a character vector naming subdirectories of 'package'
    to exclude from the list.  If 'include' is not NULL, 'exclude' is
    ignored.
  }
  \item{include}{
    either NULL or a character vector naming subdirectories of 'package'
    to exclude from the list.  If 'include' is not NULL, 'exclude is
    ignored.
  }
  \item{pattern}{
    an optional regular expression passed with the results of
    \link{system.file} to \link{dir}.  Only file names which match the
    regular expression will be returned.  This is ignored if 'recursive'
    is FALSE.
  }
  \item{recursive}{
    logical. Should the listing recurse into subdirectories?
  }
}
\details{
  1.  \code{fullPath <- system.file(package = package, lib.loc = lib.loc)}

  2.  \code{Dir <- dir(fullPath)}

  3.  Restrict Dir only to 'include' if provided and to all but
  'exclude' otherwise.

  4.  If recursive, return a list produced by \code{\link{dir}} for each
  of the subdirectories of interest determined in step 3.  Else, return
  only the list of subdirectories from step 3.

}
\value{
  If recursive, a list of the contents of the subdirectories of
  interest.  Else, a character vector of the names of the relevant
  subdirectories.
}
%\references{
%}
\author{
  Spencer Graves
}
\seealso{
  \code{\link[base]{system.file}}
  \code{\link{dir}}
  \code{\link[base]{file.path}}
}
\examples{
package.dir() # 'demo'
package.dir(recursive = TRUE) # contents of 'demo'
package.dir('nlme') # 'mlbook', 'scripts'
}
\keyword{package}
