\name{runscript}
\alias{runscript}
\title{ Run a package script }
\description{
  Run a script associated with a particular chapter
}
\usage{
runscript(x, method = c('run', 'copy', 'view', 'show', 'dir'),
          ask = TRUE, fmt = "ch\%02d.R", package = "FinTS",
          subdir = "scripts", lib.loc = NULL)
}
\arguments{
  \item{x}{
    an object to identify a file in package/subdir via
    \code{sprintf(fmt, x)}. 

    For example, the default 'fmt' translates \code{x = 2} into 'ch02.R'.  If
    no 'x' is specified, a directory of options is provided.

    CAUTION:  Under some systems like ESS (Emacs Speaks Statistics)
    under Windows, pop-up menus such as produced by \code{runscript()} may
    not work properly.
  }
  \item{method}{
    One of the following:

    \describe{
      \item{run}{
	run the desired script file, similar to \link{demo} or
	\link{example}.
      }
      \item{copy}{
	make a copy if the desired script file in the working directory,
	similar to \code{Stangle(vignette(...)[["file"]])}.
      }
      \item{view}{
	display the desired script file on R console but do not execute it.
      }
      \item{show}{
	display the desired script file using \link{file.show}
      }
      \item{dir}{
	return the directory showing only the location of the desired
	script.
      }
    }

    Partial matching is allowed.
  }
  \item{ask}{
    logical:  Should \code{par(ask=TRUE)} be called before graphical
    output happens from the script?
  }
  \item{fmt}{
    a format to be used with 'x' in \link{sprintf} to create the name of
    a file in lib.loc/package/subdir.
  }
  \item{subdir}{
    subdirectory of package containing a file of the name constructed
    via \code{sprintf(fmt, x)}.
  }
  \item{package}{
    Name of a package with subdirectory 'subdir'.
  }
  \item{lib.loc}{
    NULL or character string identifying the location where
    \verb{system.file(subdir, package, lib.loc)} will find the folder
    containing the file identified via \code{sprintf(fmt, x)}.
  }
}
\details{
  similar to \link{demo} or \link{example}
}
\value{
  the full path and filename, invisibly unless \code{method == 'dir'}
}
%\references{
%}
\author{
  Gabor Grothendieck and Spencer Graves
}
\seealso{
  \code{\link[utils]{demo}},
  \code{\link[base]{sprintf}},
  \code{\link[base]{system.file}},
  \code{\link{package.dir}},
  \code{\link{Stangle}},
  \code{\link{vignette}},
  \code{\link{example}}
}
\examples{
\dontrun{
# provide a menu
runscript()

# run R/library/FinTS/scripts/ch01.R
runscript(1)

# same as:
runscript(1, "run")

# make a copy as 'ch01.R' in the working directory
runscript(1, 'copy')

# display on console only
runscript(1, 'view')

# display using file.show
runscript(1, 'show')

# where is it?
runscript(1, 'dir')

# run R/library/nlme/scripts/afda-ch01.R
if(require(fda))
    runscript(1, fmt = "afda-ch\%02d.R", package = "fda")
}
}
\keyword{utilities}
