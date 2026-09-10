\name{accessC.wd}
\alias{accessC.wd}
\title{Get smoothed data from wavelet object (wd)}
\description{
  The smoothed and original data from a wavelet decomposition structure
  (returned from \code{\link{wd}}) are packed into a single vector in that
  structure.  This function extracts the data corresponding to a
  particular resolution level.
}
\usage{
\method{accessC}{wd}(wd, level = nlevelsWT(wd), boundary=FALSE, aspect, \dots)

}
\arguments{
  \item{wd}{
    wavelet decomposition structure from which you wish to extract the
    smoothed or original data if the structure is from a wavelet decomposition,
    or the reconstructed data if the structure is from a wavelet reconstruction.
  }
  \item{level}{
    the level that you wish to extract.  By default, this is the level with
    most detail (in the case of structures from a decomposition this is the
    original data, in the case of structures from a reconstruction this is the
    top-level reconstruction).
  }
  \item{boundary}{logical;
    if \code{TRUE} then all of the boundary correction values
    will be returned as well (note: the length of the returned vector
    may not be a power of 2).\cr
    If \code{boundary} is false, then just the coefficients will be returned.

    If the decomposition (or reconstruction)
    was done with periodic boundary conditions, this option has no effect.}
  \item{aspect}{Applies a function to the coefficients before return.
	Supplied as a text string which gets converted to a function.
	For example, "Mod" for complex-valued arguments}
\item{\dots}{any other arguments}
}
\value{A vector of the extracted data.
}
\details{
The \link{wd} (\code{\link{wr.wd}}) function produces a wavelet decomposition (reconstruction)
structure.

For decomposition, the top level contains the original data, and
subsequent lower levels contain the successively smoothed data.
So if there are \eqn{2^m} original data points, there will be m+1 levels
indexed 0,1,\dots{},m. So

  \code{accessC.wd(wdobj, level=m)}

pulls out the original data, as does

  \code{accessC.wd(wdobj)}

To get hold of lower levels just specify the level that you're interested
in, e.g.

  \code{accessC.wd(wdobj, level=2)}

gets hold of the second level.

For reconstruction, the top level contains the ultimate step in the
Mallat pyramid reconstruction algorithm, lower levels are intermediate
steps.

The need for this function is a consequence of the pyramidal structure
of Mallat's algorithm and the memory efficiency gain achieved by
storing the pyramid as a linear vector. AccessC obtains information about
where the smoothed data appears from the fl.dbase component of
an \link{wd.object}, in particular the array \code{fl.dbase$first.last.c} which
gives a complete specification of index numbers and offsets for
\code{wd.object$C}.

Note that this function is method for the generic function \code{\link{accessC}}.
When the \code{\link{wd.object}} is definitely a wd class object then you only need use the generic version of this function. 

Note that this function only gets information from \code{\link{wd}} class
objects. To insert coefficients etc. into \code{\link{wd}} structures you have to
use the \code{\link{putC}} function (or more precisely, the \code{\link{putC.wd}} method).
}

\references{
Mallat, S. G. (1989) A theory for multiresolution signal decomposition: 
the wavelet representation. \emph{IEEE Transactions on Pattern Analysis
 and Machine Intelligence} \bold{11}, 674--693.

Nason, G. P. and Silverman, B. W. (1994). The discrete wavelet transform in S. 
\emph{Journal of Computational and Graphical Statistics,} \bold{3}, 163--191. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994}
\seealso{ \code{\link{wr}}, \code{\link{wd}}, \code{\link{accessD}}, \code{\link{accessD.wd}},
  \code{\link{filter.select}}, \code{\link{threshold}}, \code{\link{putC.wd}}, \code{\link{putD.wd}}.}

\examples{
## Get the 3rd level of smoothed data from a decomposition
dat <- rnorm(64)
accessC(wd(dat), level=3)
}
\keyword{manip}
\author{G P Nason}
