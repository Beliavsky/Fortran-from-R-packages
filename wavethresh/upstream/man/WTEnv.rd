\name{WTEnv}
\alias{WTEnv}
\title{Environment that exists to store intermediate calculations for re-use within the same R session.}
\description{Environment that stores results of long calculations so that they can be made available for immediate reuse.
} 
\details{This environment is created on package load by wavethresh. The results of some intermediate calculations get stored in here (notably by
\code{\link{PsiJ}}, \code{\link{PsiJmat}} and \code{\link{ipndacw}}). The
reason for this is that the calculations are typically lengthy and it saves
wavethresh time to search the \code{WTEnv} for pre-computed results.
For example, \code{\link{ipndacw}} computes matrices of various orders.
Matrices of low order form the upper-left corner of matrices of higher order
so higher order matrix calculations can make use of the lower order instances.

A similar functionality was present in wavethresh in versions 4.6.1 and
prior to this. In previous versions computations were saved in the users
current data directory. However, the user was never notified about this
nor permission sought. 

The environment \code{WTEnv} disappears when the package disappears
and the R session stops - and results of all intermediate calculations
disappear too. This might not matter if you never use the larger objects
(as it will not take much time to recompute). 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{ipndacw}}, \code{\link{PsiJ}}, \code{\link{PsiJmat}}}

\examples{
#
# See what it is
#
WTEnv
#<environment: 0x102fc3830>
#
# Compute something that uses the environment
#
fred <- PsiJ(-5)
#
# Now let's see what got put in
#
ls(envir=WTEnv)
#[1] "Psi.5.10.DaubLeAsymm"
}
\keyword{algebra}
\author{G P Nason}
