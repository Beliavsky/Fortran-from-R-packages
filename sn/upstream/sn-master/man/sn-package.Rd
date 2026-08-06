%  file sn/man/sn-package.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2021 Adelchi Azzalini
%---------------------
\name{sn-package}
\docType{package}
\encoding{UTF-8}
\alias{sn-package}
\alias{SN}

\title{Package \pkg{sn}: overview, background and history}

\description{
  The \pkg{sn} package  provides facilities to  define and manipulate 
  probability distributions of the skew-normal (\acronym{SN}) family and 
  some related ones, notably the skew-\eqn{t} (\acronym{ST}) and 
  the unified skew-normal (\acronym{SUN}) families. 
  For a number of these families, statistical methods are provided, 
  to perform data fitting and model diagnostics, 
  in the univariate and the multivariate case.
}
%----------
\section{Overview of the package structure and commands}{
A separatate document is entirely dedicated to the presentation of 
the package structure and its basic functions; 
see the \code{\link[=overview-sn]{package overview}}. 
} % end section

%----------
 
\section{Background information and references}{
The package adopts the terminology, notation and general framework of the
monograph by Azzalini and Capitanio (2014). 
This matching constitutes a reason for the numerous references to the book 
in the documentation of the package. 

An additional reason for referring to that monograph instead of the original
research papers is that the book provides a relatively not-so-formal account 
of material which has been elaborated in a number of publications, sometimes
very technical, or re-elabotated over a few papers or possibly mixing 
the information of key interest with other  material. 
In other words, the motivation behind this policy is readability,
not indulgence in self-citation.

When one or a few original sources appeared to deliver the required
information in a compact and accessible form, they have been cited directly. 
In any case, the cited sections of the book include bibliographic notes
which refer back to the original sources.
} % end section
%----------
\section{A bit of history}{%
The first version of the package was written in 1997, 
and it was uploaded on \acronym{CRAN} in 1998.
Subsequent versions have evolved  gradually up to version 0.4-18 in May 2013. 
 
In January 2014, version 1.0-0 has been uploaded to \acronym{CRAN}.
This represented a substantial re-writing of the earlier \sQuote{version 0.x},
developed in broad connection with the book by Azzalini and Capitanio (2014). 
Differences between the \sQuote{version 0} and the \sQuote{version 1} series 
are radical; they concern the core computational and graphical part as well 
as the user interface. 
Since version 1.0-0, the S4 protocol for classes and methods has been adopted.

After various versions 1.x-y, version 2.0.0 has appeared in March 2021,
providing support for the \acronym{SUN} distribution.

Additional information on the evolution of the package is provided in
\code{NEWS} file, accessible from the package documentation index page.
}


%----------

\section{Backward compatibility versus \sQuote{version 0.4-18}}{%
There is a partial backward compatibility of newer version versus
\sQuote{version 0-4.18} of the package.
Some functions of the older version would work as before
with virtually no change; a wider set arguments is now allowed. Functions
\code{\link{dsn}}, \code{\link{dst}}, \code{\link{dmsn}} and alike fall in
this category: in some cases, the names of the arguments have been altered, 
but they work as before if called with unnamed arguments; similar cases are
\code{\link{msn.mle}}, \code{\link{sn.cumulants}} and \code{\link{T.Owen}}.
Notice, however, that \code{\link{msn.mle}} and other fitting functions have
effectively been subsumed into the more comprehensive fitting function
\code{\link{selm}}.

A second group of functions will work with little or even minimal changes.
Specific examples are functions \code{sn.mle} and \code{st.mle} which have
become \code{\link{sn.mple}} and \code{\link{st.mple}}, with some additional
arguments (again, one can achieve the same result via \code{\link{selm}}). 
Another example is constitude by the group of functions \code{dp.to.cp},
\code{cp.to.dp} and  \code{st.cumulants.inversion}, which have been replaced 
by the more general functions \code{\link{dp2cp}} and \code{\link{cp2dp}}; 
one only needs to pay attention to conversion from 3rd and 4th order 
cumulants  to their standardized form in connection with the replacement of 
\code{st.cumulants.inversion}.

Finally, some functions are not there any longer, with no similarly-working
functions in the new version.  This is the case of \code{sn.mle.grouped}
and \code{st.mle.grouped} for maximum likelihood estimation from
grouped data, that is, data recorded as intervals and corresponding 
frequencies.
} % end section

\section{Requirements}{
  \R version 2.15-3 or higher, plus packages \pkg{mnormt}, 
  \pkg{numDeriv} and \pkg{quantreg},   in addition to standard
  packages (\pkg{methods}, \pkg{graphics}, \pkg{stats4},  etc.)
}

\section{Version}{
The command \code{citation("sn")} indicates, among other information,
the running version of the package.
The most recent version of the package can be obtained from
the web page: \url{http://azzalini.stat.unipd.it/SN/}
which also provides related material.

From the above-indicated web page, one can also obtain the package 
\sQuote{sn0} which is essentially the last \sQuote{version 0} 
(that is, 0.4-18) with suitable renaming of  certain  ingredients.
This allows to have both the current and the old package installed 
at the same time.
}

\section{Author}{Adelchi Azzalini.
% Dipart. Scienze Statistiche, Università di Padova, Italia.
Please send comments, error reports \emph{et cetera} to the author, 
whose web page is \url{http://azzalini.stat.unipd.it/}.
}

\section{Licence}{
This package and its documentation are usable under the terms of 
the \dQuote{GNU General Public License} version 3 or version 2,
as you prefer; a copy of them is available from 
\url{https://www.R-project.org/Licenses/}.
  
While the software is freely usable, it would be appreciated
if a reference is inserted in publications or other work
which makes use of it. For the appropriate way of referencing it, 
see the command \code{citation("sn")}. 
}

\seealso{\code{\link[=overview-sn]{package-overview}} }

\references{
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.
}

\keyword{multivariate}
\keyword{distribution}
\keyword{univar}
\keyword{regression}
\concept{skew-elliptical distribution}
\concept{skew-normal distribution}
\concept{skew-t distribution}
\concept{symmetric distribution}
\concept{symmetry-modulated distribution}
\concept{unified skew-normal distribution}
