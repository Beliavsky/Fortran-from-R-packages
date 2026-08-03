% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/CLA/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{CLA}}% MM: look into ../svn-log-from.all
\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.96-3 (2024-07-29)}{
  \subsection{Tweaks}{
    \itemize{
      \item \file{../DESCRIPTION} now uses \code{Authors@R}.
    }
  }
}

\section{Changes in version 0.96-2 (2021-12-14)}{
  \subsection{Tweaks}{
    \itemize{
      \item \file{../tests/SP500-ex.R} tweaks to pass on 'M1mac' (Mac "Big Sur" M1 proc. arm64)
      \item \file{../tests/SP500-ex.R}: new main tweak: set all weights < 1e-10 to zero .
    }
  }
}

\section{Changes in version 0.96-1 (2020-09-19)}{
  \subsection{Tweaks}{
    \itemize{
      \item Should finally work on Solaris 32-bit
    }
  }
}

\section{Changes in version 0.96-0 (2020-09-06)}{
  \subsection{New Features}{
    \itemize{
      \item \code{CLA(mu, covar, ..)} gets new option \code{check.cov} which by default
      is \code{TRUE} and hence checks if the covariance matrix \code{covar}
      is positive definite.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{CLA()} for small cases may have "empty" new lambda even
      when there free weights left.  This is now caught, with a warning unless
      new option \code{check.f} is set to false.  See new example \code{mc3}.

      \item the platform dependent tests in \file{../tests/SP500-ex.R} now
      should also work on Solaris 32-bit.
    }
  }
}

\section{Changes in version 0.95-1 (2019-05-09)}{
  \subsection{Testing}{
    \itemize{
      \item modularized partly platform dependent testing into
      \code{moreSessionInfo()} in \file{./xtraR/platform-sessionInfo.R};
      updated to newest version from robustbase (so to detect KH's debian
      OpenBLAS setup).
    }
  }
}

\section{Changes in version 0.95-0 (2019-03-13)}{
  \subsection{New Features}{
    \itemize{
      \item \emph{Not} back-compatible \emph{changed} result: The very
      first \dQuote{turning point} from the algorithm with (\code{lambda
	= NA} etc) is no longer stored as part of the result, compatibly
      with the litterature and other implementations.

      \item Added data set \code{muS.10ex} orginally from Markowitz and Todd.

      \item Added \file{NEWS} file.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item needed \code{, drop=FALSE} in internal \code{getMatrices()},
      e.g., for the \code{*10ex} example.

      \item Test works with old R 3.3.x.

      \item In some cases such as \code{muS.10ex}, the last turning
      point with \code{lambda=0} was not computed previously.

      \item adapt \file{man/findSig.Rd} example to R 3.6.0's new default
      \code{sample()} kind.

      \item Run the \file{tests/*.R} with R versions that use one
      OpenBLAS, ATLAS, or MKL (and skip some tests on these).
    }
  }
}

\section{Changes in version 0.90-1 (2018-02-05)}{
  \subsection{New Features}{
    \itemize{
      \item Add nice \code{plot()} method base on Shi's \code{MS_plot()}
      \item Literature references in \code{?CLA}
      \item Vectorized \code{findSig()} and \code{findMu()}.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{findSig.1()} now works.
      \item test also works in \R configured with \command{--no-long-double}.
    }
  }
}

\section{Changes in version 0.90-0 (2018-01-15)}{
  \subsection{New Features}{
    \itemize{
      \item First release to CRAN
    }
  }
}

