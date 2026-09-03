% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/Rmpfr/inst/NEWS.Rd"))
\name{NEWS}
\title{Rmpfr News}
\encoding{UTF-8}

%% -- not yet in tar ball --
      %% \item Start vignette about \code{gamma()} inaccuracies and possible
      %% remedies (is \bold{unfinished}).

\section{Changes in version 1.1-2 [2025-10-21, r437]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{\link{.mpfrSizeof()}} providing the small integer
      sizes of the platform's MPFR library configuration.
    }
  }
 \subsection{BUG FIXES}{
    \itemize{
      \item For \code{mpfr(x, precBits=*)} it is now documented that
      \code{precBits} must not be larger than largest integer, i.e.,
      \eqn{2^31 - 1 = 2147483648}.   An error message for larger values
      is now easier to understand.

      \item Correct \file{configure.ac}, \file{Makevars}, etc needed when
      MPFR exponent size does not match limb size (as on r-universe Windows),
      fixing R-forge Rmpfr bug #6865 by \I{Mikael Jagan}.
    }
  }
}


\section{Changes in version 1.1-1 [2025-07-18, r430]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item added \code{dchisq()} for completeness, simply calling
      \code{dgamma(x, df/2, scale = 2)}.
    }
  }
 \subsection{BUG FIXES}{
    \itemize{
      \item \code{dgamma()} accidentally used double precision instead of
      "mpfr" in some cases.
    }
  }
}

\section{Changes in version 1.1-0 [2025-05-08, r426]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{pgamma()}, simply via \code{igamma()}, to become
      more visible to probability and statistics.\cr
      Note that in some regions of its domain, it is currently severely
      limited (both slow and inaccurate)
      as by its implementation in MPFR, see the \sQuote{Warning} in \code{?igamma}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item The default for our \code{dbinom()}'s \code{useLog} argument
      was determined alone by \code{any(abs(x) > 1e6)}; it is now also
      true, for smaller \code{x} and large \code{size} or extreme \code{prob}, e.g., in
      \code{dbinom(7L, mpfr(1e10, 96), prob = 1/4)}.

      Notably a bad thinko has been fixed which lead to only double precision accuracy
      results for \code{useLog = TRUE}.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item Requiring MPFR >= 3.2.0 such that \code{gamma_inc()} is
      available for our \code{igamma()} and \code{pgamma()}.
    }
  }
}

\section{Changes in version 1.0-0 [2024-11-15, r422]}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item The low-level bug fix about conversion has slightly changed the
      ABI and the low-level \code{validity} method on Windows: Non-finite
      mpfr-numbers, e.g., \code{NA}, \code{Inf}, now have different
      mantissa slot \code{@ d}, now starting with \code{0L} instead of
      previously \code{-1L}.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{num2bigq(x)} finds \dQuote{small} denominator
      \code{bigq} / \sQuote{bigRational} approximation to numeric or <mpfr> \code{x}.
      Basically a \CRANpkg{MASS}\code{::fractions()} version for
      \CRANpkg{Rmpfr} and \CRANpkg{gmp}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{format(mpfr())} no longer warns.

      \item in \code{formatDec()}: use \code{p*log(2, 10)} instead of
      \code{log(2^p, 10)}!

      \item low-level C fixes in \file{src/convert.c} thanks to analysis
      and proposition by Mikael Jagan.  Entails update of \code{validity()}
      for \code{"mpfr1"}.

      \item add \code{\\link[pkg]{...}} where needed.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item new \file{../tests/special-fun-dgamma.R} file; partly from \file{..../special-fun-ex.R}:
      accuracy checking for more accurate \code{stirlerr()} in \R 4.4.0 and later.
    }
  }
}

\section{Changes in version 0.9-5 [2024-01-20, r407]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \code{.mpfr2bigq(m)} transforms <mpfr>-number vectors to
      big rational aka \code{bigq} ones (from package \CRANpkg{gmp}).

      \item New low-level utility functions \code{.mpfr2d()} and \code{.mpfr2i()}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{x == y} and other comparisons of two \code{"mpfr"}
      operands now return \code{NaN} when an operand is \code{NaN} (mpfr
      numbers are never \code{NA}).

      \item \code{dbinom(x, *)} and \code{dnbinom(x, *)} now use the
      precision of \code{x} when it is an (integer valued) \code{"mpfr"}
      object, and then notice or even error when precision is lost as
      \code{x} is coerced to integer.

      \item (erange-related comments in examples)

      \item \code{chooseMpfr(a, n)} gives \code{0} when \code{n < 0}, as
      \R's \code{choose()} and when \code{a} is integer valued, may use
      the equivalent of \code{chooseMpfr(a, a-n)} when \code{a-n < n} to be
      more accurate and faster.

      \item \code{.mpfr2bigz(m)} now also works for really large \code{m}.
    }
  }
}

\section{Changes in version 0.9-4 [2023-12-04, r399]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed Windows-only (long = 32 bit) bug; adapt mpfr1-validity() check.
      \item allow `.Platform$endian != "little"  to 'work' in validity().
      \item format ("\%ld" etc) fixes in \file{src/utils.c}
      \item fix "lost braces"
    }
  }
}

\section{Changes in version 0.9-3 [2023-07-27, r394]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{any()} and \code{all()} no longer damage \R's internal
      FALSE or TRUE, fixing R-forge Rmpfr bug #6764 by Andrew Manderson.
    }
  }
}

\section{Changes in version 0.9-2 [2023-04-21, r392]}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{formatMpfr()} and hence the \code{format()} method for
      \code{"mpfr"} etc now uses a \code{scientific = <number>} the same as
      \pkg{base} \code{format.default()};  accidentally it had the
      \emph{sign} of the number interpreted differently and failed to use
      e.g., \code{getOption("scipen")} correctly.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{log(<mpfr>, base=10)} now works.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{igamma()}, \code{hypot()}, and \code{atan2()} no longer
      only return 53 bit (or more generally, mpfr default precision)
      results.

      \item New \code{max2_prec()} utility instead of \code{imax2()},
      also curing compilation warning.

      \item Internal (exported) objects are now (somewhat) documented, as
      eventually required by \R.
    }
  }
}

\section{Changes in version 0.9-1 [2023-01-30, r386]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{ldexpMpfr(f, E)} now returns an \code{"mpfr"} vector
      (instead of a \code{"mpfr1"}).
    }
  }
}

\section{Changes in version 0.9-0 [2023-01-16, r382]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Add the \code{extendInt="*"} option to \code{unirootR()}
      adopting its addition in base \R \code{uniroot()}.
      \cr
      Also add an optional \code{verbDigits} argument only used in case
      \code{verbose} is true.

      \item New \code{qnormI()}, the normal quantile function, computed via
      inversion of (our potentially arbitrarily accurate) \code{pnorm()},
      using \code{unirootR()} root finding.

      \item \code{sapplyMpfr()} gets \code{drop_1_ = TRUE} option,
      \emph{changing} default behaviour to return an \code{"mpfr"} vector
      instead of a 1-column \code{"mpfrMatrix"}.
    }
  }
  \subsection{BUG FIXES}{
     \itemize{
       \item \code{sapplyMpfr()} now gets the correct dimension when it returns
       an \code{"mpfrMatrix"} (or \code{"*array"}).

       \item Fix \code{\\eqn{(\\gamma)}{}} - as requested by KH (June 9).

       \item \code{.mpfr_gmp_numbbits()} is documented now (notably as I
       see it is only 32 on Winbuilder (!)).
    }
  }
}

\section{Changes in version 0.8-9 [2022-06-02, r363]}{
  \subsection{BUG FIXES}{
     \itemize{
      \item Embarrassing thinko in \code{dpois} (in \code{useLog=TRUE}
      case) fixed.
    }
  }
}

\section{Changes in version 0.8-8 [2022-06-01, r362]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Our \dQuote{mpfr-ized} \code{all.equal()} is now based on
      hidden \code{all.equalMpfr()} which uses a smart default
      \code{tolerance} for all methods, and shows much less digits if there
      are differences (compatibly to \code{all.equal.numeric()}).
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Our \code{dpois(x, lambda, *)} now works \R-compatibly for
      \code{lambda=+Inf} giving 0 or \code{-Inf} when \code{log=TRUE}.

      \item \code{formatMpfr(x, scientific=FALSE)} now \dQuote{works}.
      % tweaking default of 'maybe.full': true also if(isFALSE(scientific))
    }
  }
  \subsection{Misc}{
    \itemize{
      \item \file{man/*.Rd}: get rid of some \code{\\} escapes (needed in past; wrong now)
    }
  }
}

\section{Changes in version 0.8-7 [2021-10-27, r353]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item In \file{Ops.c}'s \code{R_mpfr_mod()}, no longer allocate
      but never use nor free \code{rr} - fixing YAL (yet another leak).
    }
  }
}

\section{Changes in version 0.8-6 [2021-10-08, r351]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item In \file{convert.c}'s \code{R_mpfr_frexp()}, do
      \code{mpfr_clear(*)}, fixing a valgrind-detectable leak.

      \item additionally use single mpfr_exp_t and its pointer (and assign
      to R allocated vector inside main loop, fixing a valgrind
      "unitialized value" case.
    }
  }
}

\section{Changes in version 0.8-5 [2021-10-05, r349]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New mpfr-ized \code{dt(x, df, ..)}; noncentrality \code{ncp}
      not yet supported.

      \item New arithmetic functions \code{frexpMpfr()} and
      \code{ldexpMpfr()}, corresponding to C's (and CRAN package
      \CRANpkg{DPQ}'s) \code{ldexp()} and \code{frexp()} functions.

      \item \code{sapplyMpfr()} now also returns \code{"mpfrMatrix"} or
      \code{"mpfrArray"} when appropriate.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Improved \code{`[.Ncharacter`} method for \code{formatHex()} etc.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item No longer include \file{Rdefines.h} as it is somewhat deprecated.

      \item update both \file{configure\{,.ac\}} (for \command{autoconf 2.71}).
    }
  }
}

\section{Changes in version 0.8-4 [2021-03-24, r341]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{dpois(x, *)} now also gets an explicit optional argument
      \code{useLog} instead of just switching to log-scale by the same
      test as defines \code{useLog}'s default.

      \item The default method of \code{mpfr(r)} now also works when
      \code{r} is a \code{list} of \code{"mpfr1"} objects, e.g., resulting
      from \code{Vectorize(.)} or similar applied to mpfr-vectors.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item fixed problem detected by clang-UBSAN testing in \file{src/Ops.c}
      (from arithmetic bug fix in 0.8-3). % ../tests/arith-ex.Rout_clang-UBSAN
    }
  }
}

\section{Changes in version 0.8-3 [2021-03-22, r340]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{dnbinom(x, *)} and \code{dbinom(x, *)} now also work when
      \code{x} is too large to be coerced to \code{integer} \emph{and} they
      get a new optional argument \code{useLog} (with a smart default) in
      order to choose log-scale computation also for \code{log=FALSE}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item For arithmetic (\code{`+`}, \code{`*`}, etc), the check to see
      if a numeric can be validly coerced to a long has been amended such
      that 9223372036854775808 is no longer accidentally coerced into
      \emph{negative}.

      \item \code{dpois(x, lambda)} now works via log-scale in case
      \code{exp(-lambda)} or \code{lambda^x} would under- or overflow even
      for mpfr-numnbers.

      Analogously, \code{dbinom()} and \code{dnbinom()} work via log-scale
      in case even mpfr-arithmetic would overflow.
    }
  }
}

\section{Changes in version 0.8-2 [2020-11-10, r337]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The workhorse for all matrix multiplication (\code{\%*\%},
      \code{crossprod()}, and \code{tcrossprod()}), i.e., \code{.matmult.R(x, y, *)}
      (not exported), gets new optional arguments \code{fPrec = 1} and
      \code{precBits}, which defaults to the maximum of \code{getPrec(x)}
      and \code{getPrec(y)}.
      To get correct \code{crossprod()} and tcrossprod() generics (with a
      formal \code{\dots} argument), need a new \CRANpkg{gmp} release, as
      we get these two generics from \pkg{gmp}.

      \item \code{matmult(x,y)} is identical to \code{x \%*\% y}, but
      \code{matmult()} has further optional arguments \code{fPrec} and
      \code{precBits}, see above.

      \item New \code{is.mpfr(x)} function; simple, with fast pre-test.

      \item In \code{pbetaI(q, ..)}, when \code{q} is a bigrational, i.e.,
      of class \code{"bigq"} (package \CRANpkg{gmp}), the computations are
      \emph{exact} now, using big rational arithmetic.

      \item New \code{dnbinom()} function (with corresponding new
      \code{\link{conflicts}()} with the \pkg{stats} package base \R
      function), providing an \code{"mpfr"}-number version of the negative
      binomial probabilities.

      \item \code{.mpfr_erange_set()} now can set \emph{both} exponent
      range limits simultaneously, and now returns invisibly \code{TRUE} if
      the change succeeded.

      \item New \code{log1mexp()} and \code{log1pexp()}, as "from"
      \CRANpkg{copula}, notably as the vignette has been here, and authored
      by me.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{mpfr(mm)} (and similar) now works for \code{"bigq"} or
      \code{"bigz"} matrices, thanks to a report by Jerry Lewis.

      \item Arithmetic and other \code{"\link{Ops}"} between \code{"mpfr"}
      and \code{"bigq"} aka bigrational numbers, now use the
      \dQuote{inherent} precision of the bigrational.

      \item \code{chooseMpfr(a, n)} and similar now \dQuote{work} when
      \code{n} is of length zero \emph{and} when an experimental version of
      \code{stopifnot(*, allow.logical0=FALSE)} is used.

      \item Our \code{cbind()} and \code{rbind()} methods with
      \code{signature = "mNumber"} now keep and construct column and row
      names as the corresponding base functions, obeying
      \code{deparse.level}.

      \item Fixed "not-yet"-as-cran \file{NOTE}
      \emph{Undeclared packages \pkg{dfoptim}, \pkg{pracma} in Rd xrefs}.
    }
  }
}

\section{Changes in version 0.8-1 [2020-01-14, r323]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Provide dummy C function in the case MPFR library is older than
      3.2.0, thanks to Brian Ripley.  Should work around check \code{NOTE}s.
    }
  }
}

\section{Changes in version 0.8-0 [2019-12-05, r321]}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Formatting incl \code{print()}ing by default uses a \code{"+"}
      after the exponential character (\code{"e"} by default).
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Provide new low-level utilities \code{.mpfr_formatinfo()},
      \code{.mpfr2exp()}, and \code{.mpfr_erange_is_int()}.

      \item Renamed low-level utility functions to use \code{_} instead of
      \code{.} keeping initial \code{"."}, e.g.,
      \code{.mpfr.gmp.numbbits()} to \code{.mpfr_gmp_numbbits()}.

      \item \code{formatMpfr()} gets a new optional \code{decimal.plus = TRUE} which
      adds a \code{"+"} before positive exponents when in exponential (aka
      \dQuote{scientific}) representation.

      The \code{mpfr} and \code{mpfrArray} \code{print()} methods get
      a corresponding \code{decimal.plus} argument with a default that can
      be set by \code{options(Rmpfr.print.decimal.plus = *)} to allow
      strict back compatibility where needed.

      \item For MPFR (C library) version >= 3.2.0 (not by default in Fedora 30!),
      provide the \emph{incomplete} gamma function \code{igamma(a,x)}
      which is closely related to \code{pgamma(x,a)}, see help page.% and DO
      % provide an MPFR pgamma() version, too
      %  ==> then even MPFR pchisq() with *non*centrality maybe!

      \item Now also export S3 method of \code{unique()} for \code{"mpfr"},
      such that \pkg{base} \code{factor(<mpfr>)} \dQuote{works}.% -> ../man/mpfr-class.Rd
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{formatMpfr()} and hence all \code{print()}ing suffered
      from an integer overflow bug with very large (base 2) exponents.

      \item \code{.mpfr2str(x, *)} is no longer dependent on the
      \emph{order} of the elements in \code{x}; consequently
      \code{format()} and \code{print()} may use less digits in case the
      precision decreases along \code{x}.

      \item adapt to new C compiler default behavior \code{-fno-common},
      using \code{extern #include <Syms.h>} in most \file{*.c} files.
    }
  }
}

\section{Changes in version 0.7-3 [2019-08-06, r305]}{
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item Decreased the default for \code{max.digits} to 999.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Provide \code{dgamma(x, shape)} version, e.g., for small shape
      parameter where most of the mass is on very small \code{x} not
      representable as double precision numbers.

      \item Low-level formatting function \code{.mpfr2str()} --- called by
      \code{format()} and hence \code{print()} methods for \code{"mpfr"}
      objects --- now finally obeys its \code{maybe.full} argument when it
      is \code{FALSE}, internally in C's \code{mpfr2str()}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{pnorm(<large>, log.p=TRUE)} no longer underflows much too
      early, thanks to reports by Jerry Lewis.

      \item \code{print.mpfrArray()} now also uses a finite
      \code{max.digits} default, preventing, e.g., \code{cbind(x, y)} to
      use too many digits.
    }
  }
}

\section{Changes in version 0.7-2 [2019-01-18, r299]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{str(<mpfr>)} no longer calls \code{formatMpfr(x, digits, *)}
      with a \code{digits} \emph{vector} of the same length as \code{x}
      (which never worked correctly).

      \item \code{seqMpfr(1, length.out=8)} now works correctly.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{unirootR()} gets an option to \emph{not} warn on
      non-convergence.

      \item Provide a \code{summary()} method for \code{"mpfr"} numbers
      closely modeled after \code{summary.default} for numeric.

      \item \code{mpfr(NULL)} now works, returning \code{mpfr(logical())}.

      \item a simple \code{sapplyMpfr()} function, showing how to work
      around the fact that \code{sapply()} does typically not work with
      \code{"mpfr"} numbers.
    }
  }
}


\section{Changes in version 0.7-1 [2018-07-24, r291]}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{formatMpfr()} for large low-precision numbers now uses
      scientific representation, fixing the bug RMH was reporting March
      17 already.

      \item \code{outer()} is \dQuote{imported} from base, so it
      behaves as an \pkg{Rmpfr} function which dispatches e.g., when
      calling \code{\link{tcrossprod}()}.
    }
  }
}

\section{Changes in version 0.7-0 [2018-01-12, r284]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{.mpfr2list()} and \code{mpfrXport()} gain an option
      \code{names} (for nicer output).

      \item \code{formatMpfr()} and the \code{print()} method get a new
      option \code{max.digits} with default \code{9999} for the print
      methods, to limit the number of digits printed for high precision
      numbers.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item For non-\dQuote{regular} mpfr numbers, the \code{d} slot in
      the \code{"mpfr1"} representation is now empty instead of
      \dQuote{random}.   This also eliminates valgrind warnings about
      uninitialized values in C.
    }
  }
}
\section{Changes in version 0.6-2 [2017-05-29, r264], never on CRAN}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The S3 classes \code{"Hcharacter"} and \code{"Bcharacter"}
      resulting from \code{formatHex()} and \code{formatBin()} now
      \dQuote{inherit from} \code{"character"} formally.

      \item They also got a \code{`[`} method, so subsetting should
      work, including for \code{array}s of these.

      \item The \code{"mpfr"} method of \code{str()} gains option
      \code{internal}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item when \code{print()}ing mpfr numbers, the result no longer sometimes
      loses the last digit.

      \item \code{dnorm()} now works correctly with mpfr numbers;
      similarly \code{dbinom()} and \code{dpois()} should work in all
      cases, now.

      \item in \file{NAMESPACE}, also \code{exportMethods(apply)}, so
      \pkg{SNscan} works.

      \item \code{print(formatHex(..))}, \code{formatBin()} and
      \code{formatDec()} now look better and are more correct; the first
      two get a new option \code{expAlign} indicating to use the same
      number of digits for exponents (in \dQuote{scientific} cases).

      Notably, \code{mpfr(formatBin(mpx))} works for more \code{mpx}
      objects (of class \code{"mpfr"}).

      \item \code{format(mpfr(3,7), digits = 1, base = 2)} no longer
      crashes (from inside MPFR).

      \item \code{formatDec(mpfr(NA, 7))} now works.

      \item For non-\dQuote{regular} mpfr numbers, the \code{d} slot in
      the \code{"mpfr1"} representation is now empty instead of
      \dQuote{random}.   This also eliminates valgrind warnings about
      uninitialized values in C.
    }
  }
}

\section{Changes in version 0.6-1 [2016-11-15, r249]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{head()} and \code{tail()} methods for \code{"mpfrMatrix"}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item C-level \code{mpfr2str()} no longer calls S_realloc() with
      wrong "old size" (thanks to Bill Dunlap).
      \item \code{c()} now also works when its result is a length-0
      \code{"mpfr"} object.
    }
  }
}

\section{Changes in version 0.6-0 [2015-12-04, r237]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{mpfr()} now is S3 generic with several methods, notably a
      \code{"mpfr"} method allowing to change precision or rounding mode.

      \item \code{mpfr()}, \code{formatMpfr()}, etc, now work with bases from 2
      to 62 (using digits, upper and lower case ASCII letters, \code{62
      == 10 + 2*26} characters), as this has been in MPFR since version
      3.0.0 (see \code{\link{mpfrVersion}}), which is hence (implicitly)
      required for \code{base} greater than 36.

      \item \code{formatMpfr()} gets a new argument \code{base = 10} and can
      be used to produce in other bases, notably binary (\code{base = 2})
      or hexadecimal (\code{base = 16}).

      \item \code{str(<mpfr>, ....)} is now based on \code{formatMpfr()}
      and nicely shows numbers also out of the double precision range.

      Further, it now chooses a smart default for optional argument
      \code{vec.len}.

      \item \code{matrix(mp, ..)} now also works when \code{mp} is of
      class \code{"mpfr"}.
      \item new matrix \code{norm()} for several \code{kind}s.

      \item new functions \code{formatHex()} and \code{formatBin()}
      thanks to Rich Heiberger.

      \item \code{mpfr(x)} also works as \emph{inverse} of
      \code{formatBin} and \code{formatHex}.

      \item \code{roundMpfr()} and mathematical functions such as
      \code{jn}, or \code{chooseMpfr()} get new optional argument
      \code{rnd.mode} passed to the corresponding MPFR function.

      \item \code{median(x)}, \code{mean(x, trim)} for \code{trim > 0}
      now work fine for \code{"mpfr"} x, and \code{quantile(x, *)} no
      longer needs \code{names=FALSE} to avoid a warning.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{pnorm(.)}, \code{j0()} and similar special functions
      now preserve \code{mpfrMatrix} and \code{mpfrArray} classes.

      \item similarly, \code{is.finite()} etc keep the
      \code{dim()}ensionality for \code{"mpfrArray"} arguments.

      \item \code{mpfr("0xabc", base=16)} and \code{mpfr("0b101",
	base=2)} and corresponding \code{getPrec()} now give the correct
      precBits instead of too many.

      \item \code{str(<0-length mpfr>)} now works correctly.
    }
  }
}



\section{Changes in version 0.5-7 [2014-11-27, r205]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{as.integer()} now rounds \dQuote{to zero} as for
      regular \R numbers (it accidentally did round \dQuote{to nearest}
      previously).
    }
  }
}

\section{Changes in version 0.5-6 [2014-09-05, r203]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item experimental \code{mpfrImport()}, \code{mpfrXport()}
      utilities -- as we found cases, where save() \code{"mpfr"} objects
      were \emph{not} portable between different platforms.

      \item \code{as(*,"mpfr")} now also supports rounding mode
      \code{"A"} (\dQuote{\bold{A}way from zero}).

      \item Several hidden low level utilities also get a
      \code{rnd.mode} option.
    }
  }
}

\section{Changes in version 0.5-5 [2014-06-19, r190]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The result of \code{integrateR()} now prints even if a
      warning happened.

      \item \code{pbetaI(x, a,b)}, the arbitrarily accurate
      \code{pbeta()} computation for \emph{integer} \eqn{a} and \eqn{b},
      now works for larger \eqn{(a,b)}.

      \item Newly providing \code{mpfr}-versions of \code{dbinom()},
      \code{dpois()}, and \code{dnorm()}.

      \item New utility functions \code{mpfr_default_prec()},
      \code{.mpfr.minPrec()}, etc, to get, check, set default exponent
      ranges and precision.

      \item New \code{sinpi()}, \code{cospi()} etc, notably for \R >= 3.0.1.
    }
  }
}

\section{Changes in version 0.5-4 [2013-10-22, r173]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}


%% The first CRAN release
\section{Changes in version 0.1-5 [2009-08-06]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item First CRAN release on 'Publication:' 2009-08-14 20:24:02

      \item new \code{pmin()} and \code{pmax()}, improving \code{seqMpfr()}.

      \item new \code{"Math"} and \code{"Math2"} group methods, notably for
      \code{round()} and \code{signif()}.

      \item \code{as(. , "integer")} now works (via C \code{mpfr2i}).
    }
  }
  \subsection{More details for old versions up to Feb. 2015:}{
    \itemize{
      \item See file \file{ChangeLog} %% ../ChangeLog <<<
    }
  }
}
