% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("C:/projects/statproj/R/rrcov/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{rrcov}}
\encoding{UTF-8}

\section{CHANGES in rrcov VERSION 1.7-7 (2025-04-10)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item R/datadoc - emoved a pipe in an example for ionosphere data set to avoid dependence on R 4.1.0
      \item R/*: added "#' export" to many functions to satisfy roxygen2
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.7-6 (2024-08-17)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item R/plot-utils - a bug fixed in function .rrpairs()
      \item man/*: Correct Rd file(s) with Rd \link{} targets missing package anchors
      \item src/sest.cpp: Calloc() and Free() replaced by R_Calloc() and R_Free()
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.7-5 (2024-01-30)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item data/Cars.rda - data set added
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item rrcov-utils.Rd: (mail from Kurt Hornik from 22.01.2024)
      \item fish.Rd: URL corrected
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.7-4 (2023-06-22)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item DESCRIPTION: ISBN changed to <urn:isbn:978-3838108148> 
        (see mail from Kurt Hornik from 06.06.2023)
      \item tests/thubert.Rout.save, tlda.Rout.save, tmcd4.Rout.save, tqda.Rout.save: 
        Follow up to changes in robustbase 0.99-0
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.7-3 (2023-05-11)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item data/ionosphere.rda - data set added
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item R/CovRobust, plot_utils: 'xyqqchi2' and 'xydistance' plots, fixed to use of parameter id.n
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.7-2 (2022-10-22)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed -  warning: invalid UTF-8 in comment [-Winvalid-utf8] in src/rrcov_init.c
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.7-1 (2022-07-18)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item NOTHING NEW YET
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed - maxcond passed to the low level functions
      \item Fixed errors "Found if() conditions comparing class() to string" in detmrcd.R and plot_utils.R
      \item Fixed NOTE "Found the following math rendering problems" in PcaHubert.Rd, see mail from Kurt Hornik from 11 August 2022
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.7-0 (2022-04-25)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Fix for USE_FC_LEN_T becoming the default in 4.2.0
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fortran warnings fixed in ds11.f and fsada.f
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.6-2 (2022-02-08)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Examples for PCA/adjustment for skewed data and PCA/percentage of explained 
        variance added
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item minor differences in tests for PcaProj() on some platforms fixed
      \item scoreplot() corrected to show the labels of the samples
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.6-1 (2022-01-21)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Data set Fruit added: fruit.rda
      \item URLs in Rd files replaced by DOIs to fix for the migration of the 
        www.jstatsoft.org to a new editorial system (see mail from Achim Zeileis from 06.10.2021)
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a problem when showing the percentage of explained variance in 
        summary() of all PCA functions when k is chosen to be less than the number of variables in the 
        input data matrix (k < p). The returned object contains now additionally 
        the rank of the input matrix, the original eigenvalues (of all variables) 
        and the original total variance, if available.
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.6-0 (2021-08-27)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item PcaHubert: option for adjusted outlyingness for skewed data added
      \item Data set Computer Hardware added: machines.rda
      \item Data set Wolves added: wolves.rda
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed some URLs, particularly the reference to javasoft.org
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.5-5 (2020-07-31)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item PcaHubert will crash if X is 1-dimensional and mcd=FALSE: fixed.
      \item Fixed 'noLD' issues in tlda.R and okg4.R
      \item Fixed a problem with wrong scores in PcaProj() (reported by Matthieu Lesnoff <matthieu.lesnoff@gmail.com>)
      \item Fixed a problem with nsamp="exact" or nsamp="best" in CovMve(), CovSest() (reported by Claudio Agostinelli) - 
        these functions, differently from CovMcd, should not take non-numeric 'nsamp'
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.5-2 (2020-01-14)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item  Added parameter 'control' to Linda - to select the robust location and 
        covariance estimator to use in LDA. Now any estimator derived from class 'Cov'
        can be used, even such that are not in 'rrcov'. Return this parameter in the 
        returned S4 object.
      \item Linda returns now the location and covariance estimator used as 'covobj'. This
        is useful for controlling cross-validation, for example.
      \item Linda and LdaClassic use generalized inverse if the common covariance matrix is singular.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an issue the 'predict' function.
      \item Removed the dependence on packages 'cluster' and 'ellipse'.
    }
  }
}
\section{CHANGES in rrcov VERSION 1.4-9 (2019-11-25)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item  Added data set \code{diabetes}; data set \code{soil} from package \code{rrcovHD} moved here.
      \item  Linear and quadratic discriminant analysis can use the MRCD estimates.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an issue with CovControlMcd(nsamp="deterministic") - this would not work, 
        because nsamp was defined in the class definition as "numeric". Now it is "Cnumeric" - union of "character" and "numeric'.
      \item Corrected the code for Minimum Regularized Covariance Determinant estimator (MRCD) - 
        \code{CovMrcd()} - the step of adjusting the eignevalues in r6pack() is excluded 
        now because it has no effect when n > p.
    }
  }
}
\section{CHANGES in rrcov VERSION 1.4-7 (2018-11-14)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added Minimum Regularized Covariance Determinant estimator (MRCD) 
        (Boudt et al. 2018) - 
        functions \code{CovMrcd()} and \code{CovControlMrcd()}
        
      \item  Added data set \code{octane}; data set \code{olitos} from package rrcovHD moved here.
      
      \item The 'pairs' plot is now available for classical covariance matrix
    }
  }
}
