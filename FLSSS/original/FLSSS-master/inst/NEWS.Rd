\name{NEWS}
\alias{NEWS}
\title{FLSSS updates}


\section{Changes in FLSSS version 9.0.0 (20220415)}{
  \itemize{
    \item Implemented multidimensional exact subset sum in arbitrary precision and magnitude, \code{arbFLSSS()}, and its distributed computing counterpart, \code{decomposeArbflsss()} and \code{arbFLSSSobjRun()}. A simple Bloom filter is added as an option to accelerate the algorithms.
  }
}


\section{Changes in FLSSS version 8.6.0 (20200905)}{
  \itemize{
    \item Implemented the distributed computing infrastructure, \code{decomposeMflsss()} and \code{mFLSSSobjRun()}, for multidimensional subset sum.
    \item Changed the timer from measuring total CPU time to measuring the elapsed wall time for multithreaded functions.
  }
}


\section{Changes in FLSSS version 8.5.6 (20191027)}{
  \itemize{
    \item Bug fix. A typo in the code related to integer type selection resulted in crash of the multidimensional subset sum miner given superset size > 127.
  }
}


\section{Changes in FLSSS version 8.5.5 (20190709)}{
  \itemize{
    \item Bug reported and fixed for \code{FLSSS::FLSSS()}. At a certain step, the algorithm uses a relative error threshold to decide if two floats are equal. Due to limited numeric precision, in an unusual instance the algorithm ended up comparing two subset sums 1e-14 and 1.2e-14. These two sums should both equal 0 theoretically (superset has negative elements). A wrong decision was thus made upon the significant relative error between 1e-14 and 1.2e-14. Bug fix consists of shifting superset \code{v} to nonnegative domain and scaling error bound \code{ME} to 1 before mining. The rest functions are not affected by this incident.
    \item Added namespaces to suppress link-time optimization warnings.
  }
}


\section{Changes in FLSSS version 8.5 (20190106)}{
  \itemize{
    \item Bug fix. The previous branch-and-bound knapsack solver did not initialize containers with sufficient sizes. In rare cases a pointer would be realized at the end of the container.
    \item Added a genetic algorithm with local heuristics for generalized assignment problem.
  }
}


\section{Changes in FLSSS version 8.3 (20190106)}{
  \itemize{
    \item Implemented different algorithms for the binary knapsack and the generalized assignment problems. Functions prefixed with 'aux' (auxiliary) are new.
    \item Future updates will dedicate to implementing existing foundational or state-of-the-art algorithms for combinatorial problems of interest. The implementation will focus on computing speed as usual.
  }
}


\section{Changes in FLSSS version 7.7 (20181121)}{
  \itemize{
    \item Integerization for the generalized assignment solver is canceled. A new compact representation for cost matrix is implemented and accelerates mining speed.
    \item Bug fix to the output format of \code{GAP()} given no solution exists.
    \item Package vignette is pushed to arXiv.
  }
}

