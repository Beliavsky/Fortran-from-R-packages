\name{parma-package}
\alias{parma-package}
\alias{parma}
\title{The parma package}
\description{
Portfolio Allocation and Risk Management. Models and Methods for scenario and 
moment based optimization of portfolios.
}
\details{
\tabular{ll}{
Package: \tab parma\cr
Type: \tab Package\cr
Version: \tab 1.5-2\cr
Date: \tab 2014-07-09\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
Imports: \tab nloptr, Rglpk, quadprog\cr
Suggests: \tab Rsymphony, truncnorm, timeSeries\cr}

The portfolio allocation and risk managament applications (parma) package 
contains a unique set of methods and models for the optimal allocation of 
capital in financial portfolios. It uniquely represents certain discontinuous 
problems using their smooth approximation counterparts and implements fractional 
based programming for the direct optimization of risk-to-reward ratios. In 
combination with the rmgarch package, it enables the confident solution to 
scenario based optimization problems using such risk and deviation measures as 
Mean Absolute Deviation (MAD), Variance (EV), Minimax, Conditional Value at 
Risk (CVaR), Conditional Drawdown at Risk (CDaR) and Lower Partial Moments (LPM). 
In addition, it implements moment based optimization for use with the quadratic 
EV problem, and a higher moment CARA utility expansion using the coskewness and 
cokurtosis matrices generated from the GO-GARCH with affine GH or NIG 
distributions. Benchmark relative optimization (tracking error) is also 
implemented as are basic mixed integer cardinality constraints. Finally, for 
non-convex problem formulations such as the upper to lower partial moments 
function, global optimization methods using a penalty based method are 
available. The key functions in the package are \code{\link{parmaspec}} 
which defines the optimization setup, and \code{\link{parmasolve}} which solves
the problem given a chosen representation and solver. A portfolio frontier 
function is implemented in \code{\link{parmafrontier}}, utility optimization in
\code{\link{parmautility}} and a custom translation of the cmaes global 
optimization solver of Hansen (2006) with full features is implemented in 
\code{\link{cmaes}}.
}
\section{How to cite this package}{
  Whenever using this package, please cite as\cr

 \preformatted{@Manual{Galanos_2014,
 author       = {Alexios Galanos and Bernhard Pfaff},
 title        = {{parma}: Portfolio Allocation and Risk Management Applications.},
 year         = {2014},
 note 	      = {R package version 1.5-1.},}}
}
\section{License}{
  The releases of this package is licensed under GPL version 3.
}
\author{Alexios Galanos and Bernhard Pfaff}
\references{
Charnes, A. and Cooper, W. 1962, Programming with linear fractional functionals,
\emph{Naval Research Logistics Quarterly}, \bold{9}, 181--186.\cr
Dinkelbach, W. 1967, On nonlinear fractional programming, \emph{Management Science},
\bold{13(7)}, 492--498.\cr
Fishburn, P.C. 1977, Mean-risk analysis with risk associated with below-target 
returns, \emph{The American Economic Review}, \bold{67(2)}, 116-126.\cr
Galanos, A. 2012, Higher Moment Models for Risk and Portfolio Management, 
Thesis (submitted) \emph{Cass Business School}.\cr
Hansen, N. 2006, The CMA Evolution Strategy: A Comparing Review, \emph{Towards 
a New Evolutionary Computation (Studies in Fuzziness and Soft Computing)}, 
\bold{192}, 75--102.\cr
Holthausen, D. 1981, A risk-return model with risk and return measured as 
deviations from a target return, \emph{The American Economic Review}, \bold{71},
182--188.\cr
Konno, H. and Yamazaki, H. 1991, Mean-absolute deviation portfolio optimization 
model and its applications to Tokyo stock market, \emph{Management Science},
\bold{37(5)}, 519--531.\cr
Markowitz, H. 1952, Portfolio selection, \emph{The Journal of Finance}, \bold{7(1)},
77--91.\cr
Rockafellar, R.T. and Uryasev, S. and Zabarankin, M., 2006, Generalized deviations 
in risk analysis, \emph{Finance and Stochastics}, \bold{10(1)}, 51--74.\cr
Stoyanov, S.V. and Rachev, S.T. and Fabozzi, F.J. 2007, Optimal financial portfolios,
\emph{Applied Mathematical Finance}, \bold{14(5)}, 401--436.\cr
}
