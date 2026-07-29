#' smoots: A package for data-driven nonparametric estimation of the trend and
#' its derivatives in equidistant time series.
#'
#' The \code{smoots} package provides different applicable functions for the
#' estimation of the trend or its derivatives in equidistant time series.
#' The main functions include an automated bandwidth selection method for time
#' series with short-memory errors. With package version 1.1.0 several
#' functions for forecasting as well as linearity tests were added.
#'
#'@section Functions (version 1.0.0):
#'The \code{smoots} functions are either meant for calculating nonparametric
#'estimates of the trend of a time series or its derivatives.
#'
#'\code{msmooth} is the central function of the package. It allows
#'the user to conduct a local polynomial regression of the trend based on
#'an optimal bandwidth that is obtained by an iterative plug-in algorithm.
#'There are also different algorithms implemented concerning the inflation rate
#'and other factors that can be chosen from (see also: \code{\link{msmooth}}).
#'
#'\code{dsmooth} is a function that calculates the derivatives of the
#'trend after obtaining the optimal bandwidth by an iterative plug-in
#'algorithm. The estimates are obtained for rescaled time points on the
#'interval \eqn{[0, 1]} (see also: \code{\link{dsmooth}}).
#'
#'\code{tsmooth} is similar to \code{msmooth} as it also calculates the trend
#'of the series. Instead of using the name of a predefined algorithm that
#'settles the inflation rate (and other factors), these factors can be manually
#'and individually adjusted as arguments in the function (see also:
#'\code{\link{tsmooth}}).
#'
#'\code{gsmooth} is a standard smoothing function that applies the local
#'polynomial regression method. A bandwidth can be chosen freely. The estimates
#'are obtained for rescaled time points on the interval \eqn{[0, 1]}
#'(see also: \code{\link{gsmooth}}).
#'
#'\code{knsmooth} is a standard smoothing function that applies the kernel
#'regression method. A bandwidth can be chosen freely
#'(see also: \code{\link{knsmooth}}).
#'
#'@section Added Functions (version 1.1.0):
#'With the publication of the package version 1.1.0, new functions were added.
#'These include functions for forecasting and functions for testing linearity
#'of the deterministic trend.
#'
#'\code{rescale} helps rescaling the estimates of the derivatives of the
#'nonparametric trend function, because the estimates are obtained for rescaled
#'time points on the interval \eqn{[0, 1]} and not for the actual time points
#'(see also: \code{\link{rescale}}).
#'
#'\code{critMatrix} is a quick tool for the calculation of information criteria
#'for ARMA(\eqn{p,q}) models with different order combinations \eqn{p} and
#'\eqn{q}. The function returns a matrix with the obtained values of the
#'selected criterion for the different combinations of \eqn{p} and \eqn{q}
#'(see also: \code{\link{critMatrix}}).
#'
#'\code{optOrd} is useful in combination with \code{critMatrix}. It
#'reads a matrix equal in structure to the ones returned by \code{critMatrix}
#'and returns the optimal orders p and q. Furthermore, additional restrictions
#'for the selection can be imposed (see also: \code{\link{optOrd}}).
#'
#'\code{normCast} provides means to obtain point forecasts as well as
#'forecasting intervals for a given series under the assumption that it follows
#'an ARMA(\eqn{p,q}) model and that its innovations are normally distributed
#'(see also: \code{\link{normCast}}).
#'
#'\code{bootCast} can also be used to calculate point forecasts and forecasting
#'intervals, if the series of interest follows an ARMA(\eqn{p,q}) model.
#'However, the main difference is that this function should be applied, if the
#'distribution of the innovations is unknown or explicitly non-Gaussian,
#'because the underlying bootstrap process does not need a distribution
#'assumption. In spite of this advantage, it also needs significantly more
#'computation time. Therefore, with version 1.1.1, the bootstrap can now be
#'conducted in parallel for an improved performance (see also:
#'\code{\link{bootCast}}).
#'
#'\code{trendCast} uses a \code{smoots} object that is the output of a trend
#'estimation and calculates point forecasts of the trend. Forecasting intervals
#'are omitted for reasons that are explained within the function's
#'documentation (see also: \code{\link{trendCast}}).
#'
#'\code{modelCast} calculates the point forecasts and forecasting intervals of
#'a trend-stationary series. Based on a \code{smoots} object that is the output
#'of a trend estimation, \code{trendCast} is applied to the estimated trend,
#'whereas either \code{normCast} or \code{bootCast} is applied to the
#'residuals (see also: \code{\link{modelCast}}).
#'
#'\code{rollCast} is a backtesting function for Semi-ARMA models. A given
#'series is divided into in-sample and out-of-sample observations, where the
#'in-sample is used to fit a Semi-ARMA model. One-step rolling forecasts and
#'the corresponding forecasting intervals are then obtained for the
#'out-of-sample observations. The quality of the model is then assessed via
#'a comparison of the true out-of-sample observations with the point forecasts
#'and forecasting intervals. Different quality criteria are calculated and
#'returned (see also: \code{\link{modelCast}}).
#'
#'\code{confBounds} calculates the confidence bounds of the estimated trend
#'or of the estimated derivative of the trend at a predefined confidence
#'level. A graphic of the confidence bounds alongside a previously chosen
#'constant or parametric illustration of the estimated series is created. With
#'this plot it can be tested graphically if the deterministic trend is constant
#'or if it follows a parametric polynomial model. Furthermore, it can also be
#'tested, if the derivatives of the trend are constant (see also:
#'\code{\link{confBounds}}).
#'
#'@section Datasets:
#'The package includes four datasets: \code{gdpUS} (see also:
#'\code{\link{gdpUS}}) that has data concerning the quarterly US GDP between
#'Q1 1947-01 and Q2 2019, \code{tempNH} (see also: \code{\link{tempNH}}) with
#'mean monthly Northern Hemisphere temperature changes from 1880 to 2018,
#'\code{dax} (see also: \code{\link{dax}}) with daily financial time series
#'data of the German stock index (DAX) from 1990 to July 2019 and \code{vix}
#'(see also: \code{\link{vix}}) with daily financial time series data of the
#'CBOE Volatility Index (VIX) from 1990 to July 2019.
#'
#'@section License:
#'The package is distributed under the General Public License v3
#'([GPL-3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))).
#'
#'@references
#' Feng, Y., Gries, T. and Fritz, M. (2020). Data-driven
#' local polynomial for the trend and its derivatives in economic time
#' series. Journal of Nonparametric Statistics, 32:2, 510-533.
#'
#' Feng, Y., Gries, T., Letmathe, S. and Schulz, D. (2019). The smoots package
#' in R for semiparametric modeling of trend stationary time series. Discussion
#' Paper. Paderborn University. Unpublished.
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}
#'
#' @docType package
#' @name smoots
#' @aliases smoots-package
#' @useDynLib smoots
#' @importFrom Rcpp sourceCpp
NULL
