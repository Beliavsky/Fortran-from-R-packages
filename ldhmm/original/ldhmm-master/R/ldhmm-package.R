#' ldhmm: A package for HMM using lambda distribution.
#'
#' The ldhmm package provides the core class and functions to calculate
#' Hidden Markov Model (HMM) using lambda distribution framework.
#' The main goal is to provide a theoretically solid foundation to explore 
#' the return time-series in the financial market, where the normal distribution 
#' is not adequate due to the leptokurtic nature of the data. 
#' Major features in the S&P 500 index, such as regime identification, volatility clustering, 
#' and anti-correlation between return and volatility, can be extracted from HMM cleanly. 
#' Univariate symmetric lambda distribution is essentially a location-scale family 
#' of power-exponential distribution. Such distribution is suitable for describing 
#' highly leptokurtic time series obtained from the financial market. 
#' 
#' The main change compared to a normal-distribution based HMM is to add the third paramter
#' \code{lambda} to describe the kurtosis level of the distribution. When \code{lambda} is one, the model
#' converges back to a normal-distribution based HMM (e.g. using depmixS4 package).
#' The ability to optimize kurtosis brings the model output to be more consistent with the data.
#' In particular, for daily data, the level of kurtosis is quite high. This puts
#' the normal distribution in great disadvantage. This problem is solved by using the
#' lambda distribution. 
#'
#' @author Stephen H-T. Lihn
#'
#' @docType package
#' @name ldhmm-package
#' @import xts methods graphics moments parallel
#'
#' @references Walter Zucchini, Iain L. MacDonald, Roland Langrock (2016). 
#'             "Hidden Markov Models for Time Series, An Introduction Using R."
#'             Second Edition. CRC Press.
#' 
NULL

# Some areas of this package require multi-core capability
cores <- switch( Sys.info()[['sysname']],
    Windows = 1,
    Linux   = parallel::detectCores(),
    Darwin  = parallel::detectCores(),
    parallel::detectCores()
    )

if (is.null(getOption("mc.cores"))) {
    options("mc.cores"=cores)
}

# end

