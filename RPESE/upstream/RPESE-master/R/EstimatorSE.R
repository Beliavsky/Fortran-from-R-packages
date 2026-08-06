#'
#' @import RPEIF
#' @importFrom stats cov quantile SSD arima lm coef ar na.omit fft approxfun density
#' @importFrom boot boot tsboot
#' @importFrom xts is.xts
#' @importFrom zoo zoo
#'
#' @title Wrapper Function for Standard Errors Estimates Functions
#'
#' @description \code{EstimatorSE} computes the standard error for specified risk and performance measures.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param estimator.fun Risk or performance measure to compute estimates of standard errors.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of:
#' \code{"IFiid"}, \code{"IFcor"}, \code{"IFcorAdapt"}, \code{"IFcorPW"},
#' \code{"BOOTiid"}, \code{"BOOTcor"}, or \code{"none"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through
#' a robust filter. Default if FALSE.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
#' @param d.GLM.EN Order of the polynomial for the Exponential or Gamma fitting. Default polynomial order of 5.
#' @param freq.include Frequency domain inclusion criteria. Must be one of "All" (default), "Decimate" or "Truncate."
#' @param freq.par Percentage of the frequency used if \code{"freq.include"} is "Decimate" or "Truncate." Default is 0.5.
#' @param a First adaptive method parameter.
#' @param b Second adaptive method parameter.
#' @param return.coef Boolean variable to indicate whether the coefficients of the Exponential or Gamma fit are returned. Default is FALSE.
#' @param ... Additional parameters.
#'
#' @return A vector standard error estimates.
#'
#' @export
#'
#' @author Xin Chen, \email{chenx26@uw.edu}
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @examples
#' # Examples using PerformanceAnalytics data
#' if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
#'   # Loading data
#'   data(edhec, package = "PerformanceAnalytics")
#'   # Changing the data colnames
#'   names(edhec) <- c("CA", "CTA", "DIS", "EM", "EMN",
#'                     "ED", "FIA", "GM", "LS", "MA",
#'                     "RV", "SS", "FOF")
#'   # Computing the standard errors for
#'   # the three influence functions based approaches
#'   EstimatorSE(edhec[,"CA"], se.method = c("IFcor"),
#'               cleanOutliers = FALSE,
#'               fitting.method = c("Exponential", "Gamma")[1])
#' }
#'
EstimatorSE <- function(data,
                        estimator.fun = c("DSR", "ES","ESratio","LPM",
                                          "Mean","OmegaRatio","RachevRatio","robMean",
                                          "SD","SemiSD","SR","SoR",
                                          "VaR","VaRratio")[1],
                        se.method = c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[1],
                        cleanOutliers = FALSE,
                        fitting.method = c("Exponential", "Gamma")[1], d.GLM.EN = 5,
                        freq.include = c("All", "Decimate", "Truncate")[1], freq.par = 0.5,
                        a = 0.3, b = 0.7,
                        return.coef = FALSE,
                        ...){

  # Available estimator functions
  estimators.available <- c("DSR", "ES","ESratio","LPM",
                            "Mean","OmegaRatio","RachevRatio","robMean",
                            "SD","SemiSD","SR","SoR",
                            "VaR","VaRratio")

  # Checking if the specified risk measure is available
  if(!(estimator.fun %in% estimators.available))
    stop("The specified estimator function is not available.")

  # Available SE methods
  se.available <- c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")
  # Checking if the standard error method is available
  if(!(se.method %in% se.available))
    stop("The specified standard error method is not available.")

  # Available fitting methods
  fitting.available <- c("Exponential", "Gamma")
  # Checking if the standard error method is available
  if(!(fitting.method %in% fitting.available))
    stop("The specified fitting method is not available.")

  # Available frequency inclusion
  frequency.available <- c("All", "Decimate", "Truncate")
  # Checking if the standard error method is available
  if(!(freq.include %in% frequency.available))
    stop("The specified frequency inclusion criteria is not available.")
  # Checking frequency parameter
  if(any(!is.numeric(freq.par), length(freq.par) != 1, freq.par<0, freq.par>1))
    stop("The frequency parameter must be a numerical value between 0 and 1.")

  # Storing RPE point estimate function
  RPE_fun <- switch(estimator.fun,
                    DSR = DSR,
                    ES = ES,
                    ESratio = ESratio,
                    LPM = LPM,
                    Mean = Mean,
                    OmegaRatio = OmegaRatio,
                    SD = SD,
                    SemiSD = SemiSD,
                    SoR = SoR,
                    SR = SR,
                    RachevRatio = RachevRatio,
                    robMean = robMean,
                    VaR = VaR,
                    VaRratio = VaRratio,
                    stop("The estimator.fun specified is not implemented yet, please contact Anthony Christidis (anthony.christidis@stat.ubc.ca) or submit an issue at the github repository.")
  )

  # Storing RPE IF function
  RPE.IF_fun <- switch(estimator.fun,
                       DSR = IF.DSR,
                       ES = IF.ES,
                       ESratio = IF.ESratio,
                       LPM = IF.LPM,
                       Mean = IF.Mean,
                       OmegaRatio = IF.OmegaRatio,
                       RachevRatio = IF.RachevRatio,
                       robMean = IF.robMean,
                       SD = IF.SD,
                       SemiSD = IF.SemiSD,
                       SoR = IF.SoR,
                       SR = IF.SR,
                       VaR = IF.VaR,
                       VaRratio = IF.VaRratio,
                       stop("The estimator.fun specified is not implemented yet, please contact Anthony Christidis (anthony.christidis@stat.ubc.ca) or submit an issue at the github repository.")
  )

  SE.out <- switch(se.method,
    none = NULL,
    IFiid = SE.xts(data, SE.IF.iid, RPE_fun, RPE.IF_fun),
    IFcor = SE.xts(data, SE.IF.cor, RPE_fun, RPE.IF_fun,
                   prewhiten = FALSE, cleanOutliers = cleanOutliers, fitting.method = fitting.method, d.GLM.EN = d.GLM.EN,
                   freq.include = freq.include, freq.par = freq.par,
                   return.coef = return.coef,
                   ...),
    IFcorPW = SE.xts(data, SE.IF.cor, RPE_fun, RPE.IF_fun,
                     prewhiten = TRUE, cleanOutliers = cleanOutliers, fitting.method = fitting.method, d.GLM.EN = d.GLM.EN,
                     freq.include = freq.include, freq.par = freq.par,
                     return.coef = return.coef,
                     ...),
    IFcorAdapt = list(cor = SE.xts(data, SE.IF.cor, RPE_fun, RPE.IF_fun, prewhiten = FALSE,
                                 cleanOutliers = cleanOutliers, fitting.method = fitting.method, d.GLM.EN = d.GLM.EN,
                                 freq.include = freq.include, freq.par = freq.par,
                                 return.coef = return.coef,
                                 ...),
                      corPW = SE.xts(data, SE.IF.cor, RPE_fun, RPE.IF_fun, prewhiten = TRUE,
                                   cleanOutliers = cleanOutliers, fitting.method = fitting.method, d.GLM.EN = d.GLM.EN,
                                   freq.include = freq.include, freq.par = freq.par,
                                   return.coef = return.coef,
                                   ...)),
    BOOTiid = SE.xts(data, SE.BOOT.iid, RPE_fun, ...),
    BOOTcor = SE.xts(data, SE.BOOT.cor, RPE_fun, ...)
  )

  # Formatting and returning output results
  SE.out <- EstimatorSE.format(SE.out = SE.out, data = data, se.method = se.method, a = a, b = b, return.coef = return.coef)
  return(SE.out)
}










