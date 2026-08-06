#'
#' @import RPEIF
#'
#' @title Standard Error Estimate for Rachev Ratio of Returns
#'
#' @description \code{RachevRatio.SE} computes the standard error of the Rachev ratio of the returns.
#'
#' @param data Data of returns for one or multiple assets or portfolios.
#' @param alpha Lower tail probability.
#' @param beta Upper tail probability.
#' @param se.method A character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One or a combination of:
#' \code{"IFiid"} (default), \code{"IFcor"}, \code{"IFcorPW"}, \code{"IFcorAdapt"} (default),
#' \code{"BOOTiid"} or \code{"BOOTcor"}.
#' @param cleanOutliers Boolean variable to indicate whether the pre-whitenning of the influence functions TS should be done through
#' a robust filter. Default if FALSE.
#' @param fitting.method Distribution used in the standard errors computation. Should be one of "Exponential" (default) or "Gamma".
#' @param d.GLM.EN Order of the polynomial for the Exponential or Gamma fitting. Default polynomial order of 5.
#' @param freq.include Frequency domain inclusion criteria. Must be one of "All" (default), "Decimate" or "Truncate."
#' @param freq.par Percentage of the frequency used if \code{"freq.include"} is "Decimate" or "Truncate." Default is 0.5.
#' @param corOut Return correlation of the returns or the influence function transformed returns.
#' Must be one of "retCor", "retIFCor" or "none" (default).
#' @param return.coef Boolean variable to indicate whether the coefficients of the penalized GLM fit are returned. Default if FALSE.
#' @param ... Additional parameters.
#'
#' @return A vector or a list depending on \code{se.method}.
#'
#' @export
#'
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
#'   # the two influence functions based approaches
#'   RachevRatio.SE(edhec, se.method = c("IFiid","IFcorAdapt"),
#'                  cleanOutliers = FALSE,
#'                  fitting.method = c("Exponential", "Gamma")[1])
#' }
#'
RachevRatio.SE <- function(data, alpha = 0.1, beta = 0.1,
                           se.method = c("IFiid","IFcor","IFcorAdapt","IFcorPW","BOOTiid","BOOTcor")[c(1,4)],
                           cleanOutliers = FALSE, fitting.method = c("Exponential", "Gamma")[1], d.GLM.EN  =  5,
                           freq.include = c("All", "Decimate", "Truncate")[1], freq.par = 0.5,
                           corOut  =  c("none", "retCor","retIFCor", "retIFCorPW")[1],
                           return.coef = FALSE,
                           ...){

  # Point estimate
  if(is.null(dim(data)) || ncol(data) == 1)
    point.est <- RachevRatio(data, alpha = alpha, beta = beta) else
      point.est <- apply(data, 2, function(x) RachevRatio(x, alpha = alpha, beta = beta))

    # SE Computation
    if(is.null(se.method)){
      return(point.est)
    } else{
      SE.out <- list(RachevRatio = point.est)
      for(mymethod in se.method){
        SE.out[[mymethod]] <- EstimatorSE(data, estimator.fun = "RachevRatio", alpha = alpha, beta = beta,
                                          se.method = mymethod,
                                          cleanOutliers = cleanOutliers,
                                          fitting.method = fitting.method, d.GLM.EN = d.GLM.EN,
                                          freq.include = freq.include, freq.par = freq.par,
                                          return.coef = return.coef,
                                          ...)
      }

      # Adding the correlations to the list
      SE.out <- Add_Correlations(SE.out = SE.out, data = data, cleanOutliers = cleanOutliers, corOut = corOut, IF.func = IF.RachevRatio, ...)

      # Returning the output
      return(SE.out)
    }
}
