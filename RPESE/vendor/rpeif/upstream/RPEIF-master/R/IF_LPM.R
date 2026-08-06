#' @title Influence Function - Lower Partial Moment (LPM)
#' 
#' @description \code{IF.LPM} returns the data and plots the shape of either the IF or the IF TS for the LPM
#'
#' @param returns Returns data of the asset or portfolio. This can be a numeric or an xts object.
#' @param evalShape Evaluation of the shape of the IF risk or performance measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param nuisPars Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param const Constant threshold.
#' @param order Order of LPM. Can only take values 1 or 2.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether outliers are cleaned with a robust location and scale estimator.
#' @param cleanMethod Robust method used to clean outliers from the TS. Default choice is "locScaleRob". 
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of LPM.
#' 
#' @details 
#' For further details on the usage of the \code{nuisPars} argument, please refer to Section 3.1 for the \code{RPEIF} vignette.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Plot of IF with nuisance parameter with return value
#' # This example works without any additional packages
#' outIF <- IF.LPM(returns = NULL, evalShape = TRUE, 
#'                 retVals = NULL, nuisPars = NULL,
#'                 IFplot = TRUE, IFprint = TRUE)
#'
#' # Examples using PerformanceAnalytics data
#' if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
#'   data(edhec, package = "PerformanceAnalytics")
#'   colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                       "GM", "LS", "MA", "RV", "SS", "FoF") 
#'   
#'   # Plot of IF a specified TS 
#'   outIF <- IF.LPM(returns = edhec[,"CA"], evalShape = TRUE, 
#'                   retVals = seq(-0.1, 0.1, by = 0.001), nuisPars = NULL,
#'                   IFplot = TRUE, IFprint = TRUE)
#'   
#'   # Computing the IF of the returns (with prewhitening) with a plot of IF TS
#'   outIF <- IF.LPM(returns = edhec[,"CA"], evalShape = FALSE, 
#'                   retVals = NULL, nuisPars = NULL,
#'                   IFplot = TRUE, IFprint = TRUE,
#'                   prewhiten = FALSE)
#' }
#'
IF.LPM <- function(returns = NULL, evalShape = FALSE, retVals = NULL, nuisPars = NULL, k = 4,
                   IFplot = FALSE, IFprint = TRUE,
                   const = 0, order = 1, prewhiten = FALSE, ar.prewhiten.order = 1,
                   cleanOutliers = FALSE, cleanMethod = c("locScaleRob")[1], eff = 0.99, 
                   ...){
  
  # Checking input data
  DataCheck(returns = returns, evalShape = evalShape, retVals = retVals, nuisPars = nuisPars, k = k,
            IFplot = IFplot, IFprint = IFprint,
            prewhiten = prewhiten, ar.prewhiten.order = ar.prewhiten.order,
            cleanOutliers = cleanOutliers, cleanMethod = cleanMethod, eff = eff)
  
  # Checking input for const
  if(!inherits(const, "numeric"))
    stop("const should be numeric")
  
  # Checking input for order
  if(!inherits(order, "numeric")){
    stop("order should be numeric")
  } else if(!(order %in% 1:2)) {
    stop("order should take be either 1 or 2.")
  }

  # Evaluation of nuisance parameters
  nuisPars <- NuisanceData(nuisPars)
  
  # Storing the dates
  if(xts::is.xts(returns))
    returns.dates <- zoo::index(returns)
  
  # Adding the robust filtering functionality
  if(cleanOutliers){
    temp.returns <- robust.cleaning(returns, cleanMethod, eff)
    if(xts::is.xts(returns))
      returns <- xts::xts(temp.returns, returns.dates) else
        returns <- temp.returns
  }
  
  # Plot for shape evaluation
  if(evalShape){
    IFvals <- EvaluateShape(estimator = "LPM",
                            retVals = retVals, returns = returns, k = k, nuisPars = nuisPars,
                            IFplot = IFplot, IFprint = IFprint, const = const, order = order)
    if(IFprint)
      return(IFvals) else{
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt)) 
        stop() 
      }
  }

  # Computing the IF vector for LPM (order = 1 case)
  IF.LPM.vector <- IF.LPM.fn(x = returns, returns = returns, const = const, order = order)
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.LPM.vector <- as.numeric(arima(x = IF.LPM.vector, order = c(ar.prewhiten.order,0,0), include.mean = TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.LPM.vector <- xts::xts(IF.LPM.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.LPM.vector, type = "l", main = "LPM Estimator Influence Function Transformed Returns"))
  }
  
  # Returning the IF vector for LPM (order = 1 case)
  if(xts::is.xts(returns))
    return(xts::xts(IF.LPM.vector, returns.dates)) else
      return(IF.LPM.vector)

}

# ==================
# Function for LPM
# ==================

LPM <- function(returns, const  =  0, order  =  1, ...){
  
  # Computing the LPM
  return(1/length(returns)*sum((const-returns[returns <= const])^order))
}


