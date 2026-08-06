#' @title Influence Function - Expected Shortfall (ES)
#' 
#' @description \code{IF.ES} returns the data and plots the shape of either the IF or the IF TS for the ES
#'
#' @param returns Returns data of the asset or portfolio. This can be a numeric or an xts object.
#' @param evalShape Evaluation of the shape of the IF risk or performance measure if TRUE. Otherwise, a TS of the IF of the provided returns is computed.
#' @param retVals Values used to evaluate the shape of the IF.
#' @param nuisPars Nuisance parameters used for the evaluation of the shape of the IF (if no returns are provided).
#' @param k Range parameter for the shape of the IF (the SD gets multiplied k times).
#' @param IFplot If TRUE, the plot of the IF shape or IF TS of the returns is produced.
#' @param IFprint If TRUE, the data for the IF shape or the IF TS of the returns is returned.
#' @param alpha.ES Tail Probability.
#' @param prewhiten Boolean variable to indicate if the IF TS is pre-whitened (TRUE) or not (FALSE).
#' @param ar.prewhiten.order Order of AR parameter for the pre-whitening. Default is AR(1).
#' @param cleanOutliers Boolean variable to indicate whether outliers are cleaned with a robust location and scale estimator.
#' @param cleanMethod Robust method used to clean outliers from the TS. Default choice is "locScaleRob". 
#' @param eff Tuning parameter for the normal distribution efficiency for the "locScaleRob" robust data cleaning.
#' @param ... Additional parameters.
#'
#' @return Influence function of the ES.
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
#' outIF <- IF.ES(returns = NULL, evalShape = TRUE, 
#'                retVals = NULL, nuisPars = NULL,
#'                IFplot = TRUE, IFprint = TRUE)
#'
#' # Examples using PerformanceAnalytics data
#' if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
#'   data(edhec, package = "PerformanceAnalytics")
#'   colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
#'                       "GM", "LS", "MA", "RV", "SS", "FoF") 
#'   
#'   # Plot of IF a specified TS 
#'   outIF <- IF.ES(returns = edhec[,"CA"], evalShape = TRUE, 
#'                  retVals = seq(-0.1, 0.1, by = 0.001), nuisPars = NULL,
#'                  IFplot = TRUE, IFprint = TRUE)
#'   
#'   # Computing the IF of the returns (with prewhitening) with a plot of IF TS
#'   outIF <- IF.ES(returns = edhec[,"CA"], evalShape = FALSE, 
#'                  retVals = NULL, nuisPars = NULL,
#'                  IFplot = TRUE, IFprint = TRUE,
#'                  prewhiten = FALSE)
#' }
#'
IF.ES <- function(returns = NULL, evalShape = FALSE, retVals = NULL, nuisPars = NULL, k = 4,
                  IFplot = FALSE, IFprint = TRUE,
                  alpha.ES = 0.05, prewhiten = FALSE, ar.prewhiten.order = 1,
                  cleanOutliers = FALSE, cleanMethod = c("locScaleRob")[1], eff = 0.99, 
                  ...){
  
  # Checking input data
  DataCheck(returns = returns, evalShape = evalShape, retVals = retVals, nuisPars = nuisPars, k = k,
            IFplot = IFplot, IFprint = IFprint,
            prewhiten = prewhiten, ar.prewhiten.order = ar.prewhiten.order,
            cleanOutliers = cleanOutliers, cleanMethod = cleanMethod, eff = eff)
  
  # Checking input for alpha.ES
  if(!inherits(alpha.ES, "numeric")){
    stop("alpha.ES should be numeric")
  } else if(any(alpha.ES < 0, alpha.ES > 1)) {
    stop("alpha.ES should be a numeric value between 0 and 1.")
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
    IFvals <- EvaluateShape(estimator = "ES",
                            retVals = retVals, returns = returns, k = k, nuisPars = nuisPars,
                            IFplot = IFplot, IFprint = IFprint, alpha.ES = alpha.ES)
    if(IFprint)
      return(IFvals) else{
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt)) 
        stop() 
      }
  }
  
  # Storing the dates
  if(xts::is.xts(returns))
    returns.dates <- zoo::index(returns)
  
  # IF Computation
  IF.ES.vector <- IF.ES.fn(x = returns, returns = returns, alpha.ES = alpha.ES)
  
  # Adding the pre-whitening functionality  
  if(prewhiten)
    IF.ES.vector <- as.numeric(arima(x = IF.ES.vector, order = c(ar.prewhiten.order,0,0), include.mean = TRUE)$residuals)
  
  # Adjustment for data (xts)
  if(xts::is.xts(returns))
    IF.ES.vector <- xts::xts(IF.ES.vector, returns.dates)
  
  # Plot of the IF TS
  if(isTRUE(IFplot)){
    print(plot(IF.ES.vector, type = "l", main = "ES Estimator Influence Function Transformed Returns", ylab = "IF"))
  }
  
  # Stop if no printing of the TS
  if(!IFprint){
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt)) 
    stop() 
  }
  
  # Returning the IF vector for the ES
  if(xts::is.xts(returns))
    return(xts::xts(IF.ES.vector, returns.dates)) else
      return(IF.ES.vector)
}