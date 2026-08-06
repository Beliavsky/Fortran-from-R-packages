#' @title Formatted Output for Standard Errors Functions in RPESE
#'
#' @description \code{printSE} returns a formatted output from standard error functions from RPESE.
#'
#' @param SE.data Standard error estimates output from RPESE functions.
#' @param round.digit Number of digits for rounding.
#' @param round.out Round data (TRUE) with round.digit number of digits. Default is TRUE.
#'
#' @return A data frame with formatted output from standard error functions from \code{RPESE}.
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
#'   # the two influence functions based approaches
#'   ES.out <- ES.SE(edhec, se.method = c("IFiid","IFcor"),
#'                   cleanOutliers = FALSE,
#'                   fitting.method = c("Exponential", "Gamma")[1])
#'   # Print the output
#'   printSE(ES.out)
#' }
#'
printSE <- function(SE.data, round.digit = 3, round.out = TRUE){

  # Printing formatting
  list.names <- names(SE.data)
  row
  SE.data.df <- as.numeric(data.frame(t(SE.data[[1]])))
  for(i in 2:length(list.names)){
    if(list.names[i] %in% c("IFiid","IFcor", "IFcorAdapt", "IFcorPW", "BOOTiid", "BOOTcor"))
      SE.data.df <- cbind(SE.data.df, SE.data[[i]]$se) else
        if(list.names[i] %in% c("none", "retCor","retIFCor", "retIFCorPW"))
          SE.data.df <- cbind(SE.data.df, SE.data[[i]]$out)
  }
  colnames(SE.data.df) <- list.names
  rownames(SE.data.df) <- names(SE.data[[1]])

  # Rounding of the output of value returned
  if(round.out){
    return(round(SE.data.df, digits = round.digit))
  } else
    return(SE.data.df)
}
