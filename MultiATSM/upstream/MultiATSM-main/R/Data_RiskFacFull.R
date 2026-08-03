#' @title Data: Full set of risk factors - Candelon and Moura (2024, JFEC)
#'
#' @description Full set of risk factors data used in Candelon and Moura (2024, JFEC)
#' @name RiskFacFull
#' @aliases RiskFacFull
#' @docType data
#' @usage data("RiskFacFull")
#' @format matrix containing the full risk factors: (i) global unspanned factors (global economic activity and global inflation); (ii) domestic unspanned factors (economic activity and inflation); and (iii) domestic spanned factors (level, slope, and curvature).
#' Economic system is formed by Brazil, China, Mexico and Uruguay. The data have monthly frequency and span the period from June/2004 to January/2020.
#' @source
#' \describe{
#' \item{Global unspanned factor}{ See \code{data("GlobalMacro")} for a detailed data description.}
#' \item{Domestic unspanned factor}{See \code{data("DomMacro")} for a detailed data description.}
#' \item{Domestic spanned factor}{First three principal components of each country set of bond yields. See \code{data("Yields")} for a detailed data description.}
#' }
#' @references Candelon, B. and Moura, R. (2024) "A Multicountry Model of the Term Structures of Interest Rates with a GVAR". (Journal of Financial Econometrics)
#' @keywords Risk Factors
NULL
