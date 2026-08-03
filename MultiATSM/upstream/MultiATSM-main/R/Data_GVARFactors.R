#' @title Data: Risk Factors for the GVAR - Candelon and Moura (2024, JFEC)
#'
#' @description Risk factors data used in the GVAR-ATSM from  Candelon and Moura (2024, JFEC)
#' @name GVARFactors
#' @aliases GVARFactors
#' @docType data
#' @usage data("GVARFactors")
#' @format List of risk factors organized for GVAR estimation. It includes global unspanned factors (economic activity, inflation) and domestic factors—both unspanned (economic activity, inflation) and spanned (level, slope, curvature) with their starred counterparts.
#' The dataset covers Brazil, China, Mexico, and Uruguay at a monthly frequency from June 2004 to January 2020.
#' @source
#' \describe{
#' \item{Global unspanned factors}{ See \code{data("GlobalMacro")} for a detailed data description.}
#' \item{Domestic unspanned factors}{See \code{data("DomMacro")} for a detailed data description.}
#' \item{Domestic spanned factors}{First three principal components of each country set of bond yields. See \code{data("Yields")} for a detailed data description.}
#' \item{Domestic star factors}{Weighted average of foreign factors. See \code{\link{Transition_Matrix}} for the computation of weights.}
#' }
#' @references Candelon, B. and Moura, R. (2024) "A Multicountry Model of the Term Structures of Interest Rates with a GVAR". (Journal of Financial Econometrics)
#' @keywords Risk Factors GVAR
NULL
