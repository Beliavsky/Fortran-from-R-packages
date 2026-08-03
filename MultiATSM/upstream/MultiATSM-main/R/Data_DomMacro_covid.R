#' @title Data: Risk Factors for the GVAR - Candelon and Moura (2023)
#'
#' @description Domestic risk factors data used in the GVAR models - Candelon and Moura (2023)
#' @name DomMacro_covid
#' @aliases DomMacro_covid
#' @docType data
#' @usage data("DomMacro_covid")
#' @format A matrix of country-specific risk factors (inflation, output growth, CDS, and COVID-19 reproduction rate) for Brazil, India, Mexico, and Russia.
#' The data have weekly frequency and span the period from March 22, 2020, to September 26, 2021.
#' @source
#' \describe{
#' \item{Inflation}{Monthly CPI (from OECD) interpolated to daily data (spline), converted to weekly year-over-year changes, and detrended <https://www.oecd.org/en/data/indicators/inflation-cpi.html>}
#' \item{Output growth}{ Detrended weekly estimate of GDP year-over-year growth derived from the OECD Weekly Tracker index <https://web-archive.oecd.org/sections/weekly-tracker-of-gdp-growth/index.htm>}
#' \item{CDS}{ 5-year maturity CDS. Simulated data constructed using Bloomberg bond yield series.}
#' \item{COVID-19 reproduction rate}{detrended R rate from the Our World in Data database <https://ourworldindata.org/coronavirus>}
#' }
#' @references Candelon, B. and Moura, R. (2023) "Sovereign yield curves and the COVID-19 in emerging markets". (Economic Modelling)
#' @keywords domestic risk factors
NULL
