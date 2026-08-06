# ------------------------------------------------------------------------------

#' Constant Forward Curve
#'
#' A dataset containing 2 percent continuously compounded annual interest rate
#' for illustration purposes.
#'
#' @format A vector with 360 elements:
#' \describe{
#'   \item{rate}{discount rate}
#'   ...
#' }
"cForwardCurve"

# ------------------------------------------------------------------------------

#' Fund Map for 10 Funds
#'
#' A dataset containing a default mapping from five indices to ten different
#' funds.
#'
#' @format A matrix with 10 rows and 5 columns:
#' \describe{
#'   \item{index name}{name for each index}
#'   \item{fund number}{proportion of fund allocated to a particular index}
#'   ...
#' }
"fundMap"

# ------------------------------------------------------------------------------

#' Covariance Matrix for 5 Indices
#'
#' A dataset containing the covariance matrix among the returns of five indices.
#'
#' @format A matrix with 5 rows and 5 columns:
#' \describe{
#'   \item{index number}{number for each index}
#'   ...
#' }
"mCov"

# ------------------------------------------------------------------------------

#' 5 Indices for 10 Scenarios over 360 Months
#'
#' A dataset containing a 3D array, number of scenarios (10) by
#' number of indices (5) by number of time steps (360), of Black-Scholes return
#' factors for each index in each of time step and each of scenario.
#'
#' @format A 3D array with dimensions 10x360x5:
#' \describe{
#'   \item{scenario}{scenario number}
#'   \item{month}{month since valuation date}
#'   \item{index number}{monthly return for a particular index in one scenario
#'   one month}
#'   ...
#' }
"indexScen"

# ------------------------------------------------------------------------------

#' Historical Scenario Dates
#'
#' A dataset containing the dates at which historical returns for different
#' indices were observed.
#'
#' @format A vector with 175 elements:
#' \describe{
#'   \item{date}{each observation date of the historical scenarios}
#'   ...
#' }
"histDates"

# ------------------------------------------------------------------------------

#' Historical Index Scenario for 5 Indices over 175 Months
#'
#' A dataset containing a matrix, number of indices (5) by
#' number of time steps (175), of observed historical returns
#' for each index in each of time step in the past.
#'
#' @format A data frame with dimensions 175 rows and 10 columns:
#' \describe{
#'   \item{FIXED}{historical return for index "FIXED" in one month}
#'   \item{INT}{historical return for index "INT" in one month}
#'   \item{MONEY}{historical return for index "MONEY" in one month}
#'   \item{SMALL}{historical return for index "SMALL" in one month}
#'   \item{US}{historical return for index "US" in one month}
#'   ...
#' }
#' @source \url{http://www.math.uconn.edu/~gan/software.html}
#' @section Remark:
#' These historical index scenarios were assessed on 2008-09-12
"histIdxScen"

# ------------------------------------------------------------------------------

#' Index Names
#'
#' A dataset containing names for each index.
#'
#' @format A vector with 5 elements:
#' \describe{
#'   \item{name}{name of the index}
#'   ...
#' }
"indexNames"

# ------------------------------------------------------------------------------

#' Mortality Rate for Male and Female from Ages 5 to 115
#'
#' A dataset containing the mortality rates for male and female from ages 5
#' to 115 (table IAM 1996 from the Society of Actuaries).
#'
#' @format A data frame with 110 rows and 3 columns:
#' \describe{
#'   \item{age}{individual's age}
#'   \item{male}{mortality of a male at a particular age ranging from 5 to 115}
#'   \item{female}{mortality of a female at a particular age ranging
#'   from 5 to 115}
#'   ...
#' }
#' @source \url{https://mort.soa.org}
"mortTable"

# ------------------------------------------------------------------------------

#' Swap Rates across 30 Years
#'
#' A dataset containing US swap rates for various maturities.
#'
#' @format A vector with 8 elements:
#' \describe{
#'   \item{rate}{swap rate}
#'   ...
#' }
#' @source \url{http://www.federalreserve.gov}
#' @section Remark:
#' These swap rates were assessed on 2016-02-08
"swapRate"

# ------------------------------------------------------------------------------

#' A Randomly Generated Pool of Variable Annuities
#'
#' A dataset containing information of the policy and the policy holder.
#'
#' @format A data frame with 19 row and 45 columns:
#' \describe{
#'   \item{recordID}{Unique identifier of the policy}
#'   \item{survivorShip}{Positive weighting number}
#'   \item{gender}{Gender of the policyholder}
#'   \item{productType}{Product type}
#'   \item{issueDate}{Issue date}
#'   \item{matDate}{Maturity date}
#'   \item{birthDate}{Birth date of the policyholder}
#'   \item{currentDate}{Current date}
#'   \item{baseFee}{M&E (Mortality & Expense) fee}
#'   \item{riderFee}{Rider fee}
#'   \item{rollUpRate}{Roll-up rate}
#'   \item{gbAmt}{Guaranteed benefit}
#'   \item{gmwbBalance}{GMWB balance}
#'   \item{wbWithdrawalRate}{Guaranteed withdrawal rate}
#'   \item{withdrawal}{Withdrawal so far}
#'   \item{fundNum1}{Fund number of the 1st investment fund}
#'   \item{fundNum2}{Fund number of the 2nd investment fund}
#'   \item{fundNum3}{Fund number of the 3rd investment fund}
#'   \item{fundNum4}{Fund number of the 4th investment fund}
#'   \item{fundNum5}{Fund number of the 5th investment fund}
#'   \item{fundNum6}{Fund number of the 6th investment fund}
#'   \item{fundNum7}{Fund number of the 7th investment fund}
#'   \item{fundNum8}{Fund number of the 8th investment fund}
#'   \item{fundNum9}{Fund number of the 9th investment fund}
#'   \item{fundNum10}{Fund number of the 10th investment fund}
#'   \item{fundValue1}{Fund value of the 1st investment fund}
#'   \item{fundValue2}{Fund value of the 2nd investment fund}
#'   \item{fundValue3}{Fund value of the 3rd investment fund}
#'   \item{fundValue4}{Fund value of the 4th investment fund}
#'   \item{fundValue5}{Fund value of the 5th investment fund}
#'   \item{fundValue6}{Fund value of the 6th investment fund}
#'   \item{fundValue7}{Fund value of the 7th investment fund}
#'   \item{fundValue8}{Fund value of the 8th investment fund}
#'   \item{fundValue9}{Fund value of the 9th investment fund}
#'   \item{fundValue10}{Fund value of the 10th investment fund}
#'   \item{fundFee1}{Fund management fee of the 1st investment fund}
#'   \item{fundFee2}{Fund management fee of the 2nd investment fund}
#'   \item{fundFee3}{Fund management fee of the 3rd investment fund}
#'   \item{fundFee4}{Fund management fee of the 4th investment fund}
#'   \item{fundFee5}{Fund management fee of the 5th investment fund}
#'   \item{fundFee6}{Fund management fee of the 6th investment fund}
#'   \item{fundFee7}{Fund management fee of the 7th investment fund}
#'   \item{fundFee8}{Fund management fee of the 8th investment fund}
#'   \item{fundFee9}{Fund management fee of the 9th investment fund}
#'   \item{fundFee10}{Fund management fee of the 10th investment fund}
#'   ...
#' }
"VAPort"
