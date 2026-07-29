#'Calculates the accrued interest with 30-by-360, day convention.
#'@details
#'There are different conventions used in bond markets to count days for valuation purposes. The two most common day-count conventions are actual-by-actual and 30-by-360. The 30-by-360 day-count convention often is used on corporate bonds. It assumes that each month has 30 days and that a full year has 360 days. Therefore, for this method, there are assumed to be 42 days between 15 May and 27 June: 15 days between 15 May and 30 May and 27 days between 1 June and 27 June. There are assumed to be 180 days (for being used as days between two semi-annual coupon payments) in the six-month period between 15 May and 15 November. So, the accrued interest at a 4.375 percent, for a semi-annual payment corporate bond is 0.510417 per 100 of par value (Adams & Smith, 2019).
#'Based on this, the method \code{aiRoundedDaysConv()} is developed to compute the accrued interest with 30-by-360 day convention. So, \code{aiRoundedDaysConv()} gives accrued interest with 30-by-360 day convention for values passed to its five arguments. Here, \code{cPmt} is dollar value of Coupon Payment, \code{bfrStlDt} is date of interest payment before settlement date, \code{stlDt} is Settlement Date, \code{elpsMnths} denotes number of months elapsed from day of last coupon payment to the \code{bfrStlDt}(so in this example number of months from last coupon payment to June 15, 2019), and \code{daysBtwnCpns} denotes number of days between two coupon payments. The output is rounded off to six decimal places.
#'@param cPmt A number.
#'@param bfrStlDt A character.
#'@param stlDt A character.
#'@param elpsMnths A number.
#'@param daysBtwnCpns A number.
#'@return Input values to four arguments  \code{cPmt} , \code{bfrStlDt}  \code{stlDt} and \code{elpsMnths}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'aiRoundedDaysConv(4.375, as.Date("2019-6-15"), as.Date("2019-6-27"),1,180)
#'@export
aiRoundedDaysConv <- function(cPmt,bfrStlDt,stlDt,elpsMnths,daysBtwnCpns){
  daysFromLastCoupon = stlDt - bfrStlDt
  acruredInterest <- (cPmt/2)*(as.numeric(daysFromLastCoupon+(30*elpsMnths))/daysBtwnCpns)
  (acruredInterest = round(acruredInterest, digits=6))
}


