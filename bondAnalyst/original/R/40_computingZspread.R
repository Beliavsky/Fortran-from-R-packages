#'Calculates Z-Spread.
#'@details
#'In fixed-income security analysis, it is important to understand why bond prices and yields-to-maturity change. In this context, it is useful to separate a yield-to-maturity into two components: the benchmark and the spread. The benchmark yield for a fixed-income security with a given time-to-maturity is the base rate which is often a government bond yield. The spread is the difference between the yield-to-maturity and the benchmark. So a Z-spread (zero-volatility spread) is based on the all benchmark spot rates taken from the spot curve. It is the constant spread that is added to each spot rate such that the present value of the cash flows matches the price of the bond (Adams & Smith, 2019).
#'According to information provided by Adams and Smith (2019), the method \code{computingZspread()} is developed to calculate Z-Spread based on the entire benchmark spot curve. Here, \code{coupons} stands for dollar value of the Coupon Payment, \code{mv} is maturity value of the bond, \code{bondPv} is present value or price of the bond, \code{n} is number of years, and \code{spots} is vector of spot rates taken from the spots curve. An output with the value 0.0234 means Z-Spread of 234 bps.
#'@param coupons A number.
#'@param mv A number.
#'@param bondPv A number.
#'@param n A number.
#'@param spots A vector.
#'@return Input values to five arguments  \code{coupons} ,\code{mv},\code{bondPv}, \code{n} and \code{spots}.
#'@importFrom Rdpack reprompt
#'@importFrom stats uniroot
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'computingZspread(coupons=5,mv=100, bondPv=92.38,n=3,spots=c(0.0486,0.0495,0.0565))
#'computingZspread(coupons=50,mv=1000,bondPv=923.8,n=3,spots=c(0.0486,0.0495,0.0565))
#'computingZspread(coupons= 6,mv=100,bondPv=100.75,n=2,spots=c(0.021,0.03635))
#'@export
computingZspread <- function(coupons, mv, bondPv, n, spots){
  fnToFindRoot = function(z) {
    cfv <- mv + (coupons)
    fact <- (cfv/((1 + spots[n]+ z)^n)) - bondPv
    for (i in 1:(n-1)){
      fact <- fact + ((coupons )/((1 + spots[i]+ z)^i))
    }
    return(fact)
  }
    xmin <- uniroot(fnToFindRoot, lower=0, upper=1, extendInt = "yes")$root
   (xmin = round(xmin, digits=4))
}


