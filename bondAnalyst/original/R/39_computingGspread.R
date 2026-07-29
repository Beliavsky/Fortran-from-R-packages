#'Calculates the G-Spread which is the spread between the yields-to-maturity on the corporate bond and that of government bond having the same maturity.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{computingGspread()} is developed to calculate G-Spread for given values of yields-to-maturity on the corporate bond and that of government bond having the same maturity. Here, \code{ytmCorpBond} stands for yields-to-maturity on the corporate bond and \code{ytmBenchGovtBond} denotes yields-to-maturity on government bond with the same maturity. An output with the value 0.02327 means G-Spread of 232.7 bps.
#'@param ytmCorpBond A number.
#'@param ytmBenchGovtBond A number.
#'@return Input values to two arguments  \code{ytmCorpBond} and \code{ytmBenchGovtBond}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'computingGspread(ytmCorpBond=0.05932, ytmBenchGovtBond=0.03605)
#'@export
computingGspread <- function(ytmCorpBond, ytmBenchGovtBond){
  (gspread_val <- ytmCorpBond- ytmBenchGovtBond)
}

