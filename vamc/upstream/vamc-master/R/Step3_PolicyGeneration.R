# ------------------------------------------------------------------------------
# ----- Step3_PolicyGeneration.R -----------------------------------------------
# ------------------------------------------------------------------------------

#' Generate Portfolio at Inception
#' 
#' @description Generate a portfolio of VA contracts at inception based on given attribute
#' ranges and investment fund information.
#'
#' @param birthDayRng A vector of two strings in 'YYYY-MM-DD' of birthday range.
#' @param issueRng A vector of two strings in 'YYYY-MM-DD' of issue date range.
#' @param matRng A vector of two integers, range of policy maturity.
#' @param acctValueRng A vector of two doubles, range of initial account values.
#' @param femPct A double, percentage of female policyholders in the
#'   portfolio.
#' @param fundFee A vector of doubles, fees charged by each fund in bps.
#' @param baseFee A double, base fee for all funds in bps.
#' @param prodPct A vector of non-negative doubles, proportions of rider types.
#' @param prodType A vector of strings, names of different rider types.
#' @param riderFee A vector of doubles, rider fees for different riders in bps.
#' @param rollUpRate A vector of doubles, roll up rates for different rider
#'   types in bps.
#' @param withdrawalRate A vector of doubles, withdrawal rates for different
#'   rider types in bps.
#' @param numPolicy An integer, number of each type of policies to be generated.
#' @return Outputs a data frame of 45 columns of attributes in an annuity
#'   contract.
#' @examples
#' genPortInception(c("1980-01-01", "1990-01-01"), c("2001-08-01", "2014-01-01"),
#' c(15, 30), c(5e4, 5e5), 0.4, c(30, 50, 60, 80, 10, 38, 45, 55, 47, 46),
#' 200, rep(1 / 4, 4), c("WBRP", "WBRU", "WBSU", "DBWB"),
#' riderFee = c(25, 35, 35, 50), rep(5, 4), rep(5, 4), 100)
#' \dontrun{
#' genPortInception()
#' }
#' @importFrom stats runif
#' @export
genPortInception <- function(birthDayRng = c("1950-01-01", "1980-01-01"),
                             issueRng = c("2001-08-01", "2014-01-01"),
                             matRng = c(15, 30), acctValueRng = c(5e4, 5e5),
                             femPct = 0.4,
                             fundFee = c(30, 50, 60, 80, 10,
                                         38, 45, 55, 47, 46),
                             baseFee = 200,
                             prodPct = rep(1 / 19, 19),
                             prodType = c("DBRP", "DBRU", "DBSU", "ABRP",
                                          "ABRU", "ABSU", "IBRP", "IBRU",
                                          "IBSU", "MBRP", "MBRU", "MBSU",
                                          "WBRP", "WBRU", "WBSU", "DBAB",
                                          "DBIB", "DBMB", "DBWB"),
                             riderFee = c(25, 35, 35, 50, 60, 60,
                                          60, 70, 70, 50, 60, 60,
                                          65, 75, 75, 75, 85, 75, 90),
                             rollUpRate = rep(5, 19),
                             withdrawalRate = rep(5, 19),
                             numPolicy = 10){
    if (matRng[1] < 0 || matRng[1] > matRng[2]) {
      stop("Please input a positive and increasing maturity range")
    }
    if (acctValueRng[1] < 0 || acctValueRng[1] > acctValueRng[2]) {
      stop("Please input a positive and increasing account value range")
    }
    if (femPct < 0 || femPct > 1) {
      stop("Please input a female percentage between 0 and 1")
    }
    if (sum(prodPct) != 1 || (TRUE %in% (prodPct < 0))) {
      stop("Please input a vector of valid product percentages sum up to 1")
    }
    if ((length(prodType) == length(rollUpRate) &&
         length(rollUpRate) == length(riderFee) &&
         length(riderFee) == length(withdrawalRate)) == FALSE){
      msg <- gsub("[\r\n]", "", "Please input prodType, rollUpRate, riderFee,
 and withdrawalRate in the same length")
      stop(msg)
    }
    if (numPolicy < 1){
      stop("Please input numPolicy > 0")
    }

    numFund <- length(fundFee)
    numType <- length(prodType)
    portSize <- numPolicy * numType

    # recordID, survivorShip, gender, productType, baseFee
    recordID <- 1:portSize
    survivorShip <- rep(1, portSize)
    gender <- sample(c("F", "M"), portSize, replace = T,
                     prob = c(femPct, 1 - femPct))
    productType <- rep(prodType, each = numPolicy)

    # Random issue dates
    issueDate1 <- as.Date(issueRng[1])
    issueDate2 <- as.Date(issueRng[2])
    # sample random dates within range then set them to beginning of month
    issueDate <- sample(seq.Date(issueDate1, issueDate2, by = "day"), portSize)
    issueDate <- as.Date(cut(issueDate, "month"))

    # Random maturity dates
    matDate <- issueDate
    matYear <- sample(seq(matRng[1], matRng[2]),
                      size = portSize, replace = TRUE)
    datePOS <- as.POSIXlt(matDate)
    datePOS$year <- datePOS$year + matYear
    matDate <- as.Date(datePOS)

    # Random birth dates
    birDate1 <- as.Date(birthDayRng[1])
    birDate2 <- as.Date(birthDayRng[2])
    # sample random dates within range then set them to beginning of month
    birthDate <- sample(seq.Date(birDate1, birDate2, by = "day"), portSize)
    birthDate <- as.Date(cut(birthDate, "month"))

    currentDate <- issueDate    # at inception, current date is the issue date
    baseFee <- rep(baseFee, portSize) / 1e4    # constant M&E base fee

    # riderFee, rollUpRate, wbWithdrawalRate all depend on product type
    names(riderFee) <- prodType
    riderFee <- riderFee[productType] / 1e4
    names(rollUpRate) <- prodType
    rollUpRate <- rollUpRate[productType] / 1e4
    names(withdrawalRate) <- prodType
    wbWithdrawalRate <- withdrawalRate[productType] / 1e4

    # Initial values at inception
    gbAmt <- rep(0, portSize)
    gmwbBalance <- rep(0, portSize)
    withdrawal <- rep(0, portSize)

    portfolio <- data.frame(recordID, survivorShip, gender, productType,
                            issueDate, matDate, birthDate, currentDate,
                            baseFee, riderFee, rollUpRate, gbAmt, gmwbBalance,
                            wbWithdrawalRate, withdrawal)

    # Fund info calculation
    # -- fundNum, unclear what the usage is at the moment
    fundNum <- matrix(rep(1:numFund, each = portSize), nrow = portSize)
    colnames(fundNum) <- paste0("fundNum", sprintf("%d", 1:numFund))

    # -- fundValue, at inception, equally allocate account value to
    #    selected funds
    fundValue <- matrix(0, nrow = portSize, ncol = numFund)
    colnames(fundValue) <- paste0("fundValue", sprintf("%d", 1:numFund))

    acctVal <- runif(portSize, min = acctValueRng[1], max = acctValueRng[2])
    # no. of funds to be invested
    nSelected <- ceiling(runif(portSize, max = numFund))

    for (i in 1:portSize) {
        selectedFunds <- sample(1:numFund, nSelected[i], replace = FALSE)
        fundValue[i, selectedFunds] <- acctVal[i] / nSelected[i]
    }

    # -- fundFee, same for all policies
    fundFee <- matrix(rep(fundFee, each = portSize), nrow = portSize) / 1e4
    colnames(fundFee) <- paste0("fundFee", sprintf("%d", 1:numFund))

    return(data.frame(cbind(portfolio, fundNum, fundValue, fundFee)))
}
