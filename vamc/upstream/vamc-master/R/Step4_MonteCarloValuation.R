# ------------------------------------------------------------------------------
# ----- Step4_MonteCarloValuation.R --------------------------------------------
# ------------------------------------------------------------------------------

#' @importFrom utils read.csv
readMortTable <- function(dirFileName, header = T){
    # reads mortality table from .csv file specified by dirFileName
    # colname indicates whether column names are included in the file
    # By default, the first column is age x, the second and third columns are
    # mortality rates qx for male & female
    # Inputs:
    # -- dirFileName: string of directory and filename
    # -- colname: binary indicator for whether column name is given in the file
    # Output:
    # -- mortTable: a dataframe with three columns, x and qx for male & female
    # -- ALSO need to overwrite existing default mortTable in the workspace

    mortTable <- read.csv(file = dirFileName, header = header, sep = ",")
    colnames(mortTable) <- c("age", "F", "M")

    # mortTable is a dataframe with three columns, x and qx for male & female
    return(mortTable)
}

# ------------------------------------------------------------------------------
#' Calculate Mortality Factors
#' 
#' @description Calculates the mortality factors (t - 1)px q(x + t - 1) and tpx required to
#' valuate the inPolicy. Extract gender, age (birth date & current date),
#' valuation date (current date), and maturity date from inPolicy, mortality
#' rates from mortTable.
#'
#' @param inPolicy A vector containing 45 attributes of a VA policy,
#'   usually a row of a VA portfolio dataframe.
#' @param mortTable A dataframe with three columns of doubles representing the
#'   mortality table.
#' @param dT A double of stepsize in years; dT = 1 / 12 would be monthly.
#' @return Outputs a two-column data frame of doubles of mortFactors (t - 1)px
#'   q(x + t - 1) and tpx.
#' @examples
#' exPolicy <- VAPort[1, ]
#' calcMortFactors(exPolicy, mortTable, dT = 1 / 12)
#' @export
calcMortFactors <- function(inPolicy, mortTable, dT = 1 / 12){

    birDate <- as.POSIXlt(inPolicy[1, "birthDate"])
    curDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    matDate <- as.POSIXlt(inPolicy[1, "matDate"])
    if (curDate < birDate) stop("Current date is prior to birth date.")
    if (matDate < birDate) stop("Maturity date is prior to birth date.")
    if (matDate < curDate) stop("Maturity date is prior to current date.")

    # no. of months sine birth to current date,
    # all dates are first of the month
    monBirCur <- 12 * (curDate$year - birDate$year) +
      (curDate$mon - birDate$mon)
    # no. of months sine birth to maturity date,
    # all dates are first of the month
    monMatCur <- 12 * (matDate$year - birDate$year) +
      (matDate$mon - birDate$mon)
    ageCur <- monBirCur / 12              # age at current date
    x0 <- floor(ageCur)                 # integer age at current date
    t0 <- ageCur - x0              # fractional age at current date
    xT <- floor(monMatCur / 12)     # integer age after numStep periods

    gender <- ifelse(inPolicy[1, "gender"] == "F", "female", "male")

    # annual mortality rates from current age to maturity age
    maxRng <- c(x0:mortTable$age[length(mortTable$age)])
    ageRngLen <- xT - x0 + 1
    annualQ <- mortTable[mortTable$age %in% maxRng, gender]
    # append some 1's if xMat > max age in mortality table
    annualQ <- c(annualQ, rep(1, as.numeric(ageRngLen > length(maxRng)) *
                                (ageRngLen - length(annualQ))))

    p <- rep(1)     # placeholder for the p factors, default 1 for 0px
    pq <- rep(1)    # placeholder for the pq factors

    # loop over for numStep with stepsize dT
    xCur <- x0
    tCur <- t0
    dP <- 1
    step <- 1
    while (dP > 0.00001) {
        baseQ <- annualQ[xCur - x0 + 1]   # base qx in current loop
        # prob. of death within next dT at age xCur+tCur
        dtQxt <- (dT * baseQ) / (1 - tCur * baseQ)
        # prob. of surviving the current dT
        p[step + 1] <- p[step] * (1 - dtQxt)
        dP <- p[step + 1]
        pq[step] <- p[step] * dtQxt # prob. of death in the current dT

        # update age for next period
        # age at current date
        ageCur <- ageCur + dT
        # integer age at current date
        xCur <- floor(ageCur + 1e-10)
        # fractional age at current date
        tCur <- ifelse(abs(ageCur - xCur) < 1e-10, 0, ageCur - xCur)
        step <- step + 1
    }

    return(mortFactors = data.frame(pq = pq, p = p[-1]))
}

# ------------------------------------------------------------------------------
# Remarks:
# 1. projectDBRP is 1 of the 19 cash flow projection functions
# 2. Each cash flow projection function follows calculations in
# Appendix A of the paper
# 3. User-defined VA policies should have an associated projectXXXX function
# Project a DBRP policy whose attributes are specified in inPolicy based on
# one fund scenario oneFundScen and steplength dT.
projectDBRP <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
      (length(p) >= numStep) &&
      (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV * inPolicy[1, "riderFee"] -
        dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
    LA[step] <- 0
    RC[step] <- dFee[step]
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectDBRU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)
    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV * inPolicy[1, "riderFee"] -
        dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary: annual roll-up
      inPolicy[1, "gbAmt"] <- inPolicy[1, "gbAmt"] *
        (1 + inPolicy[1, "rollUpRate"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
    LA[step] <- 0
    RC[step] <- dFee[step]
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectDBSU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon){
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
    LA[step] <- 0
    RC[step] <- dFee[step]
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectABRP <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV

    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }

    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # renewal
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      if (inPolicy[1, "gbAmt"] < dAV[step]) inPolicy[1, "gbAmt"] <- dAV[step]
      else {
        LA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
        numFund <- length(inPolicy[1, grep("fundValue", colnames(inPolicy))])
        if (dAV[step] > 0.00001) {
          inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
            inPolicy[1, grep("fundValue", colnames(inPolicy))] *
            inPolicy[1, "gbAmt"] / dAV[step]
        }else {
            inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
              inPolicy[1, "gbAmt"] / numFund
          }
        }
      }
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + (as.POSIXlt(inPolicy[1, "currentDate"]) -
                                    as.POSIXlt(inPolicy[1, "issueDate"]))
    inPolicy[1, "issueDate"] <- inPolicy[1, "currentDate"]
    inPolicy[1, "matDate"] <- as.Date(newDate)

    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectABRU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary: annual roll-up
      inPolicy[1, "gbAmt"] <- inPolicy[1, "gbAmt"] *
        (1 + inPolicy[1, "rollUpRate"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }

    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # renewal
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      if (inPolicy[1, "gbAmt"] < dAV[step]) inPolicy[1, "gbAmt"] <- dAV[step]
      else {
        LA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
        numFund <- length(inPolicy[1, grep("fundValue", colnames(inPolicy))])
        if (dAV[step] > 0.00001) {
          inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
            inPolicy[1, grep("fundValue", colnames(inPolicy))] *
            inPolicy[1, "gbAmt"] / dAV[step]
        }else {
          inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
            inPolicy[1, "gbAmt"] / numFund
        }
      }
    }
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + (as.POSIXlt(inPolicy[1, "currentDate"]) -
                                    as.POSIXlt(inPolicy[1, "issueDate"]))
    inPolicy[1, "issueDate"] <- inPolicy[1, "currentDate"]
    inPolicy[1, "matDate"] <- as.Date(newDate)

    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectABSU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon){
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }

    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # renewal
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]){
      if (inPolicy[1, "gbAmt"] < dAV[step]) inPolicy[1, "gbAmt"] <- dAV[step]
      else {
        LA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
        numFund <- length(inPolicy[1, grep("fundValue", colnames(inPolicy))])
        if (dAV[step] > 0.00001){
          inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
            inPolicy[1, grep("fundValue", colnames(inPolicy))] *
            inPolicy[1, "gbAmt"] / dAV[step]
        }else {
          inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
            inPolicy[1, "gbAmt"] / numFund
        }
      }
    }
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + (as.POSIXlt(inPolicy[1, "currentDate"]) -
                                    as.POSIXlt(inPolicy[1, "issueDate"]))
    inPolicy[1, "issueDate"] <- inPolicy[1, "currentDate"]
    inPolicy[1, "matDate"] <- as.Date(newDate)

    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectIBRP <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # calculate guarateed price of an annuity with payments of $1 per annum
  # with r = 5%
  ag <- 0
  dP <- 1
  nY <- 0
  r <- 0.05
  while (dP > 0.00001) {
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    ag <- ag + dP * exp(-r * nY)
    nY <- nY + 1
  }
  # calculate market price of an annuity with payments of $1 per annum
  # beginning at time T
  aT <- 0
  dP <- 1
  nY <- 0
  while (dP > 0.00001) {
    dFR <- 0
    if (nY * 12 < numStep) dFR <- df[nY * 12 + 1]
    else dFR <- df[numStep - 1]
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    aT <- aT + dP * exp(-dFR * nY)
    nY <- nY + 1
  }
  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon){
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV

    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      LA[step] <- max(inPolicy[1, "gbAmt"] * aT / ag - dAV[step], 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation

  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectIBRU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # calculate guarateed price of an annuity with payments of $1 per annum
  # with r = 5%
  ag <- 0
  dP <- 1
  nY <- 0
  r <- 0.05
  while (dP > 0.00001) {
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    ag <- ag + dP * exp(-r * nY)
    nY <- nY + 1
  }
  # calculate market price of an annuity with payments of $1 per annum
  # beginning at time T
  aT <- 0
  dP <- 1
  nY <- 0
  while (dP > 0.00001){
    dFR <- 0
    if (nY * 12 < numStep) dFR <- df[nY * 12 + 1]
    else dFR <- df[numStep - 1]
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    aT <- aT + dP * exp(-dFR * nY)
    nY <- nY + 1
  }
  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary: annual roll-up
      inPolicy[1, "gbAmt"] <- inPolicy[1, "gbAmt"] *
        (1 + inPolicy[1, "rollUpRate"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      LA[step] <- max(inPolicy[1, "gbAmt"] * aT / ag - dAV[step], 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation

  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectIBSU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # calculate guarateed price of an annuity with payments of $1 per annum
  # with r = 5%
  ag <- 0
  dP <- 1
  nY <- 0
  r <- 0.05
  while (dP > 0.00001) {
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    ag <- ag + dP * exp(-r * nY)
    nY <- nY + 1
  }
  # calculate market price of an annuity with payments of $1 per annum
  # beginning at time T
  aT <- 0
  dP <- 1
  nY <- 0
  while (dP > 0.00001) {
    dFR <- 0
    if (nY * 12 < numStep) dFR <- df[nY * 12 + 1]
    else dFR <- df[numStep - 1]
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    aT <- aT + dP * exp(-dFR * nY)
    nY <- nY + 1
  }
  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      LA[step] <- max(inPolicy[1, "gbAmt"] * aT / ag - dAV[step], 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation

  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectMBRP <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon){
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV

    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      LA[step] <- max(inPolicy[1, "gbAmt"] - dAV, 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectMBRU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary: annual roll-up
      inPolicy[1, "gbAmt"] <- inPolicy[1, "gbAmt"] *
        (1 + inPolicy[1, "rollUpRate"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      LA[step] <- max(inPolicy[1, "gbAmt"] - dAV, 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectMBSU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon){
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]){
      LA[step] <- max(inPolicy[1, "gbAmt"] - dAV, 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectWBRP <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }
  numFund <- dim(oneFundScen)[2]

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon){
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV

      #update annual withdrawal
      dWAG <- inPolicy[1, "gbAmt"] * inPolicy[1, "wbWithdrawalRate"]
      dWA <- min(dWAG, inPolicy[1, "gmwbBalance"])
      inPolicy[1, "gmwbBalance"] <- inPolicy[1, "gmwbBalance"] - dWA
      inPolicy[1, "withdrawal"] <- inPolicy[1, "withdrawal"] + dWA

      if (inPolicy[1, "gmwbBalance"] < 0.0001) inPolicy[1, "gbAmt"] <- 0

      dAV[step] <- max(0, dAV[step] - dWA)
      if (dAV[step] > 0.00001) {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
          inPolicy[1, grep("fundValue", colnames(inPolicy))] *
          (dAV[step] / (dAV[step] + dWA))
      } else {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <- rep(0, numFund)
      }
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- max(0, dWA - dAV[step])
    RC[step] <- dFee[step]
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}
# ------------------------------------------------------------------------------
projectWBRU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }
  numFund <- dim(oneFundScen)[2]

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary: annual roll-up
      inPolicy[1, "gbAmt"] <- inPolicy[1, "gbAmt"] *
        (1 + inPolicy[1, "rollUpRate"])
      #update annual withdrawal
      dWAG <- inPolicy[1, "gbAmt"] * inPolicy[1, "wbWithdrawalRate"]
      dWA <- min(dWAG, inPolicy[1, "gmwbBalance"])
      inPolicy[1, "gmwbBalance"] <- inPolicy[1, "gmwbBalance"] - dWA
      inPolicy[1, "withdrawal"] <- inPolicy[1, "withdrawal"] + dWA

      if (inPolicy[1, "gmwbBalance"] < 0.0001) inPolicy[1, "gbAmt"] <- 0

      dAV[step] <- max(0, dAV[step] - dWA)
      if (dAV[step] > 0.00001) {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
          inPolicy[1, grep("fundValue", colnames(inPolicy))] *
          (dAV[step] / (dAV[step] + dWA))
      } else {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <- rep(0, numFund)
      }
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- max(0, dWA - dAV[step])
    RC[step] <- dFee[step]
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectWBSU <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }
  numFund <- dim(oneFundScen)[2]

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
      #update annual withdrawal
      dWAG <- inPolicy[1, "gbAmt"] * inPolicy[1, "wbWithdrawalRate"]
      dWA <- min(dWAG, inPolicy[1, "gmwbBalance"])
      inPolicy[1, "gmwbBalance"] <- inPolicy[1, "gmwbBalance"] - dWA
      inPolicy[1, "withdrawal"] <- inPolicy[1, "withdrawal"] + dWA

      if (inPolicy[1, "gmwbBalance"] < 0.0001) inPolicy[1, "gbAmt"] <- 0

      dAV[step] <- max(0, dAV[step] - dWA)
      if (dAV[step] > 0.00001) {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
          inPolicy[1, grep("fundValue", colnames(inPolicy))] *
          (dAV[step] / (dAV[step] + dWA))
      } else {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <- rep(0, numFund)
      }
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- 0
    LA[step] <- max(0, dWA - dAV[step])
    RC[step] <- dFee[step]
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectDBAB <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }

    DA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
    LA[step] <- 0
    RC[step] <- dFee[step]
    # renewal
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      if (inPolicy[1, "gbAmt"] < dAV[step]) inPolicy[1, "gbAmt"] <- dAV[step]
      else {
        LA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
        numFund <- length(inPolicy[1, grep("fundValue", colnames(inPolicy))])
        if (dAV[step] > 0.00001) {
          inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
            inPolicy[1, grep("fundValue", colnames(inPolicy))] *
            inPolicy[1, "gbAmt"] / dAV[step]
        } else {
          inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
            inPolicy[1, "gbAmt"] / numFund
        }
      }
    }
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + (as.POSIXlt(inPolicy[1, "currentDate"]) -
                                    as.POSIXlt(inPolicy[1, "issueDate"]))
    inPolicy[1, "issueDate"] <- inPolicy[1, "currentDate"]
    inPolicy[1, "matDate"] <- as.Date(newDate)

    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectDBIB <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # calculate guarateed price of an annuity with payments of $1 per annum
  # with r = 5%
  ag <- 0
  dP <- 1
  nY <- 0
  r <- 0.05
  while (dP > 0.00001) {
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    ag <- ag + dP * exp(-r * nY)
    nY <- nY + 1
  }
  # calculate market price of an annuity with payments of $1 per annum
  # beginning at time T
  aT <- 0
  dP <- 1
  nY <- 0
  while (dP > 0.00001) {
    dFR <- 0
    if (nY * 12 < numStep) dFR <- df[nY * 12 + 1]
    else dFR <- df[numStep - 1]
    # fetch p at j+1th step
    if (nY * 12 < length(p)) dP <- p[nY * 12 + 1]
    else dP <- 0
    aT <- aT + dP * exp(-dFR * nY)
    nY <- nY + 1
  }
  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon){
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
    } else{
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      LA[step] <- max(inPolicy[1, "gbAmt"] * aT / ag - dAV[step], 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation

  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectDBMB <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
    LA[step] <- 0
    RC[step] <- dFee[step]
    # at maturity date
    if (inPolicy[1, "matDate"] == inPolicy[1, "currentDate"]) {
      LA[step] <- max(inPolicy[1, "gbAmt"] - dAV, 0)
    }
    # update currentDate
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
projectDBWB <- function(inPolicy, oneFundScen, dT = 1 / 12, pq, p, df){
  # projection will always start from the current date for numStep periods
  numStep <- dim(oneFundScen)[1]
  if (((length(pq) >= numStep) &&
       (length(p) >= numStep) &&
       (length(df) >= numStep)) == FALSE) {
    stop("df, pq, and p must have length > numStep from oneFundScen")
  }
  numFund <- dim(oneFundScen)[2]

  # Initialize DA, LA, RC
  DA <- rep(0, numStep)
  LA <- rep(0, numStep)
  RC <- rep(0, numStep)

  # These two variables are defined within the function only
  dAV <- rep(0, numStep)
  dFee <- rep(0, numStep)

  # The main for loop to do the calculations, hard to parallelize
  for (step in 1:numStep) {
    # one suggestion for partial account evolution to avoid fund_helper
    dPartialAV <- inPolicy[1, grep("fundValue", colnames(inPolicy))] *
      oneFundScen[step, ] *
      (1 - inPolicy[1, grep("fundFee", colnames(inPolicy))] * dT)

    #
    if (as.POSIXlt(inPolicy[1, "currentDate"])$mon ==
        as.POSIXlt(inPolicy[1, "issueDate"])$mon) {
      raw <- sum(dPartialAV)
      dFee[step] <- sum(raw * inPolicy[1, "riderFee"])

      #update individual fund value
      dPartialAV <- dPartialAV - dPartialAV *
        inPolicy[1, "riderFee"] - dPartialAV * inPolicy[1, "baseFee"]
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
      # anniversary
      inPolicy[1, "gbAmt"] <- max(dAV[step], inPolicy[1, "gbAmt"])
      #update annual withdrawal
      dWAG <- inPolicy[1, "gbAmt"] * inPolicy[1, "wbWithdrawalRate"]
      dWA <- min(dWAG, inPolicy[1, "gmwbBalance"])
      inPolicy[1, "gmwbBalance"] <- inPolicy[1, "gmwbBalance"] - dWA
      inPolicy[1, "withdrawal"] <- inPolicy[1, "withdrawal"] + dWA

      if (inPolicy[1, "gmwbBalance"] < 0.0001) inPolicy[1, "gbAmt"] <- 0

      dAV[step] <- max(0, dAV[step] - dWA)
      if (dAV[step] > 0.00001) {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <-
          inPolicy[1, grep("fundValue", colnames(inPolicy))] *
          (dAV[step] / (dAV[step] + dWA))
      } else {
        inPolicy[1, grep("fundValue", colnames(inPolicy))] <- rep(0, numFund)
      }
    } else {
      dAV[step] <- sum(dPartialAV)
      inPolicy[1, grep("fundValue", colnames(inPolicy))] <- dPartialAV
    }
    # update policy info
    DA[step] <- max(0, inPolicy[1, "gbAmt"] - dAV[step])
    LA[step] <- max(0, dWA - dAV[step])
    RC[step] <- dFee[step]
    newDate <- as.POSIXlt(inPolicy[1, "currentDate"])
    newDate$mon <- newDate$mon + 1
    inPolicy[1, "currentDate"] <- as.Date(newDate)
  }

  # return death benefit, living benefit, and risk charge
  # outPolicy is used in aging, but not in valuation
  pq <- pq[1:numStep] * df[1:numStep]
  p <- p[1:numStep] * df[1:numStep]
  return(list(DA = pq %*% DA, LA = p %*% LA, RC = p %*% RC,
              outPolicy = inPolicy))
}

# ------------------------------------------------------------------------------
#' Valuate One Policy
#' 
#' @description Valuate a VA policy specified in inPolicy based on the simulated fund
#' scenarios fundScen. The time step length is specified in dT and the
#' discount rate for each period is specified in df.
#' @param inPolicy A vector containing 45 attributes of a VA policy,
#'   usually a row of a VA portfolio dataframe.
#' @param mortTable A dataframe with three columns of doubles representing the
#'   mortality table.
#' @param fundScen A numScen-by-numStep-by-numFund array of doubles of
#'   return factors (i.e., exp(mu_t dt)) in each period.
#' @param dT A double of stepsize in years; dT = 1 / 12 would be monthly.
#' @param df A vector of doubles of risk-free discount rates of different tenor
#'   (not forward rates), should have length being numStep.
#' @return Outputs a list of doubles of policyValue, the average discounted
#'   payoff of the VA, and riskCharge, the average discounted risk charges.
#' @examples
#' fundScen <- genFundScen(fundMap, indexScen)[1, , ]
#' exPolicy <- VAPort[1, ]
#' valuateOnePolicy(exPolicy, mortTable, fundScen, 1 / 12, cForwardCurve)
#' @export
valuateOnePolicy <- function(inPolicy, mortTable, fundScen, dT = 1 / 12, df){
  if (inPolicy[1, "matDate"] <= inPolicy[1, "currentDate"]) {
    return(policyValue = 0)
  }

  # If the policy is inforce, select the risk project function to calculate
  # the death benefits, living benefits, and risk charges for each period
  Type <- inPolicy[1, "productType"]

  # projectFun is projectXXXX where XXXX is the product type of the policy
  projectFun <- get(paste0("project", Type))

  # In valuation, project from current date to maturity date
  startDate <- as.POSIXlt(inPolicy[1, "currentDate"])
  endDate <- as.POSIXlt(inPolicy[1, "matDate"])
  numStep <- 12 * (endDate$year - startDate$year) +
    (endDate$mon - startDate$mon)

  scenDim <- dim(fundScen)
  # Check if only one scenario is provided
  if (length(scenDim) == 2){
    numScen <- 1
  } else {
    numScen <- scenDim[1]
  }

  # Calculate actuarial discount factors
  mortFactors <- calcMortFactors(inPolicy, mortTable, dT)
  pq <- matrix(mortFactors[, "pq"], nrow = 1)  # enforced row vector
  p <- matrix(mortFactors[, "p"], nrow = 1)  # enforced row vector

  # Placeholders for results
  DA <- c()
  LA <- c()
  RC <- c()

  # this for loop is parallelizable
  for (scen in 1:numScen) {
    if (numScen == 1) {
      curScen <- matrix(fundScen[1:numStep, ], nrow = numStep)
      VABenefits <- projectFun(inPolicy, curScen, dT, pq, p, df)
      DA[scen] <- VABenefits$DA
      LA[scen] <- VABenefits$LA
      RC[scen] <- VABenefits$RC
    } else {
    # curScen is numStep-by-numFund
      curScen <- matrix(fundScen[scen, 1:numStep, ], nrow = numStep)
      VABenefits <- projectFun(inPolicy, curScen, dT, pq, p, df)
      DA[scen] <- VABenefits$DA
      LA[scen] <- VABenefits$LA
      RC[scen] <- VABenefits$RC
    }
  }

  # Calculate discounted payoffs for all scenarios (1-by-numScen vectors)
  DAs <- sum(DA)
  LAs <- sum(LA)
  RCs <- sum(RC)

  return(list(policyValue = (DAs + LAs) / numScen,
              riskCharge = RCs / numScen))
}

# ------------------------------------------------------------------------------
#' Age One Policy
#'
#' @description Age a VA policy specified in inPolicy from currentDate (specified in
#' inPolicy) to targetDate. The againg scenario is given in fundScen.
#' The time step length is specified in dT.
#' Here we input a rather irrelevant parameter df to "hack" for a more
#' flexible user-defined projection function.
#' @param inPolicy A vector containing 45 attributes of a VA policy,
#'   usually a row of a VA portfolio dataframe.
#' @param mortTable A dataframe with three columns of doubles representing the
#'   mortality table.
#' @param fundScen A numScen-by-numStep-by-numFund array of doubles of
#'   return factors (i.e., exp(mu_t dt)) in each period.
#' @param scenDates A vector containing strings in the format of "YYYY-MM-DD"
#'   of dates corresponding to each period in fundScen.
#' @param dT A double of stepsize in years; dT = 1 / 12 would be monthly.
#' @param targetDate A string in the format of "YYYY-MM-DD" of valuation date
#'   of the portfolio.
#' @param df A vector of doubles of risk-free discount rates of different tenor
#'   (not forward rates), should have length being numStep.
#' @return Outputs a vector containing 45 attributes of a VA policy, where
#'   currentDate, gbAmt, GMWBbalance, withdrawal, & fundValue could be updated
#'   as a result of aging. Usually a row of a VA portfolio dataframe.
#' @examples
#' exPolicy <- VAPort[1, ]
#' targetDate <- "2016-01-01"
#' histFundScen <- genFundScen(fundMap, histIdxScen)
#' ageOnePolicy(exPolicy, mortTable, histFundScen, histDates, dT = 1 / 12,
#' targetDate, cForwardCurve)
#' \dontrun{
#' targetDate <- "2001-01-01"
#' histFundScen <- genFundScen(fundMap, histIdxScen)
#' ageOnePolicy(exPolicy, mortTable, histFundScen, histDates, dT = 1 / 12,
#' targetDate, cForwardCurve)
#' }
#' \dontrun{
#' exPolicy <- VAPort[1, ]
#' exPolicy[1, c("currentDate", "issueDate")] <- c("2001-01-01", "2001-01-01")
#' histFundScen <- genFundScen(fundMap, histIdxScen)
#' ageOnePolicy(exPolicy, mortTable, histFundScen, histDates, dT = 1 / 12,
#' targetDate, cForwardCurve)
#' }
#' @section Note:
#' Target date MUST be PRIOR to the last date of historical scenario date,
#' Current date MUST be LATER than the first date of historical scenario date.
#' @export
ageOnePolicy <- function(inPolicy, mortTable, fundScen, scenDates,
                         dT = 1 / 12, targetDate, df){
  scenStartDate <- as.POSIXlt(scenDates[1])
  targetDate <- as.Date(targetDate)
  if (targetDate <= inPolicy[1, "currentDate"]) {
    warning("No aging required. Original policy returned.")
    return(inPolicy)
  }
  if (targetDate > max(scenDates)) {
    msg <- "No aging performed. Target date beyond the last date of
historical scenario date 2016-02-01"
    msg <- gsub("[\r\n]", " ", msg)
    stop(msg)
  }
  if (min(scenDates) > inPolicy[1, "currentDate"]) {
    msg <- "No aging performed. Current date before the first date of
historical scenario date 2001-08-01"
    msg <- gsub("[\r\n]", " ", msg)
    stop(msg)
  }
  rollEndDate <- as.POSIXlt(targetDate)
  rollEndDate$mon <- rollEndDate$mon - 1
  # If the policy is inforce, select the risk project function to calculate
  # the death benefits, living benefits, and risk charges for each period
  type <- inPolicy[1, "productType"]

  # projectFun is projectXXXX where XXXX is the product type of the policy
  projectFun <- get(paste0("project", type))

  # In valuation, project from current date to maturity date
  startDate <- as.POSIXlt(inPolicy[1, "currentDate"])
  endDate <- as.POSIXlt(inPolicy[1, "matDate"])
  numStep <- 12 * (endDate$year - startDate$year) +
    (endDate$mon - startDate$mon)

  agingStartIndx <- 12 * (startDate$year - scenStartDate$year) +
    (startDate$mon - scenStartDate$mon)
  agingEndIndx <- 12 * (rollEndDate$year - scenStartDate$year) +
    (rollEndDate$mon - scenStartDate$mon)
  rngScenIndx <- c(agingStartIndx:agingEndIndx)
  curScen <- fundScen[rngScenIndx, ]

  # Calculate actuarial discount factors
  mortFactors <- calcMortFactors(inPolicy, mortTable, dT)
  actFactors <- mortFactors * df[1:numStep]
  pq <- matrix(actFactors[, "pq"], nrow = 1)  # enforced row vector
  p <- matrix(actFactors[, "p"], nrow = 1)  # enforced row vector
  VABenefits <- projectFun(inPolicy, curScen, dT, pq, p, df)

  # Update different fields in the output policy
  outPolicy <- VABenefits$outPolicy

  return(outPolicy)
}

# ------------------------------------------------------------------------------
#' Valuate a Portfolio
#' 
#' @description Valuate a portfolio VA policies specified in each curPolicy of inPortfolio
#' based on the simulated fund scenarios fundScen.
#' The time step length is specified in dT and the discount rate for each period
#' is specified in df.
#' @param inPortfolio A dataframe containing numPolicy rows and 45 attributes
#'   of each VA policy.
#' @param mortTable A dataframe with three columns of doubles representing the
#'   mortality table.
#' @param fundScen A numScen-by-numStep-by-numFund array of doubles of
#'   return factors (i.e., exp(mu_t dt)) in each period.
#' @param dT A double of stepsize in years; dT = 1 / 12 would be monthly.
#' @param df A vector of doubles of risk-free discount rates of different tenor
#'   (not forward rates), should have length being numStep.
#' @return Outputs a list of doubles of portVal, the sum of average discounted
#'   payoff of the VAs in inPortfolio, portRC, the sum of average discounted
#'   risk charges of the VAs in inPortfolio, and vectors of doubles of these
#'   average discounted values for each policy.
#' @examples
#' fundScen <- genFundScen(fundMap, indexScen)[1, , ]
#' valuatePortfolio(VAPort[1:2, ], mortTable, fundScen, 1 / 12, cForwardCurve)
#' @export
valuatePortfolio <- function(inPortfolio, mortTable, fundScen, dT = 1 / 12, df){
  numPolicy <- nrow(inPortfolio)
  vecVal <- rep(0, numPolicy)
  vecRC <- rep(0, numPolicy)
  for (i in 1:numPolicy){
    curPolicy <- inPortfolio[i, ]
    outTemp <- valuateOnePolicy(curPolicy, mortTable, fundScen, dT, df)
    vecVal[i] <- outTemp$policyValue
    vecRC[i] <- outTemp$riskCharge
  }
  return (list(portVal = sum(vecVal), portRC = sum(vecRC),
               vecVal = vecVal, vecRC = vecRC))
}

# ------------------------------------------------------------------------------
#' Age a Portfolio
#' 
#' @description Age a portfolio of VA policies specified in each inPolicy of inPortfolio from
#' currentDate (specified in inPolicy) to targetDate. The againg scenario is
#' given in fundScen. The time step length is specified in dT.
#' Here we input a rather irrelevant parameter df to "hack" for a more flexible
#' user-defined projection function.
#' @param inPortfolio A dataframe containing numPolicy rows and 45 attributes
#'   of each VA policy.
#' @param mortTable A dataframe with three columns of doubles representing the
#'   mortality table.
#' @param fundScen A numScen-by-numStep-by-numFund array of doubles of
#'   return factors (i.e., exp(mu_t dt)) in each period.
#' @param scenDates A vector containing strings in the format of "YYYY-MM-DD"
#'   of dates corresponding to each period in fundScen.
#' @param dT A double of stepsize in years; dT = 1 / 12 would be monthly.
#' @param targetDate A string in the format of "YYYY-MM-DD" of valuation date of
#'   the portfolio.
#' @param df A vector of doubles of risk-free discount rates of different tenor
#'   (not forward rates), should have length being numStep.
#' @return Outputs a dataframe containing numPolicy rows and 45 attributes of
#'   each VA policy, where currentDate, gbAmt, GMWBbalance, withdrawal,
#'   & fundValue of each policy could be updated as a result of aging.
#' @examples
#' targetDate <- "2016-01-01"
#' histFundScen <- genFundScen(fundMap, histIdxScen)
#' agePortfolio(VAPort[1:2, ], mortTable, histFundScen, histDates, dT = 1 / 12,
#' targetDate, cForwardCurve)
#' \dontrun{
#' targetDate <- "2001-01-01"
#' histFundScen <- genFundScen(fundMap, histIdxScen)
#' agePortfolio(VAPort, mortTable, histFundScen, histDates, dT = 1 / 12,
#' targetDate, cForwardCurve)
#' }
#' \dontrun{
#' VAPort[1, c("currentDate", "issueDate")] <- c("2001-01-01", "2001-01-01")
#' histFundScen <- genFundScen(fundMap, histIdxScen)
#' agePortfolio(VAPort, mortTable, histFundScen, histDates, dT = 1 / 12,
#' targetDate, cForwardCurve)
#' }
#' @section Note:
#' Target date MUST be PRIOR to the last date of historical scenario date,
#' Current date MUST be LATER than the first date of historical scenario date.
#' @export
agePortfolio <- function(inPortfolio, mortTable, fundScen, scenDates,
                         dT = 1 / 12, targetDate, df){
  numPolicy <- nrow(inPortfolio)
  for (i in 1:numPolicy){
    inPolicy <- inPortfolio[i, ]
    inPortfolio[i, ] <- ageOnePolicy(inPolicy, mortTable, fundScen, scenDates,
                                                 dT = 1 / 12, targetDate, df)
  }
  return (outPortfolio = inPortfolio)
}
