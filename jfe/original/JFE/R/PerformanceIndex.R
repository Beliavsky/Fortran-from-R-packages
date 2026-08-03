ActivePremium<-function (Ra, Rb, scale = NA) {
  Ra = xts::as.xts(Ra)
  Rb = xts::as.xts(Rb)
  if (is.na(scale)) {
    freq = xts::periodicity(Ra)
    if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
    } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
    } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
  }

  .activePremium <- function(Ra, Rb, scale) {
    merged = na.omit(merge(Ra, Rb))
    ap = Return.annualized(merged[, 1], scale = scale) - Return.annualized(merged[, 2], scale = scale)
    ap
  }

  result = sapply(Ra, .activePremium, Rb = Rb, scale = scale)
  if (length(result) == 1)
    return(result)
  else {
    dim(result) = c(ncol(Ra), ncol(Rb))
    colnames(result) = paste0("Active Premium:", colnames(Rb))
    rownames(result) = colnames(Ra)
    return(t(result))
  }
}

AdjustedSharpeRatio <-function(R,Rf = 0,FUN="StdDev"){
  R=xts::as.xts(R)
  freq = xts::periodicity(R)
  if (freq$scale %in% c("minute","hourly")) {stop("Data frequency too high")
  } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
  } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}

  SR <- SharpeRatio.annualized(R,Rf,FUN=FUN)
  results <-SR*(1+SR*(.Skewness(R)/6)-(SR^2)*.Kurtosis(R)/24)

  if (FUN=="VaR") {  rownames(results) = paste0("Adjusted VaR Sharpe(Risk = ", Rf, ")")
  } else if (FUN=="ES") {  rownames(results) = paste0("Adjusted ES Sharpe(Risk = ", Rf, ")")
  } else {  rownames(results) = paste0("Adjusted StdDev Sharpe(Risk = ", Rf, ")")}

  return(results)
}

BernardoLedoitRatio<-function (R) {
  R <- xts::as.xts(R)

  .bernardoLedoitRatio<-function(R) {
    R = na.omit(R)
    r1 = R[which(R > 0)]
    r2 = R[which(R < 0)]
    result = sum(r1)/-sum(r2)
    return(result)
  }

  results = sapply(R, .bernardoLedoitRatio)
  results <- t(results)
  colnames(results) = colnames(R)
  rownames(results) = "Bernardo and Ledoit ratio"
  return(results)

}



BurkeRatio <- function (R, Rf = 0, modified = FALSE) {
  drawdown = c()
  R0 <- R

  .burkeRatio<-function(R, Rf = 0, modified = modified) {
    number_drawdown = 0
    in_drawdown = FALSE
    peak = 1
    period = .Frequency(R)
    R=as.numeric(R)
    n = length(R)

    for (i in (2:length(R))) {
      if (R[i] < 0) {
        if (!in_drawdown) {
          peak = i - 1
          number_drawdown = number_drawdown + 1
          in_drawdown = TRUE
        }
      }
      else {
        if (in_drawdown) {
          temp = 1
          boundary1 = peak + 1
          boundary2 = i - 1
          for (j in (boundary1:boundary2)) {
            temp = temp * (1 + R[j] * 0.01)
          }
          drawdown = c(drawdown, (temp - 1) * 100)
          in_drawdown = FALSE
        }
      }
    }
    if (in_drawdown) {
      temp = 1
      boundary1 = peak + 1
      boundary2 = i
      for (j in (boundary1:boundary2)) {
        temp = temp * (1 + R[j] * 0.01)
      }
      drawdown = c(drawdown, (temp - 1) * 100)
      in_drawdown = FALSE
    }
    if (modified) {result = result * sqrt(n)
    } else {
      Rp = (prod(1 + R))^(period/length(R)) - 1
      result = (Rp - Rf)/sqrt(sum(drawdown^2))}

    return(result)
  }


  R = na.omit(xts::as.xts(R))
  result = apply(R, MARGIN = 2, .burkeRatio, Rf = Rf, modified = modified)
  result <- t(result)
  colnames(result) = colnames(R)
  if (modified) {
    rownames(result) = paste0("Modified Burke ratio (Risk free = ",Rf, ")")
  }
  else {
    rownames(result) = paste0("Burke ratio (Risk free = ",Rf, ")")
  }
  return(result)

}


DRatio <- function (R) {

  .dRatio<-function(R) {
    R = as.numeric(R)
    r1 = R[which(R > 0)]
    r2 = R[which(R < 0)]
    nd = length(r2)
    nu = length(r1)
    result = (-nd * sum(r2))/(nu * sum(r1))
    return(result)
  }

  R = xts::as.xts(na.omit(R))
  result = apply(R, MARGIN=2, .dRatio)
  result <- t(result)
  colnames(result) = colnames(R)
  rownames(result) = "D ratio"
  return(result)

}


KellyRatio<- function (R, Rf = 0) {


  .kellyRatio <- function(R, Rf) {
    xR = R-Rf
    KR = mean(xR, na.rm = TRUE)/sd(R, na.rm = TRUE)^2
    KR = KR/2
    return(KR)
  }
  R = xts::as.xts(R)
  result = sapply(R, .kellyRatio, Rf = Rf)
  dim(result) = c(1, NCOL(R))
  colnames(result) = colnames(R)
  rownames(result) = "Kelly Ratio"
  return(result)
}

#-------------------------------------------------------------




MartinRatio<-function (R, Rf = 0) {

  .martinRatio <- function(R, Rf = 0) {
    Period = .Frequency(R)
    UI = UlcerIndex(R)
    R=as.numeric(R)
    n = length(R)
    Rp = (prod(1 + R))^(Period/length(R)) - 1
    result = (Rp - Rf)/UI
  }


  R = xts::as.xts(na.omit(R))
  result = sapply(R, .martinRatio, Rf = Rf)
  result <- t(result)
  colnames(result) = colnames(R)
  rownames(result) = paste0("Martin Ratio (Rf = ", Rf, ")")
  return(result)

}



SkewnessKurtosisRatio <- function (R)
{

  .skewnessKurtosisRatio<-function(R) {
    numerator=.Skewness(R, method = "moment")
    denominator=.Kurtosis(R,method = "moment")
    result = numerator/denominator
    return(result)
  }
  R = na.omit(xts::as.xts(R))

  results = sapply(R,  .skewnessKurtosisRatio)
  results <- t(results)
  colnames(results) = colnames(R)
  rownames(results) = "SkewnessKurtosisRatio"
  return(results)

}

PainIndex<- function (R) {

  .pi <- function(R) {
    result = sum(abs(DrawdownPeak(R)))/nrow(R)
    return(result)
  }

  R = xts::as.xts(R)
  results = sapply(R, .pi)
  return(results)
}



MeanAbsoluteDeviation<-function (R) {

  .meanAbsoluteDeviation <- function(R) {
    R=as.matrix(R)
    result = sum(abs(R - mean(R)))/nrow(R)
    return(result)
  }
  R=xts::as.xts(R)
  results = apply(R, MARGIN = 2, .meanAbsoluteDeviation)
  results <- t(results)
  colnames(results) = colnames(R)
  rownames(results) = "Mean absolute deviation"
  return(results)

}

# ------------------------------------------------------------------------------
CalmarRatio <- function (R, scale = NA) {
  R = xts::as.xts(R)
  if (is.na(scale)) {
    freq = xts::periodicity(R)
    if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
    } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
    } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
  }
  result = Return.annualized(R, scale = scale)/abs(maxDrawdown(R))
  rownames(result) = "Calmar Ratio"
  return(result)
}


SterlingRatio <-function (R, scale = NA, excess = 0.1) {
  if (is.na(scale)) {
    freq = xts::periodicity(R)
    if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
    } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
    } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
  }
  #  R = xts::as.xts(R)
  annualized_return = Return.annualized(R, scale = scale)
  drawdown = abs(maxDrawdown(R) + excess)
  result = annualized_return/drawdown
  rownames(result) = paste0("Sterling Ratio (Excess = ", round(excess *100, 0), "%)")
  return(result)
}




# ------------------------------------------------------------------------------
AppraisalRatio <- function (Ra, Rb, Rf = 0, method = c("appraisal", "modified",
                                                       "alternative")) {
  method = method[1]

  .appraisalRatio <-function(Ra, Rb, Rf = 0, method = method) {
	Period=.Frequency(Ra)
	numerator=CAPM.jensenAlpha(Ra, Rb, Rf)

    switch(method, appraisal = {
      epsilon = residuals(lm(Ra~Rb))
      specifikRisk = sqrt(sum((epsilon - mean(epsilon))^2)/length(epsilon))*sqrt(Period)
      result = numerator/specifikRisk
    }, modified = {
      y=Ra-Rf
      x=Rb-Rf
      Beta=coef(lm(y~x))[2]
      result = numerator/Beta
    }, alternative = {
      y=Ra-Rf
      x=Rb-Rf
      Beta=coef(lm(y~x))[2]
      sRisk=Beta*(sd(Rb)*sqrt(Period))
      result = numerator/sRisk
    })

  }
  Ra = xts::as.xts(Ra)
  Rb = xts::as.xts(Rb)

  result = sapply(Ra, .appraisalRatio, Rb = Rb, Rf = Rf, method = method)
  result <- t(result)
  colnames(result) = colnames(Ra)
  switch(method, appraisal = {
    rownames(result) = paste0("Appraisal ratio (Risk free = ",Rf, ")")
  }, modified = {
    rownames(result) = paste0("Modified Jensen's alpha (Risk free = ",Rf, ")")
  }, alternative = {
    rownames(result) = paste0("Alternative Jensen's alpha (Risk free = ",Rf, ")")
  })
  return(result)

}


TrackingError<- function (Ra, Rb, scale = NA) {
  if (is.na(scale)) {
    freq = xts::periodicity(Ra)
    if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
    } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
    } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
  }

  .trackingError <- function(Ra,Rb,scale) {
    result=sd(Ra-Rb,na.rm = TRUE)* sqrt(scale)
  }

  Ra=xts::as.xts(Ra)
  Rb=xts::as.xts(Rb)
  result=sapply(Ra, .trackingError, Rb,scale=scale)

  if (length(result) == 1)
    return(result)
  else {
    dim(result) = c(ncol(Ra), ncol(Rb))
    colnames(result) = paste0("Tracking Error:", colnames(Rb))
    rownames(result) = colnames(Ra)
    return(t(result))
  }
}

InformationRatio<- function (Ra, Rb, scale = NA) {
  Ra = xts::as.xts(Ra)
  Rb = xts::as.xts(Rb)

  if (is.na(scale)) {
    freq = xts::periodicity(Ra)
    if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
    } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
    } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
  }

  .informationRatio <- function(Ra, Rb, scale) {
    ap = ActivePremium(Ra, Rb, scale = scale)
    te = TrackingError(Ra, Rb, scale = scale)
    IR = ap/te
    return(IR)
  }
  result = sapply(Ra, .informationRatio, Rb = Rb, scale = scale)

  if (length(result) == 1)
    return(result)
  else {
    result = matrix(result, ncol = ncol(Ra), nrow = ncol(Rb),byrow = TRUE)
    rownames(result) = paste("Information Ratio:", colnames(Rb))
    colnames(result) = colnames(Ra)
    return(result)
  }
}


TreynorRatio <- function (Ra, Rb, Rf = 0, scale = NA, modified = FALSE) {
  Ra = xts::as.xts(Ra)
  Rb = xts::as.xts(Rb)
  if (is.na(scale)) {
    freq = xts::periodicity(Ra)
    if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
    } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
    } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
  }
  Ra.ncols = NCOL(Ra)
  Rb.ncols = NCOL(Rb)
  pairs = expand.grid(1:Ra.ncols, 1:Rb.ncols)
  xRa = Ra-Rf
  xRb = Rb-Rf

  .treynorRatio <- function(xRa, xRb, scale) {
    beta = coef(lm(xRa~xRb))[2]
    TR = Return.annualized(xRa, scale = scale)/beta
    TR
  }

  .mTreynorRatio <- function(xRa, xRb, scale) {
    Period = .Frequency(xRa)
    Beta=coef(lm(xRa~xRb))[2]
    sRisk=Beta*(sd(xRb)*sqrt(Period))
    mTR = Return.annualized(xRa, scale = scale)/sRisk
    mTR
  }
  if (modified) {result = sapply(xRa, .mTreynorRatio, xRb = xRb, scale = scale)}
  else {result = sapply(xRa, .treynorRatio, xRb = xRb, scale = scale)}


  if (length(result) == 1)
    return(result)
  else {
    dim(result) = c(Ra.ncols, Rb.ncols)
    colnames(result) = paste0("Treynor Ratio:", colnames(Rb))
    rownames(result) = colnames(Ra)
    return(t(result))
  }
}


DownsideDeviation <-function (R, MAR = 0, method = c("full","subset"), potential = FALSE)
{

  method = method[1]

  .downsideDeviation <- function (R,MAR,method) {
    r = subset(R, R < MAR)
    if (!is.null(dim(MAR))) {
      if (xts::is.timeBased(zoo::index(MAR))) {
        MAR <- MAR[zoo::index(r)]
      }
      else {
        MAR = mean(MAR)
      }
    }
    switch(method, full = {
      len = length(R)
    }, subset = {
      len = length(r)
    })
    if (potential) {
      result = sum((MAR - r)/len)
    }
    else {
      result = sqrt(sum((MAR - r)^2/len))
    }
    return(result)

  }


  R = xts::as.xts(na.omit(R))
  results = apply(R, MARGIN=2, .downsideDeviation, MAR = MAR, method = method)
  results <- t(results)
  colnames(results) = colnames(R)
  if (potential)
    rownames(results) = paste0("Downside Potential (MAR = ",round(mean(MAR), 1), "%)")
  else rownames(results) = paste0("Downside Deviation (MAR = ",round(mean(MAR), 1), "%)")
  return(results)

}


OmegaSharpeRatio <-function (R, MAR = 0) {

  .omegaSharpeRatio <- function(R, MAR = 0) {

    R = na.omit(R)
    r = R[which(R > MAR)]

    if (xts::is.timeBased(zoo::index(MAR))) {
      MAR <- MAR[zoo::index(r)]
    }
    else {
      MAR = mean(MAR)
    }

    result = (.UpsideRisk(R, MAR, stat = "potential") -
                .DownsidePotential(R, MAR))/(.DownsidePotential(R,MAR))

    return(result)
  }

  R = xts::as.xts(R)
  results = sapply(R, .omegaSharpeRatio, MAR = MAR)
  results <- t(results)
  colnames(results) = colnames(R)
  rownames(results) = paste0("OmegaSharpeRatio (MAR = ", MAR, "%)")
  return(results)

}

SortinoRatio<- function (R, MAR = 0) {

  .sortinoRatio <- function(R, MAR) {
    SR = mean((R-MAR), na.rm = TRUE)/DownsideDeviation(R,MAR)
    SR
  }
  R = xts::as.xts(R)
  result = sapply(R, .sortinoRatio, MAR = MAR)
  dim(result) = c(1, NCOL(R))
  colnames(result) = colnames(R)
  rownames(result) = paste0("Sortino Ratio (MAR = ", round(mean(MAR) *100, 3), "%)")
  return(result)
}


ProspectRatio <- function (R, MAR) {

  .prospectRatio<-function(R, MAR) {
    n = nrow(R)
    SigD = DownsideDeviation(R, MAR)
    r1 = R[which(R > 0)]
    r2 = R[which(R < 0)]
    result = (sum(r1) + 2.25 * sum(r2) - MAR)/(SigD *n)
    return(result)
  }

  R <- xts::as.xts(R)
  result = sapply(R, .prospectRatio, MAR = MAR)
  result <- t(result)
  colnames(result) = colnames(R)
  rownames(result) = paste0("Prospect ratio (MAR = ", MAR,"%)")
  return(result)

}


VolatilitySkewness<- function (R, MAR = 0, stat = c("volatility", "variability")) {

  stat = stat[1]

  .volatilitySkewness<-function(R, MAR = 0,stat=stat) {

    if (!is.null(dim(MAR))) {
      if (xts::is.timeBased(zoo::index(MAR))) {
        MAR <- MAR[zoo::index(R)]
      }
      else {
        MAR = mean(MAR)
      }
    }

    switch(stat, volatility = {
      result = .UpsideRisk(R, MAR, stat = "variance")/DownsideDeviation(R,MAR)^2
    }, variability = {
      result = .UpsideRisk(R, MAR, stat = "risk")/DownsideDeviation(R, MAR)
    }, )
    return(result)
  }

  R = xts::as.xts(R)
  result = sapply(R, .volatilitySkewness, MAR = MAR,stat = stat)
  result <- t(result)
  colnames(result) = colnames(R)
  rownames(result) = paste0("VolatilitySkewness (MAR = ",MAR, "%, stat= ", stat, ")")
  return(result)

}
# ------------------------------------------------------------------------------
M2Sortino <- function (Ra, Rb, MAR = 0) {


  .m2Sortino<-function(Ra, Rb, MAR = MAR)  {

    #Period = .Frequency(Rb)
    Rp = (prod(1 + Ra))^(.Frequency(Rb)/nrow(Ra)) - 1
    SigmaD = DownsideDeviation(Ra, MAR) * sqrt(.Frequency(Rb))
    SigmaDM = DownsideDeviation(Rb, MAR) * sqrt(.Frequency(Rb))
    SR = SortinoRatio(Ra, MAR)
    result = Rp + SR * (SigmaDM - SigmaD)

    return(result)
  }
  Ra = xts::as.xts(Ra)
  Rb = xts::as.xts(Rb)

  results = sapply(Ra, .m2Sortino, Rb = Rb, MAR = MAR)
  results <- t(results)
  colnames(results) = colnames(Ra)
  rownames(results) = paste0("M2Sortino (MAR = ", MAR, ")")
  return(results)

}


# ------------------------------------------------------------------------------

SharpeRatio <- function(R, Rf = 0, alpha=0.05,FUN="StdDev", annualize=FALSE) {

  R=xts::as.xts(R)
  if (annualize) {

    SRm =SharpeRatio.annualized(R, Rf,scale = scale, FUN)

  } else {

    if (FUN =="StdDev") {

      SRm = sapply(R,mean)/sapply(R,.covRisk)
      SRm=t(as.matrix(SRm))
      rownames(SRm)=paste0("StdDev Sharpe(Rf=", Rf, "%):")

    } else if (FUN =="VaR") {

      SRm= sapply(R,mean)/-sapply(R,.varRisk, alpha)
      SRm=t(as.matrix(SRm))
      rownames(SRm)=paste0("VaR Sharpe(Rf=", Rf, "%,alpha=",alpha,")")

    } else if (FUN =="ES") {

      SRm= sapply(R,mean)/-sapply(R,.cvarRisk,alpha)
      SRm=t(as.matrix(SRm))
      rownames(SRm)=paste0("ES Sharpe(Rf=", Rf, "%,alpha=",alpha,")")}
  }

  results=SRm
  return(results)
}


SharpeRatio.annualized <-function(R, Rf = 0, alpha=0.05,scale = NA, geometric = TRUE,FUN="StdDev"){

  R=xts::as.xts(R)
  freq = xts::periodicity(R)
  if (freq$scale %in% c("minute","hourly")) {stop("Data frequency too high")
  } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
  } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}


  if (FUN =="StdDev") {
    SRm = (Return.annualized(R, scale = scale, geometric = geometric)-Rf)/( timeSeries::colSds(R)*sqrt(scale))
    SRm=t(as.matrix(SRm))
    rownames(SRm)=paste0("Annualized StdDev Sharpe(Rf=", Rf, "%):")
  }  else if (FUN =="VaR") {
    SRm= (Return.annualized(R, scale = scale)-Rf)/(-sapply(R,.varRisk,alpha))
    SRm=t(as.matrix(SRm))
    rownames(SRm)=paste0("Annualized VaR Sharpe(Rf=", Rf, "%,alpha=",alpha,")")
  }  else if (FUN =="ES") {
    SRm= (Return.annualized(R, scale = scale)-Rf)/(-sapply(R,.cvarRisk,alpha))
    SRm=t(as.matrix(SRm))
    rownames(SRm)=paste0("Annualized ES Sharpe(Rf=", Rf, "%,alpha=",alpha,")")}

  results <- SRm
  return(results)
}


PainRatio <- function (R, Rf = 0) {


  .painRatio <- function(R, Rf=Rf) {

    Period = .Frequency(R)
    PI = PainIndex(R)
    Rp = (prod(1 + R))^(Period/nrow(R)) - 1
    result = (Rp - Rf)/PI

    return(result)
  }
  R=xts::as.xts(R)
  results = sapply(R,  .painRatio, Rf = Rf)
  results <- t(results)
  colnames(results) = colnames(R)
  rownames(results) = paste0("Pain Ratio (Rf = ", Rf, ")")
  return(results)

}




table.AnnualizedReturns <-function (R, scale = NA, Rf = 0, geometric = TRUE, digits = 4)
{
  y = xts::as.xts(R)
  #  if (!is.null(dim(Rf)))
  columns = ncol(y)
  columnnames = colnames(y)

  if (is.na(scale)) {
    freq = xts::periodicity(R)
    if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
    } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="weekly") {scale=52
    } else if (freq$scale=="monthly") {scale=12} else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
  }

  for (column in 1:columns) {
    z = c(Return.annualized(y[, column, drop = FALSE], scale = scale, geometric = geometric)
          ,sd(y[, column])*sqrt(scale),SharpeRatio.annualized(y[,column, drop = FALSE], scale = scale, Rf = Rf, geometric = geometric))
    znames = c("Annualized Return", "Annualized Std Dev",
               paste0("Annualized Sharpe (Rf=", base::round(mean(Rf) *scale, 4) * 100, "%)"))

    if (column == 1) {
      resultingtable = data.frame(Value = z, row.names = znames)
    }
    else {
      nextcolumn = data.frame(Value = z, row.names = znames)
      resultingtable = cbind(resultingtable, nextcolumn)
    }
  }
  colnames(resultingtable) = columnnames
  ans = base::round(resultingtable, digits)
  ans
}


Return.annualized <-function (R, scale = NA, geometric = TRUE) {
  R=xts::as.xts(R)

  .return.annualized <- function(R, scale = scale, geometric = geometric) {
    R = na.omit(R)
    n = nrow(R)

    if (is.na(scale)) {
      freq = xts::periodicity(R)
      if (freq$scale %in% c("minute","hourly")) {stop("Data frequency is too high")
      } else if (freq$scale=="daily") {scale=252} else if (freq$scale=="monthly") {scale=12} else if (freq$scale=="weekly") {scale=52
      } else if (freq$scale=="quarterly") {scale=4} else if (freq$scale=="yearly") {scale=1}
    }

    if (geometric) {
      result = prod(1 + R)^(scale/n) - 1
    }
    else {
      result = mean(R) * scale
    }
    result
  }
  results=sapply(R, .return.annualized, scale = scale, geometric = geometric)

  return(results)
}


CAPM.jensenAlpha<- function (Ra, Rb, Rf = 0) {

  Ra = xts::as.xts(Ra)
  Rb = xts::as.xts(Rb)

  .capmJensenAlpha <-function(Ra, Rb = Rb,Rf = Rf) {
    #Period = .Frequency(Ra)
    Rp = (prod(1 + Ra))^(.Frequency(Ra)/length(Ra)) - 1
    Rpb = (prod(1 + Rb))^(.Frequency(Ra)/length(Rb)) - 1
    y=Ra-Rf
    x=Rb-Rf
    Beta=coef(lm(y~x))[2]
    result = Rp - Rf - Beta * (Rpb -Rf)
  }

  result = sapply(Ra, .capmJensenAlpha, Rb = Rb,Rf = Rf)
  result <- t(result)
  colnames(result) = colnames(Ra)
  rownames(result) = paste0("Jensen's Alpha (Risk free = ",Rf, ")")
  return(result)

}




UlcerIndex <-function (R) {

  .ui <- function(R) {
    result = sqrt(sum(DrawdownPeak(R)^2))
    result = result/sqrt(nrow(R))
    return(result)
  }

  R = xts::as.xts(R)
  results = sapply(R, .ui)
  dim(results) = c(1, NCOL(R))
  colnames(results) = colnames(R)
  rownames(results) = "Ulcer Index"
  return(results)
}



DrawdownPeak <- function (R) {

  .drawdownpeak <- function(R) {
    R=as.numeric(R)
    drawdownpeak = rep(NA,length(R))
    peak = 0
    for(i in (1:length(R))) {
      val = 1
      borne = peak+1
      for(j in (borne:i)) {
        val = val*(1+R[j]/100)
      }
      if (val > 1) {
        peak = i
        drawdownpeak[i] = 0
      } else {
        drawdownpeak[i] = (val-1)*100
      }
    }
    drawdownpeak

  }

  R = xts::as.xts(R)
  results = sapply(R, .drawdownpeak)
  colnames(results) = colnames(R)
  return(results)

}


.DownsidePotential <- function (R, MAR = 0) {

  .downsidePotential<-function(R, MAR = 0) {
    return(DownsideDeviation(R, MAR = MAR, method = "full", potential = TRUE))
  }

  R = xts::as.xts(R)
  results = sapply(R, .downsidePotential, MAR = MAR)
  results = matrix(results, nrow = 1)
  colnames(results) = colnames(R)
  rownames(results) = paste0("Downside Potential (MAR = ",round(mean(MAR), 1), "%)")
  return(results)

}


.UpsideRisk<-  function (R, MAR = 0, method = c("full", "subset"), stat = c("risk",
                                                                           "variance", "potential"))
{
  method = method[1]
  stat = stat[1]

  .upsideRisk <- function(R, MAR = 0,method,stat) {
    R = na.omit(R)
    r = R[which(R > MAR)]
    if (!is.null(dim(MAR))) {
      if (xts::is.timeBased(zoo::index(MAR))) {
        MAR <- MAR[zoo::index(r)]
      }
      else {
        MAR = mean(MAR)
      }
    }
    switch(method, full = {
      len = length(R)
    }, subset = {
      len = length(r)
    })
    switch(stat, risk = {
      result = sqrt(sum((r - MAR)^2/len))
    }, variance = {
      result = sum((r - MAR)^2/len)
    }, potential = {
      result = sum((r - MAR)/len)
    })
    return(result)
  }

  R = xts::as.xts(R)
  results = sapply(R, .upsideRisk, MAR = MAR, method = method, stat = stat)
  results <- t(results)
  colnames(results) = colnames(R)
  rownames(results) = paste0("Upside Risk (MAR = ", MAR,"%)")
  return(results)

}



.Drawdowns <- function (R, geometric = TRUE) {
  x = xts::as.xts(R)
  columns = ncol(x)
  columnnames = colnames(x)

  .colDrawdown <- function(x, geometric, ...) {
    if (geometric) { Return.cumulative = cumprod(1 + x)
    } else {Return.cumulative = 1 + cumsum(x)}

    maxCumulativeReturn = cummax(c(1, Return.cumulative))[-1]
    column.drawdown = Return.cumulative/maxCumulativeReturn - 1
  }

  for (column in 1:columns) {
    column.drawdown <- .na.skip(x[, column], FUN = .colDrawdown, geometric = geometric)
    if (column == 1)
      drawdown = column.drawdown
    else drawdown = merge(drawdown, column.drawdown)
  }

  colnames(drawdown) = columnnames
  #  drawdown = reclass(drawdown, x)
  drawdown = xts::as.xts(drawdown, x)
  return(drawdown)
}


.na.skip <- function (x, FUN=NULL, ...) {
  nx <- na.omit(x)
  fx <- FUN(nx, ... = ...)
  if (is.vector(fx)) {
    result <- xts::.xts(fx, xts::.index(x), .indexCLASS = xts::indexClass(x))
  }
  else {
    result <- merge(fx, xts::.xts(, xts::.index(x)))
  }
  return(result)
}

maxDrawdown <-function (R, geometric = TRUE, invert = TRUE) {

  .maxdrawdown<-function(R, geometric = TRUE, invert = TRUE) {
    drawdown = .Drawdowns(R, geometric = geometric)
    result = min(drawdown)
    if (invert)
      result <- -result
    return(result)
  }

  R = na.omit(xts::as.xts(R))
  result = sapply(R, .maxdrawdown, geometric = geometric,invert = invert)
  dim(result) = c(1, NCOL(R))
  colnames(result) = colnames(R)
  rownames(result) = "Worst Drawdown"
  return(result)

}

.Kurtosis<- function (R, na.rm = FALSE, method = c("excess", "moment", "fisher",
                                                   "sample", "sampleExcess"))
{
  method = method[1]

  .kurtosis<-function(x,moment) {
    x=as.numeric(x)
    n=length(x)
    if (method == "excess") {
      kurtosis = sum((x - mean(x))^4/(var(x) * (n - 1)/n)^2)/length(x) -
        3
    }
    if (method == "moment") {
      kurtosis = sum((x - mean(x))^4/(var(x) * (n - 1)/n)^2)/length(x)
    }
    if (method == "fisher") {
      kurtosis = ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 -
                                         (3 * (n - 1))/(n + 1)))/((n - 2) * (n - 3))
    }
    if (method == "sample") {
      kurtosis = sum((x - mean(x))^4/var(x)^2) * n * (n +
                                                        1)/((n - 1) * (n - 2) * (n - 3))
    }
    if (method == "sampleExcess") {
      kurtosis = sum((x - mean(x))^4/var(x)^2) * n * (n +
                                                        1)/((n - 1) * (n - 2) * (n - 3)) - 3 * (n - 1)^2/((n -
                                                                                                             2) * (n - 3))
    }
    kurtosis

  }
  R = xts::as.xts(R)
  result=sapply(R, .kurtosis,method)
  result=t(as.matrix(result))
  colnames(result) = colnames(R)
  rownames(result) = "Excess Kurtosis"

  result
}






.Skewness <-function (R, method = c("moment", "fisher", "sample")){
  method = method[1]

  .skewness<-function(x, method) {
    x=as.numeric(x)
    n = length(x)
    if (method == "moment") {
      Skewness = sum((x - mean(x))^3/sqrt(var(x) * (n -1)/n)^3)/n
    } else  if (method == "fisher") {
      if (n < 3) {Skewness = NA
      } else {Skewness = ((sqrt(n * (n - 1))/(n - 2)) * (sum(x^3)/n))/((sum(x^2)/n)^(3/2))}
    }
    else if (method == "sample") {
      Skewness = sum((x - mean(x))^3/sqrt(var(x) * (n -1)/n)^3) * n/((n - 1) * (n - 2))
    }

    Skewness
  }
  R = xts::as.xts(na.omit(R))
  result=sapply(R,  .skewness, method=method)
  result=t(as.matrix(result))
  colnames(result) = colnames(R)
  rownames(result) = "Skewness"

  result
}






# ------------------------------------------------------------------------------

.covRisk <-  function(data)  {

    x.mat <- as.matrix(data)
    # Risk:
    Std = sd(x.mat)
    names(Std) = "Cov"

    # Return Value:
    Std

  }


.varRisk <-  function(data, alpha = 0.05)  {
    # VaR:

    x.mat = as.matrix(data)
    VaR = quantile(x.mat, alpha, type = 1)
    names(VaR) <- paste0("VaR.", alpha*100, "%")

    # Return Value:
    VaR
  }

.cvarRisk <-  function(data, alpha = 0.05)  {

    # CVaR:
  x.mat = as.matrix(data)
    VaR = quantile(x.mat, alpha, type = 1)
    CVaR = c(CVaR = VaR - 0.5 * mean(((VaR-x.mat) + abs(VaR-x.mat))) / alpha)
    names(CVaR) <- paste0("CVaR.", alpha*100, "%")

    # Return Value:
    CVaR
  }

.Frequency <- function (R) {
  R = xts::as.xts(R)

  .frequency <- function(R) {
      freq = xts::periodicity(R)
      switch(freq$scale, minute = {
        stop("Data Frequency is too high")
      }, hourly = {
        stop("Data Frequency is too high")
      }, daily = {
        result = 252
      }, weekly = {
        result = 52
      }, monthly = {
        result = 12
      }, quarterly = {
        result = 4
      }, yearly = {
        result = 1
      })
      return(result)
    }

    results = sapply(R, .frequency)
    results <- t(results)
    colnames(results) = colnames(R)
    rownames(results) = "Frequency"
    return(results)

}
