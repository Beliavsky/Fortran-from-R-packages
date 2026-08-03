Monthly2Weekly <- function(Data, Economy, EcoVar) {
  Tmonths <- nrow(Data)

  # Replace dates to be the last one of the month
  FirstDayMonth <- as.Date(Data$Period)
  LastDayMonth <- lubridate::ceiling_date(FirstDayMonth, "month") - 1
  Data$Period <- LastDayMonth

  # Find daily interpolated series
  Date <- seq(Data$Period[1], Data$Period[Tmonths], by = "1 day")
  d <- length(Date)
  DaysYear <- 366 # since 2020 is a leap year

  DataDaily <- data.frame(Date)
  DataDaily[[EcoVar]] <- spline(Data$Period, Data[[Economy]], n = d)$y

  # (log) YoY growth
  X <- c()
  for (i in 1:(d - DaysYear)) {
    X[i] <- log(DataDaily[[EcoVar]][i + DaysYear]) - log(DataDaily[[EcoVar]][i])
  }

  DataDaily$YoY <- c(rep(NA, times = DaysYear), X)
  VarLabel <- paste(EcoVar, "_", Economy, sep = "")
  YoYLabel <- paste("YoY_", VarLabel, sep = "")
  names(DataDaily) <- c("Date", VarLabel, YoYLabel)


  # Generate weekly series
  WeeklySeries <- Daily2WeeklyCOVIDEveryDay(DataDaily)

  return(WeeklySeries)
}

#######################################################################################################
Daily2WeeklyCOVIDEveryDay <- function(SeriesOfInterest) {
  # GOAL: creates end of week series for COVID times for series containing only business days (market related)
  DatesWeekSundays <- seq(from = as.Date("2020-01-05"), to = as.Date("2022-01-16"), by = "7 days")

  # Uniform dates
  DatesWeekSundays <- as.Date(DatesWeekSundays)
  SeriesOfInterestDates <- as.Date(SeriesOfInterest$Date)

  # Find the indexes of Friday weekdays
  IdxDates <- match(DatesWeekSundays, SeriesOfInterestDates)
  checkNAs <- which(is.na(IdxDates == 1))

  # Generate clean VIX time series
  SeriesOfInterestweekly <- SeriesOfInterest[IdxDates, ]
  SeriesOfInterestweekly$Date <- as.Date(DatesWeekSundays)

  return(SeriesOfInterestweekly)
}
