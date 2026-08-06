# ------------------------------------------------------------------------------
# ----- Step1_YieldCurveGeneration.R -------------------------------------------
# ------------------------------------------------------------------------------

checkDCC <- function(dcc){
  if (!(dcc == "Thirty360" || dcc == "ACT360"
        || dcc == "ACT365" || dcc == "ACTACT")){
    errstr <- gsub("[\r\n]", "", 
                   sprintf("dcc should be one of Thirty360, ACT360, ACT365, or ACTACT. Given %s.", dcc))
    stop (errstr)
  }
}

checkBDC <- function(bdc){
  # Check if the business day convention bdc is a valid input ----------------
  if (bdc != "Actual" && bdc != "Following" && bdc != "Preceding"
      && bdc != "Modified_Prec" && bdc != "Modified_Foll") {
    errstr <- sprintf("Business day convention should be one of Actual, Following,
                          Preceding, Modified_Prec, or Modified_Foll.
                          Given %s.", bdc)
    stop (errstr)
  }
}
checkCalendar <- function(calendar){
  # Check if the holiday calendar is a valid input -----------------
  if (calendar != "General" && calendar != "NY"){
    errstr <- sprintf("calendar should be one of General or NY. Given %s.",
                      calendar)
    stop (errstr)
  }
}

checkHolidays <- function(holidays){
  # Check if the user-defined holidays is a valid input ------------------
  if ((!is.null(holidays)) && (!inherits(holidays, "Date")) ) {
    errstr <- sprintf("Holidays should be date type, given %s type.", typeof(holidays))
    stop (errstr)
  }
}
# ------------------------------------------------------------------------------
isLeapYear <- function(date){
  # Indicate if the given date is in a leap year  ----------------------------
  # Input  -------------------------------------------------------------------
  # -- date: string "YYYY-MM-DD" or a Date data type -------------------------
  # Output -------------------------------------------------------------------
  # -- Binary, TRUE if date is in a leap year or FALSE otherwise -------------
  
  date <- as.Date(date)
  year <- as.numeric(format(date, format = "%Y"))
  
  bres <- FALSE
  if (year %% 4 == 0 && year %% 100 != 0) {
    bres <- TRUE
  }
  if (year %% 400 == 0) {
    bres <- TRUE
  }
  
  return(bres)
}

# ------------------------------------------------------------------------------
leapYearDays <- function(dateA, dateB){
  # Calculate the number of days in leap years between dateA and dateB
  # Input
  # -- dateA: string "YYYY-MM-DD" or a Date data type
  # -- dateB: string "YYYY-MM-DD" or a Date data type
  # Output
  # -- Integer number of days in leap years between dateA and dateB
  
  dateA <- as.Date(dateA)
  dateB <- as.Date(dateB)
  
  if (dateA > dateB) {
    return(leapYearDays(dateB, dateA))
  }
  # Then we can assume datA is before datB
  yearA <- as.numeric(format(dateA, format = "%Y"))
  yearB <- as.numeric(format(dateB, format = "%Y"))
  
  days <- 0
  for (y in yearA:yearB) {
    dateTemp <- as.Date(paste(y, 1, 1, sep = "-"))
    if (isLeapYear(dateTemp)) {
      days <- days + as.Date(paste(y + 1, 1, 1, sep = "-")) - dateTemp
    }
  }
  
  if (isLeapYear(dateA)) {
    days <- days + (as.Date(paste(yearA - 1, 12, 31, sep = "-")) - dateA)
  }
  
  if (isLeapYear(dateB)) {
    days <- days + (dateB - as.Date(paste(yearB, 12, 31, sep = "-")))
  }
  
  return(as.numeric(days))
}

# ------------------------------------------------------------------------------
fracYear <- function(dateA, dateB, dcc){
  # Calculates the year fraction between dateA and dateB according to day
  # count convention dcc
  # Input
  # -- dateA: string "YYYY-MM-DD" or a Date data type
  # -- dateB: string "YYYY-MM-DD" or a Date data type
  # -- dcc: string, day count convention
  #    -- Thirty360: 30 days a month, 360 days a year
  #    -- ACT360: actual day count, 360 days a year
  #    -- ACT365: actual day count, 365 days a year
  #    -- ACTACT: actual day count, actual days in a year
  # Output
  # -- double, year fraction between two dates
  
  dateA <- as.Date(dateA)
  dateB <- as.Date(dateB)
  checkDCC(dcc)
  
  if (dateB < dateA) {
    return (fracYear(dateB, dateA, dcc))
  }
  
  yearA <- as.numeric(format(dateA, format = "%Y"))
  yearB <- as.numeric(format(dateB, format = "%Y"))
  monthA <- as.numeric(format(dateA, format = "%m"))
  monthB <- as.numeric(format(dateB, format = "%m"))
  dayA <- as.numeric(format(dateA, format = "%d"))
  dayB <- as.numeric(format(dateB, format = "%d"))
  
  if (dcc == "Thirty360") {
    dayA <- min(dayA, 30)
    if (dayB == 31 && dayA == 30) {
      dayB <- 30
    }
    frac <- (yearB - yearA) + (monthB - monthA) / 12 + (dayB - dayA) / 360
  }
  else if (dcc == "ACT360") {
    frac <- as.numeric((dateB - dateA)) / 360
  }
  else if (dcc == "ACT365") {
    frac <- as.numeric((dateB - dateA)) / 365
  }
  else if (dcc == "ACTACT") {
    leapDays <- leapYearDays(dateA, dateB)
    nonLeapDays <- abs(as.numeric(dateB - dateA)) - leapDays
    frac <- leapDays / 366 + nonLeapDays / 365
  }
  return (frac)
}

# ------------------------------------------------------------------------------
isWeekend <- function(date){
  # Indicate whether the given date is a weekend
  # Input:
  # -- date: string "YYYY-MM-DD" or a Date data type
  # Output:
  # -- binary output, TRUE if date is a weekend or FALSE otherwise
  
  date <- as.Date(date)
  w <- weekdays(date)
  
  if (w == "Saturday" || w == "Sunday"){
    return (TRUE)
  }else{
    return (FALSE)
  }
}

# ------------------------------------------------------------------------------
isHolidayNY <- function(date){
  # Indicate whether the given date is a holiday based on New York holiday
  # schedule
  # Input:
  # -- date: string "YYYY-MM-DD" or a Date data type
  # Output:
  # -- binary output, TRUE if holiday or FALSE otherwise
  
  date <- as.Date(date)
  w <- weekdays(date) # day of the week as string
  d <- as.numeric(format(date, format = "%d"))  # day in a month as numeric
  m <- as.numeric(format(date, format = "%m"))  # month in a year as numeric
  
  if (# MLK day, third Mon in Jan
    ((d >= 15 && d <= 21) && w == "Monday" && m == 1)
    # Washington's birthday, 3rd Mon in Feb
    || ((d >= 15 && d <= 21) && w == "Monday" && m == 2)
    # Mamorial day, last Mon in may
    || (d >= 25 && w == "Monday" && m == 5)
    # National day (Fri/Mon before/after)
    || ((d == 4 && m == 7) || (d == 5 && w == "Monday" && m == 7)
        || (d == 3 && w == "Friday" && m == 7))
    # Labor day, 1st Mon in Sep
    || (d <= 7 && w == "Monday" && m == 9)
    # Columbus day, 2nd Mon in Oct
    || ((d >= 8 && d <= 14) && w == "Monday" && m == 10)
    # Veteran's day
    || ((d == 11 && m == 11) || (d == 12 && w == "Monday" && m == 11)
        || (d == 10 && w == "Friday" && m == 11))
    # US thanksgiving, 4th Thu in Nov
    || ((d >= 22 && d <= 28) && w == "Thursday" && m == 11)
    # Christmas (Fri/Mon before/after)
    || ((d == 25 && m == 12) || (d == 26 && w == "Monday" && m == 12)
        || (d == 24 && w == "Friday" && m == 12))
    # New Year (Fri/Mon before/after)
    || ((d == 1 && m == 1) || (d == 2 && w == "Monday" && m == 1)
        || ((d == 31 && w == "Friday" && m == 12)))
  ){
    return (TRUE)
  } else {
    return (FALSE)
  }
}

# ------------------------------------------------------------------------------
isBusinessDay <- function(date, calendar, holidays){
  # Indicates whether date is a business day based on calendar
  # Input
  # -- date: string "YYYY-MM-DD" or a Date data type
  # -- calendar: string, indicate the desired calendar
  #    -- NY: New York holiday calendar
  #    -- General: all weekdays are business days
  # Output
  # -- Binary output, indicating whether the given date is a business day
  
  date <- as.Date(date)
  checkCalendar(calendar)
  
  if(!is.null(holidays) && date >= min(holidays) && date <= max(holidays)){
    return(date %in% holidays)
  }
  else if (calendar == "NY") {
    return(!isWeekend(date) && !isHolidayNY(date))
  } 
  else if (calendar == "General") {
    return(!isWeekend(date))
  }
}

# ------------------------------------------------------------------------------
rollDate <- function(date, bdc, calendar, holidays){
  # Roll the given date to the nearest business day based on a a given
  # business day convention bdc and holiday calendar calendar
  # Input
  # -- date: string "YYYY-MM-DD" or a Date data type
  # -- bdc: string, business day convention 
  #    -- Actual: No rolling on the date applied even if it is a non-business day
  #    -- Preceding: 1st business day before holiday
  #    -- Following: 1st business day after holiday
  #    -- Modified_Prec: Same as "Preceding" unless it belongs to a different
  #                      month, in which case 1st business day after holiday
  #    -- Modified_Foll: Same as "Following" unless it belongs to a different
  #                      month, in which case 1st business day before holiday
  # -- calendar: string, indicate the desired calendar
  #    -- NY: New York holiday calendar
  #    -- General: all weekdays are business days
  
  date <- as.Date(date)
  checkCalendar(calendar)
  checkBDC(bdc)
  
  dateTemp <- date
  while ( !bdc == "Actual" && !isBusinessDay(dateTemp, calendar, holidays)) {
    if (bdc == "Following") {
      dateTemp <- dateTemp + 1
    } else if (bdc == "Preceding") {
      dateTemp <- dateTemp - 1
    } else if (bdc == "Modified_Prec") {
      dateTemp <- rollDate(date, bdc = "Preceding", calendar, holidays)
      if (format(dateTemp, format = "%m")
          != format(date, format = "%m")){
        dateTemp <- rollDate(date, bdc = "Following", calendar, holidays)
      }
    }else if (bdc == "Modified_Foll") {
      dateTemp <- rollDate(date, bdc = "Following", calendar, holidays)
      if (format(dateTemp, format = "%m")
          != format(date, format = "%m")){
        dateTemp <- rollDate(date, bdc = "Preceding", calendar, holidays)
      }
    }
  }
  
  return(dateTemp)
}

# -----------------------------------------------------------------------------
genSchedule <- function(settleDate, freq, tenor, calendar, bdc, holidays){
  # Calculates an array of dates, spaced by freq months until tenor years.
  # All dates are business days according to holiday calendar calendar and
  # adjusted according to business day count convention bdc.
  # Input
  # -- settleDate: string "YYYY-MM-DD" or a Date data type,
  #                the settlement date
  # -- freq: integer, frequency of payment in months, e.g., 3 for quarterly
  # -- tenor: integer, number of years until maturity
  # -- bdc: string, business day convention
  #    -- Actual: No rolling on the date applied even if it is a non-business day
  #    -- Preceding: 1st business day before holiday
  #    -- Following: 1st business day after holiday
  #    -- Modified_Prec: Same as "Preceding" unless it belongs to a different
  #                      month, in which case 1st business day after holiday
  #    -- Modified_Foll: Same as "Following" unless it belongs to a different
  #                      month, in which case 1st business day before holiday
  # -- calendar: string, indicate the desired calendar
  #    -- NY: New York holiday calendar
  #    -- General: all weekdays are business days
  
  settleDate <- as.Date(settleDate)
  checkCalendar(calendar)
  checkBDC(bdc)
  
  count <- tenor * 12 / freq
  month <- paste(toString(freq), "months", sep = " ")
  schedule <- seq.Date(settleDate, length.out = count + 1, by = month)
  
  for (i in 1:length(schedule)){
    schedule[i] <- rollDate(schedule[i], bdc = bdc, calendar = calendar, holidays = holidays)
  }
  return (schedule)
}

# ------------------------------------------------------------------------------
Initialize <- function(rate, tenor, fixFreq, fixDCC, fltFreq, fltDCC,
                       calendar, bdc, curveDate, numSetDay, yieldCurveDCC, holidays){
  # Rolls the settlement date based on yield curve observation date
  # curveDateAnd and number of settlement dates numSetDay
  # the first entry of discount vector to 1.
  # Input
  # -- rate: vector of doubles of zero coupon rates
  # -- tenor: vector of integers corresponding tenors
  # -- fixFreq: integer, fixed leg frequency of payment in months
  # -- fixDCC: string, fixed leg DCC
  # -- fltFreq: integer, floating leg frequency of payment in months
  # -- fltDCC: string, floating leg DCC
  # -- bdc: string, business day convention
  #    -- Actual: No rolling on the date applied even if it is a non-business day
  #    -- Preceding: 1st business day before holiday
  #    -- Following: 1st business day after holiday
  #    -- Modified_Prec: Same as "Preceding" unless it belongs to a different
  #                      month, in which case 1st business day after holiday
  #    -- Modified_Foll: Same as "Following" unless it belongs to a different
  #                      month, in which case 1st business day before holiday
  # -- calendar: string, indicate the desired calendar
  #    -- NY: New York holiday calendar
  #    -- General: all weekdays are business days
  # -- curveDate: string "YYYY-MM-DD" or a Date data type, yield curve date
  # -- numSetDay: integer, number of settlement days from yield curve date
  # -- yieldCurveDCC: Yield curve DCC
  
  checkBDC(bdc)
  checkCalendar(calendar)
  checkDCC(fixDCC)
  checkDCC(fltDCC)
  checkDCC(yieldCurveDCC)
  
  curveDate <- as.Date(curveDate)
  
  settleDate <- curveDate
  for (i in 1:numSetDay) {
    settleDate <- rollDate(settleDate + 1, bdc = "Following",
                           calendar = calendar, holidays)
  }
  # initialize curveDate to curveDate and df to 1
  curveDate <- c(curveDate)
  df <- c(1)
  return (data.frame(curveDate, df, settleDate))
}

# ------------------------------------------------------------------------------
logLinear <- function(dateIn, obsDate, discountFac, curveDate, yieldCurveDCC){
  # Calculates the log-linearly interpolated/extrapolated discount factor at
  # dateIn given the yield curve specified by the vector of discount factors
  # discountFac and their corresponding tenors. The yield curve date is
  # curveDate and the yield curve day count convention is yieldCurveDCC.
  # Input
  # -- dateIn: string "YYYY-MM-DD" or a Date data type, date of interest
  # -- obsDate: vector of curveDate, tenors on observed yield curve
  # -- discountFac: vector of doubles, discount factors on observe yield curve
  # -- curveDate: string "YYYY-MM-DD" or a Date data type, yield curve date
  # -- yieldCurveDCC: Yield curve DCC
  
  checkDCC(yieldCurveDCC)
  dateIn <- as.Date(dateIn)
  obsDate <- as.Date(obsDate)
  curveDate <- as.Date(curveDate)
  
  N <- length(obsDate)
  if (N < 2){
    stop("Need at least two observed values to interpolate")
  }
  
  if (dateIn < min(obsDate)){
    stop("Date is ealier than the curve date")
  }
  
  # if dateIn is already in obsDate, return its equivalent discount rate
  if (dateIn %in% obsDate){
    return (discountFac[which(obsDate == dateIn)])
  }
  
  # coding convenience: obsDate in the last interval unless identified
  # otherwise
  date1 <- obsDate[N - 1]
  date2 <- obsDate[N]
  dDF1 <- discountFac[N - 1]
  dDF2 <- discountFac[N]
  
  for (i in 1:(N - 1)) {
    if (dateIn > obsDate[i] && dateIn < obsDate[i + 1]){
      date1 <- obsDate[i]
      date2 <- obsDate[i + 1]
      dDF1 <- discountFac[i]
      dDF2 <- discountFac[i + 1]
      break
    }
  }
  
  date1Frac <- fracYear(curveDate, date1, yieldCurveDCC)
  date2Frac <- fracYear(curveDate, date2, yieldCurveDCC)
  dInFrac <- fracYear(curveDate, dateIn, yieldCurveDCC)
  
  # log-scale linear interpolation of discount factors
  temp <- log(dDF2) * (dInFrac - date1Frac) + log(dDF1) *
    (date2Frac - dInFrac)
  return (exp(temp / (date2Frac - date1Frac)))
}
# ------------------------------------------------------------------------------
pvSwap <- function(rate, tenor, fixFreq, fixDCC, fltFreq, fltDCC,
                   calendar, bdc, curveDate, numSetDay, yieldCurveDCC,
                   paymentDate, discountFac, settleDate, holidays){
  # Calculates the present value of an Interest rate Swap (IRS) by subtracting
  # the present value of floating leg from that of fix leg under specified
  # day count convention and holiday calendar.
  # Input
  # -- rate: double, given swap rate (fixed leg)
  # -- tenor: integer, given tenor of swap
  # -- fixFreq: integer, fixed leg frequency of payment in months
  # -- fixDCC: string, fixed leg DCC
  # -- fltFreq: integer, floating leg frequency of payment in months
  # -- fltDCC: string, floating leg DCC
  # -- bdc: string, business day convention
  #    -- Actual: No rolling on the date applied even if it is a non-business day
  #    -- Preceding: 1st business day before holiday
  #    -- Following: 1st business day after holiday
  #    -- Modified_Prec: Same as "Preceding" unless it belongs to a different
  #                      month, in which case 1st business day after holiday
  #    -- Modified_Foll: Same as "Following" unless it belongs to a different
  #                      month, in which case 1st business day before holiday
  # -- calendar: string, indicate the desired calendar
  #    -- NY: New York holiday calendar
  #    -- General: all weekdays are business days
  # -- curveDate: string "YYYY-MM-DD" or a Date data type, yield curve date
  # -- numSetDay: integer, number of settlement days from yield curve date
  # -- yieldCurveDCC: Yield curve DCC
  # -- paymentDate: payment dates
  # -- discountFac: discount factors on payment curveDate
  # -- settleDate: settlement date
  
  checkBDC(bdc)
  checkCalendar(calendar)
  checkDCC(fixDCC)
  checkDCC(fltDCC)
  checkDCC(yieldCurveDCC)
  curveDate <- as.Date(curveDate)
  settleDate <- as.Date(settleDate)
  
  # Calculates the present value of fixed leg
  fixSchedule <- genSchedule(settleDate, fixFreq, tenor, calendar, bdc, holidays)
  fixPV <- 0
  numFix <- length(fixSchedule)
  if (numFix < 2) {
    stop("Fixed leg has only one payment date.")
  }
  for (i in 2:numFix){
    df <- logLinear(fixSchedule[i], paymentDate, discountFac,
                    curveDate, yieldCurveDCC)
    fixPV <- fixPV + fracYear(fixSchedule[i - 1],
                              fixSchedule[i], fixDCC) * rate * df
  }
  
  # Calculates the present value of floating leg
  fltSchedule <- genSchedule(settleDate, fltFreq, tenor, calendar, bdc, holidays)
  fltPV <- 0
  numFlt <- length(fltSchedule)
  if (numFlt < 2) {
    stop("Floating leg has only one payment date.")
  }
  for (i in 2:numFlt){
    df <- logLinear(fltSchedule[i], paymentDate, discountFac,
                    curveDate, yieldCurveDCC)
    dt <- fracYear(fltSchedule[i - 1], fltSchedule[i], fltDCC)
    dfr <- (logLinear(fltSchedule[i - 1], paymentDate, discountFac,
                      curveDate, yieldCurveDCC) /
              logLinear(fltSchedule[i], paymentDate, discountFac,
                        curveDate, yieldCurveDCC) - 1) / dt
    fltPV <- fltPV + dt * dfr * df
  }
  return (fixPV - fltPV)
}
# ------------------------------------------------------------------------------
solve_rate <- function(index, swapRates, tenors, fixFreq, fixDCC,
                       fltFreq, fltDCC, calendar, bdc, curveDate,
                       numSetDay, yieldCurveDCC, paymentDate,
                       discountFac, settleDate, holidays){
  # Secant method to solve for the discount rate such that the present value
  # of a IRS equals zero
  # Input
  # -- index: integer, index of the input swap rate
  # -- swapRates: vector of doubles of swap rates
  # -- tenors: vector of integers of corresponding tenors
  # -- fixFreq: integer, fixed leg frequency of payment in months
  # -- fixDCC: string, fixed leg DCC
  # -- fltFreq: integer, floating leg frequency of payment in months
  # -- fltDCC: string, floating leg DCC
  # -- bdc: string, business day convention
  # -- calendar: string, indicate the desired calendar
  # -- curveDate: string "YYYY-MM-DD" or a Date data type, yield curve date
  # -- numSetDay: integer, number of settlement days from yield curve date
  # -- yieldCurveDCC: Yield curve DCC
  # -- paymentDate: payment dates
  # -- discountFac: discount factors on payment dates
  # -- settleDate: settlement date
  # -- holidays: user-defined holidays
  
  checkBDC(bdc)
  checkCalendar(calendar)
  checkDCC(fixDCC)
  checkDCC(fltDCC)
  checkDCC(yieldCurveDCC)
  curveDate <- as.Date(curveDate)
  settleDate <- as.Date(settleDate)
  
  # pick the first initial discount rate rate0
  rate0 <- exp(-swapRates[index] *
                 fracYear(as.Date(curveDate),
                          as.Date(paymentDate[index + 1]),
                          yieldCurveDCC))
  discountFac[index + 1] <- rate0
  df0 <- pvSwap(swapRates[index], tenors[index],
                fixFreq, fixDCC, fltFreq, fltDCC,
                calendar, bdc, curveDate, numSetDay, yieldCurveDCC,
                paymentDate, discountFac, settleDate, holidays)
  # pick the second initial discount rate rate1
  rate1 <- rate0 + 0.001
  discountFac[index + 1] <- rate1
  df1 <- pvSwap(swapRates[index], tenors[index],
                fixFreq, fixDCC, fltFreq, fltDCC,
                calendar, bdc, curveDate, numSetDay, yieldCurveDCC,
                paymentDate, discountFac, settleDate, holidays)
  # apply Secant method with threshold 1e-10
  while (abs(rate1 - rate0) >= 1e-10){
    dx <- rate0 - df0 * (rate1 - rate0) / (df1 - df0)
    discountFac[index + 1] <- dx
    df <- pvSwap(swapRates[index], tenors[index],
                 fixFreq, fixDCC, fltFreq, fltDCC,
                 calendar, bdc, curveDate, numSetDay, yieldCurveDCC,
                 paymentDate, discountFac, settleDate, holidays)
    rate0 <- rate1
    df0 <- df1
    rate1 <- dx
    df1 <- df
  }
  return (discountFac[index + 1])
}
# ------------------------------------------------------------------------------

#' Build Curve
#' 
#' @description  Bootstrap discount factors from a yield curve.
#'
#' @param swapRates A vector of doubles of swap rates.
#' @param tenors A vector of integers of corresponding tenors.
#' @param fixFreq An integer of fixed leg frequency of payment in months. 
#' 
#'    Default is 6, semi-annual payments.
#' @param fixDCC A string of fixed leg day count convention from four options:
#'   "Thirty360", "ACT360", "ACT365", or "ACTACT". 
#'   
#'   Default is "Thirty360".
#' @param fltFreq An integer of floating leg frequency of payment in months. 
#'    
#'    Default is 6, semi-annual payments.
#' @param fltDCC A string of floating leg day count convention from four options:
#'   "Thirty360", "ACT360", "ACT365", or "ACTACT". 
#'   
#'   Default is "Thirty360".
#' @param bdc A string of business day convention from five options: 
#' \itemize{
#'   \item "Actual": No rolling on the date applied even if it is a non-business day
#'   \item "Preceding": 1st business day before holiday
#'   \item "Following": 1st business day after holiday
#'   \item "Modified_Prec": Same as "Preceding" unless it belongs to a different
#'                      month, in which case 1st business day after holiday
#'   \item "Modified_Foll": Same as "Following" unless it belongs to a different
#'                      month, in which case 1st business day before holiday
#'  }
#'  Default is "Actual".
#' @param calendar A string of the desired calendar convention 
#'  from two options:
#'  \itemize{
#'   \item "NY": New York holiday calendar
#'   \item "General": all weekdays are business days
#'  }
#' @param curveDate A string in the format of "YYYY-MM-DD" of yield curve date.
#' @param numSetDay An integer of settlement days from yield curve date.
#' @param yieldCurveDCC A string of yield curve day count convention from four options:
#'   "Thirty360", "ACT360", "ACT365", or "ACTACT". Default is "Thirty360".
#' @param holidays An optional vector dates of user-defined holidays. If provided, 
#'    within the given holidays range, the calendar provided in the parameter "calendar" 
#'    will not be applied; 
#'    
#'    If the date is not in the given holidays range, it will follow the calendar provided in
#'    the "calendar" parameter
#' 
#' 
#' @return Outputs a data frame of strings of discount dates and doubles of
#'   discount factors.
#' @examples
#' rate <- c(0.69, 0.77, 0.88, 1.01, 1.14, 1.38, 1.66, 2.15) * 0.01
#' tenor <- c(1, 2, 3, 4, 5, 7, 10, 30)
#' fixFreq <- 6
#' fixDCC <- "Thirty360"
#' fltFreq <- 6
#' fltDCC <- "ACT360"
#' calendar <- "NY"
#' bdc <- "Modified_Foll"
#' curveDate <- "2016-02-08"
#' numSetDay <- 2
#' yieldCurveDCC <- "Thirty360"
#' holidays <- NULL
#' buildCurve(rate, tenor, fixFreq, fixDCC, fltFreq, fltDCC, calendar, bdc,
#'            curveDate, numSetDay, yieldCurveDCC, holidays)
#' @export
buildCurve <- function(swapRates, 
                       tenors, 
                       fixFreq = 6, 
                       fixDCC = "Thirty360",
                       fltFreq = 6, 
                       fltDCC = "Thirty360", 
                       calendar = "General", 
                       bdc = c("Actual", "Preceding", "Following", "Modified_Prec", "Modified_Foll"),
                       curveDate, 
                       numSetDay, 
                       yieldCurveDCC = "Thirty360",
                       holidays = NULL){
  
  checkBDC(bdc)
  checkCalendar(calendar)
  checkDCC(fixDCC)
  checkDCC(fltDCC)
  checkDCC(yieldCurveDCC)
  curveDate <- as.Date(curveDate)
  checkHolidays(holidays)
  
  
  if (length(swapRates) != length(tenors)){
    stop("Number of discount factors is different from the number of
           observation dates")
  }
  
  # temporary output: yield curve date dates, discount factor discountFac,
  # and settlement date settleDate ---------------------------------------
  outTemp <- Initialize(swapRates, tenors, fixFreq, fixDCC,
                        fltFreq, fltDCC, calendar, bdc,
                        curveDate, numSetDay, yieldCurveDCC, holidays)
  
  settleDate <- outTemp$settleDate
  
  numTenor <- length(tenors)
  obsDate <- rep(outTemp$curveDate, numTenor + 1)
  discountFac <- rep(outTemp$df, numTenor + 1)
  
  # Set up maturity dates in the form of "YYYY-MM-DD"
  for (i in 1:numTenor){
    posixDate <- as.POSIXlt(settleDate)
    posixDate$year <- posixDate$year + tenors[i]
    obsDate[i + 1] <- as.Date(posixDate)
  }
  
  # Discount factors
  for (i in 1:numTenor){
    discountFac[i + 1] <- solve_rate(i, swapRates, tenors, fixFreq,
                                     fixDCC, fltFreq, fltDCC, calendar,
                                     bdc, curveDate, numSetDay, yieldCurveDCC,
                                     obsDate[1:(i + 1)], discountFac[1:(i + 1)],
                                     settleDate, holidays)
  }
  
  # Zero rate
  zeroRate <- mat.or.vec(nr = numTenor + 1, nc = 1)
  dayCount <- mat.or.vec(nr = numTenor + 1, nc = 1)
  for (i in 1:(numTenor + 1)){
    dayCount[i] <- fracYear(curveDate, obsDate[i], yieldCurveDCC)
    if (dayCount[i] == 0) {
      zeroRate[i] <- 0
    } else {
      zeroRate[i] <- -log(discountFac[i]) / dayCount[i]
    }
  }
  
  # Forward Curve
  forwardCurve <- mat.or.vec(nr = numTenor + 1, nc = 1)
  for (i in 1:numTenor){
    forwardCurve[i] <- (logLinear(obsDate[i], obsDate, discountFac,
                                  curveDate, yieldCurveDCC) /
                          logLinear(obsDate[i + 1], obsDate, discountFac,
                                    curveDate, yieldCurveDCC) - 1) /
      fracYear(obsDate[i], obsDate[i + 1], yieldCurveDCC)
  }
  forwardCurve[numTenor + 1] <- forwardCurve[numTenor]
  
  return (data.frame(obsDate, discountFac, zeroRate, forwardCurve, dayCount))
}
