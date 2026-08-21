#' Calculate the expected (population) survival
#'
#' For a given individual with sex, year, and age, calculate the expected (population) survival at the supplied time points based on the mortality tables.
#'
#' The follow-up time and age must be specified in days. The calendar year can be in any date format
#' (Date and POSIXt are allowed)
#'
#'
#' @param sex Either character ('male'/'female'), or integer (1/2).
#' @param year The year from which the individual is followed. Either a Date or POSIXt object. Default is as.Date('1970-01-01').
#' @param age The age from which the individual is followed. Must be in days.
#' @param ratetable a table of event rates, organized as a \code{ratetable}
#' object, such as \code{slopop}.
#' @param times The times at which the expected (population) survival should be calculated. Must be in days.
#' @return A vector containing the survival estimate at the supplied times.
#' @seealso \code{\link{expprep2}}
#' @examples
#'
#' library(relsurv)
#' # Estimate P(T>2000 days) for a newborn:
#' popsurv(sex='male', year=as.Date('1970-01-01'), age=0, ratetable=slopop, times=2000)
#' # P(T>300 days) for a 50-year old:
#' popsurv(sex='male', year=as.Date('1970-01-01'), age=50*365.241, ratetable=slopop, times=300)
popsurv <- function(sex, year=as.Date('1970-01-01'), age=0, ratetable, times){

  # year <- as.Date(paste0(year, '-01-01'), origin=as.Date('1970-01-01'))

  data <- data.frame(sex=sex, year=year, age=age)

  data$time <- max(times)+age
  data$status <- 1

  if(age!=0) data$age <- 0

  times <- as.numeric(times)

  formula <- Surv(age, time, status)~1

  rform <- suppressWarnings(rformulate(formula, data, ratetable))

  data <- rform$data

  times_arg <- times
  if(age!=0){
    times_arg <- times_arg+age
    times_arg <- unique(c(age, times_arg))
    rform$R[,"year"] <- rform$R[,"year"]-age
  }

  temp <- exp_prep(rform$R[,,drop=FALSE], rform$Y,rform$ratetable,rform$status,times=times_arg,
                             fast=FALSE, cmp=FALSE, ys=as.numeric(age), netweiDM = TRUE)
  # if(age!=0 & length(times)>1) temp$yidli[1] <- 0

  # return(exp(-sum(temp$yidli)))
  # out <- exp(-tail(temp$yidli,1))

  out <- temp$yidli

  if(age!=0){
    out <- out[2:length(out)]
  }

  out <- exp(-cumsum(out))
  names(out) <- paste0('times=', times)

  return(out)
}

# TESTING:
# library(dplyr)
# library(survival)
# library(relsurv)
# library(ggplot2)
#
# # Define objects:
# slopop <- relsurv::slopop
# # slopop[,,] <- slopop["1","1970","male"]
# tajms <- 5*365.241
# xvred <- seq(0, 95, 1)
#
# # Calculate survival at different times since age 0:
# yvred <- sapply(xvred, function(xx) popsurv(sex='male', year=as.Date('1970-01-01'), 0, slopop, times=max(1,xx*365.241)))
# yvred[1] <- 1 # popravi prvo vrednost
# # Calculate 5-year conditional survival:
# rez <- data.frame(x=xvred[1:91], y=(yvred/lag(yvred, 5))[6:96])
#
# # Calculate 5-year conditional survival directly in popsurv:
# yvred2 <- sapply(xvred, function(xx) popsurv(sex='male', year=as.Date('1970-01-01')+xx*365.241,
#                                              slopop, times=tajms, age = xx*365.241))
#
# # Check:
# plot(rez$x, rez$y, xlab='Age', ylab='5-year survival', ylim=c(0,1), type='s')
# lines(rez$x[1:length(yvred2)], yvred2, col='red', type='s')
#
# plot(rez$x, rez$y - yvred2[1:nrow(rez)], type='s')
# # plot(rez$x, -log(rez$y) + log(yvred2[1:nrow(rez)]), type='s')
#
# ggplot(rez)+
#   geom_point(aes(x,y))+
#   theme_bw()+
#   scale_x_continuous(breaks = seq(0, 100, 10))+
#   # scale_y_continuous(breaks = seq(0, 1, 0.2), limits=c(0,1))+
#   # xlab('Age')+
#   # ylab('5-year probability of dying')+
#   # ggtitle('5-year probability of dying conditional on being alive at a given age')+
#   theme(plot.title = element_text(hjust = 0.5))
