#' @title Calculate Credit Default Swap rates
#'
#' @description
#' Calculate CDS rates starting from default intensities
#'
#' @param t premium timetable.
#' @param Tj CDS maturities.
#' @param tr interest rates timetable.
#' @param int default intensities timetable.
#' @param tint intensity timetable.
#' @param r spot interest rates.
#' @param R constant premium payment.
#' @param ... further arguments on \code{cds}.
#'
#' @return
#' An object of class \code{data.frame} that contains the quantities calculated by \code{cds}
#' on Tj timetable.
#'
#' @details
#' The function \code{cds2} is based on \code{cds} but allows a more fine controll on maturities
#' and on discretization of \code{r} and \code{int}. In particular input \code{(t, tr, tint)}
#' can be of different length thanks to the function \link[stats]{approx}.
#'
#' @references
#' David Lando (2004) Credit Risk Modeling.
#'
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes
#'
#' @examples
#' cds2(t = c(1:20),Tj = c(1:20), tr = c(1:20), r = seq(0.01,0.06, len =20),
#' tint = c(1:20), int= seq(0.01,0.06, len =20))
#'
#' @export
#'
cds2 = function (t, Tj, tr, r, tint, int, R = 0.005, ...)
{
  if (!is.numeric(t))
    stop("time t must be a numeric vector")
  if (!is.numeric(int))
    stop("intensity int must be a numeric vector")
  if (!is.numeric(r))
    stop("rates r must be a numeric vector")

  new.r <- stats::approx(x = tr, y = r, xout = t, method = "linear",
                         rule = 2)$y
  new.int <- stats::approx(x = tint, y = int, xout = t, method = "linear",
                           rule = 2)$y

  cds1 <- cds(t = t, r = new.r, int = new.int, R = R, ...)
  ##Output
  Q <- stats::approx(x = t, y = cds1$Survival, xout = Tj, rule = 2)$y
  price <- stats::approx(x = t, y = cds1$Price, xout = Tj, rule = 2)$y
  premium <- stats::approx(x = t, y = cds1$PremiumLeg, xout = Tj,
                           rule = 2)$y
  protection <- stats::approx(x = t, y = cds1$ProtectionLeg,
                              xout = Tj, rule = 2)$y
  out <- data.frame(Tj = Tj, Q = Q, premium = premium, protection = protection,
                    rate = R * protection/premium, price = price)
  names(out) <- c("Tj", "Survival", "PremiumLeg", "ProtectionLeg",
                  "Rate", "Price")
  return(out)
}
