#' @import ggplot2
#' @import reshape2

fut_to_zero_eff <- function(fut) {
  (cumprod(1 + fut)) ^ (1 / (seq_along(fut))) - 1
}

disc_to_swap <- function(disc) {
  (1 - disc) / cumsum(disc)
}

swap_to_disc <- function(swap) {
  d = swap
  d[1] = 1 / (1 + swap[1])
  for (j in 2:length(swap)) {
    d[j] = max((1 - sum(d[1:(j - 1)])*swap[j]) / (1 + swap[j]), 0)
  }
  d
}

disc_to_zero_eff <- function(disc) {
  (1 / disc) ^ (1 / (seq_along(disc))) - 1
}

zero_eff_to_disc <- function(zero) {
  1 / ((1 + zero) ^ (seq_along(zero)) )
}

disc_to_zero_nom <- function(disc) {
  (1 / disc - 1) / seq_along(disc)
}

zero_nom_to_disc <- function(zero) {
  1 / (1 + zero * seq_along(zero))
}

disc_to_fut <- function(disc) {  
  exp(-diff(log(c(1,disc)))) - 1  
}

fut_to_disc <- function(fut) {
  1 / cumprod(1 + fut)
}

disc_to_german <- function(disc) {
  vapply(
    1:length(disc),
    function(i) (1 - 1 / i * sum(disc[1:i]) ) / sum( disc[1:i] * (1 - (1:i - 1 ) / i ) ),
    1)  
}

disc_to_french <- function(disc, search_interval = c(0.0001,1), tol = 1e-8) {  
  zerome = function(r,i,disc) 1/r * (1 - (1 + r) ^ (-i) ) - sum(disc[1:i])
  vapply(
    1:length(disc),
    function(i) uniroot(function(r) zerome(r,i,disc), interval = search_interval, tol = tol)$root,
    1)
}

disc_to_zero_cont <- function(disc) {
  - log(disc) / seq_along(disc)
}

zero_cont_to_disc <- function(zero) {
  exp(- zero * seq_along(zero))
}

eff_to_dir <- function(r) {
  (1 + r) ^ seq_along(r) - 1
}

dir_to_eff <- function(r) {
  (1 + r) ^ (1 / seq_along(r)) - 1
}

unscale_nom <- function(x, rate_scale) x / rate_scale

rescale_nom <- function(x, rate_scale) x * rate_scale

unscale_eff <- function(x, rate_scale) (1 + x) ^ (1 / rate_scale) - 1

rescale_eff <- function(x, rate_scale) (1 + x) ^ (rate_scale) - 1

unscale <- function(x, rate_scale, rate_type) {
  if (rate_type %in% c("zero_nom", "german", "french", "swap", "fut", "zero_cont")) {
    unscale_nom(x, rate_scale)
  } else {
    unscale_eff(x, rate_scale)
  }
}

rescale <- function(x, rate_scale, rate_type) {
  if (rate_type %in% c("zero_nom", "german", "french", "swap", "fut", "zero_cont")) {
    rescale_nom(x, rate_scale)
  } else {
    rescale_eff(x, rate_scale)
  }
}

#' @title Creates a rate curve instance
#' 
#' @param rates A rate vector
#' @param rate_type The rate type. Must be on of c("fut", "zero_nom", "zero_eff", "swap", "zero_cont)
#' @param pers The periods the rates correspond to
#' @param rate_scale In how many periods is the rate expressed.
#' For example, when measuring periods in days, and using annual rates, you should use 365. 
#' When measuring periods in months, and using annual rates, you should use 12.
#' If no scaling, use 1.
#' @param fun_d A discount factor function. fun_d(x) returns the discount factor for time x, vectorized on x
#' @param fun_r A rate function. fun_r(x) returns the EPR for time x, vectorized on x
#' @param knots The nodes used to bootstrap the rates. This is a mandatory argument if a rate function or discount function is provided
#' @param functor A function with parameters x and y, that returns a function used to interpolate
#' 
#' @note Currently a rate curve can only be built from one of the following sources
#' \enumerate{
#' \item A discount factor function
#' \item A rate function and a rate type from the following types: "fut", "zero_nom", "zero_eff", "swap" or "zero_cont
#' \item A rate vector, a pers vector and a rate type as before
#' }
#' @examples
#' rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
#' rate_curve(fun_r = function(x) rep_len(0.1, length(x)), rate_type = "swap", knots = 1:12)
#' rate_curve(fun_d = function(x) 1 / (1 + x), knots = 1:12)
#' @export
rate_curve <- function(
  rates = NULL,
  rate_type = "zero_eff",
  pers = 1:length(rates),
  rate_scale = 1,
  fun_d = NULL,
  fun_r = NULL,
  knots = seq.int(from = 1, to = max(pers), by = 1),
  functor = function(x, y) splinefun(x = x, y = y, method = "monoH.FC")) {
  if (!(
    (!is.null(fun_d) &&  !is.null(knots)) ||
      (!is.null(fun_r) && !is.null(rate_type) && !is.null(knots)) ||
      (!is.null(rates) && !is.null(rate_type))
    ))
    stop("A rate or discount function must be given to create a rate_curve")
  if (!is.null(fun_d)) {    
    r <- structure(list(), class = "rate_curve")
    r$f <- fun_d
    r$knots <- knots
    r$functor <- functor
    r$rate_scale <- rate_scale
    r
  } else if (!is.null(fun_r)) {
    x <- unscale(x = fun_r(knots), rate_scale = rate_scale, rate_type = rate_type)
    d <- do.call(what = paste0(rate_type,"_to_disc"), args = list(x))
    f <- functor(x = c(0 , knots), y = c(1 , d))    
    rate_curve(fun_d = f, knots = knots, functor = functor, rate_scale = rate_scale)
  } else if (!is.null(rates)) {
    f <- functor(x = pers, y = rates)
    rate_curve(fun_r = f, rate_type = rate_type, knots = knots, functor = functor, rate_scale = rate_scale)
  } else {
    stop("The rate_curve constructor lacks arguments")
  }  
}

get_rate_fun <- function(r, rate_type = "zero_eff") {
  stopifnot(rate_type %in% c("french", "fut", "german", "zero_eff", "zero_nom", "swap", "zero_cont"))
  d <- (r$f)(r$knots)
  y <- do.call(what = paste0("disc_to_",rate_type), args = list(d))
  f <- r$functor(x = r$knots, y = y)
  function(x) rescale(f(x), rate_scale = r$rate_scale, rate_type = rate_type)
}

#' @title Returns a particular rate or rates from a curve
#' 
#' @param r The rate_curve object
#' @param x The points in time to return
#' @param rate_type The rate type
#' 
#' @return If \code{x} is \code{NULL}, then returns a rate function of \code{rate_type} type.
#' Else, it returns the rates of \code{rate_type} type and corresponding to time \code{x}
#' @examples
#' r <- rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
#' r["zero_eff"]
#' r["swap",c(1.5, 2)]
#' @export
`[.rate_curve` <- function(r, rate_type = "zero_eff", x = NULL) {
  f <- get_rate_fun(r = r, rate_type = rate_type)
  if (is.null(x))
    f
  else
    f(x)
}

#' @title Plots a rate curve
#' 
#' @param x The rate curve
#' @param rate_type The rate types to plot, in c("french", "fut", "german", "zero_eff", "zero_nom", "swap", "zero_cont")
#' @param y_labs_perc If TRUE, the y axe is labeled with percentages 
#' @param y_labs_acc If y_labs_perc is TRUE, the accuracy for the percentages (i.e., 1 for xx\%, 0.1 for xx.x\%, 0.01 for xx.xx\%, etc) 
#' @param ... Other arguments (unused)
#' @examples
#' r <- rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
#' plot(r)
#' \dontrun{
#' plot(r, rate_type = "german")
#' plot(r, rate_type = c("french", "german"))
#' }
#' @export
plot.rate_curve <- function(x, rate_type = NULL, y_labs_perc = TRUE, y_labs_acc = NULL, ...) {
  all_rate_types = c("french", "fut", "german", "zero_eff", "zero_nom", "swap", "zero_cont")
  if (is.null(rate_type)) {
    rate_type = all_rate_types
  }
  if (any(!(rate_type %in% all_rate_types))) {
    stop("The rate type is not recognized")
  }
  df <- as.data.frame(lapply(
    X = rate_type,
    FUN = function(z)  x[z, x$knots]
    ))
  names(df) <- rate_type
  df$Time = x$knots    
  dfm <- reshape2::melt(data = df, id.vars = "Time", variable.name = "RateType", value.name = "Rate")
  x <- ggplot2::ggplot(data = dfm) +
    ggplot2::geom_line(mapping = ggplot2::aes_string(x = "Time", y = "Rate", color = "RateType"))
  if (y_labs_perc) x <- x + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = y_labs_acc))
  x
}

#' @title Calculates the present value of a cashflow
#' 
#' @param r A rate curve
#' @param cf The vector of values corresponding to the cashflow
#' @param d The periods on which the cashflow occurs. If missing, it is assumed that cf[i] occurs on period i
#' 
#' @return The present value of the cashflow
#' 
#' @examples
#' r <- rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
#' disc_value(r, cf = c(-1, 1.10), d = c(0,1))
#' disc_value(r, cf = c(-1, 1.15*1.15), d = c(0,2))
#' 
#' @export
disc_value <- function(r, cf, d = 1:length(cf)) {
  sum(r$f(d) * cf)
}