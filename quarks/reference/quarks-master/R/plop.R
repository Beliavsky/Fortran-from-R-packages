#' Profit & Loss operator function
#'
#' Calculates portfolio returns or losses by assigning weights
#'
#' @param x a numeric matrix of asset returns or losses
#' @param wts a numeric vector or matrix containing the portfolio
#' weights; portfolio value is standardized to 1 on any observation unit;
#' sum of weights should not exceed 1 (row-wise for matrices); by default
#' the portfolio is equally weighted over time and across all assets; if a
#' vector is passed to \code{wts} the portfolio is equally weighted over time
#' @param approxim controls if a first-order approximation for the calculation
#' of returns or losses is used; default is \code{1} (first-order
#' approximation is employed)
#'
#' @export
#'
#' @return Returns a list with the following elements:
#' \describe{
#' \item{pl}{Weighted portfolio returns or losses}
#' \item{wts}{Portfolio weights}
#' }
#' @examples
#' # creating portfolio
#' portfol <- cbind(SP500$price_close, DJI$price_close)
#' returns <- apply(portfol, 2, function(x) diff(log(x)))
#' # defining weights and applying the P&L operator function
#' wts <- c(0.4, 0.6)
#' portret <- plop(returns, wts = wts, approxim = 1)
#' portloss <- plop(-returns, wts = wts, approxim = 1)
#' plot.ts(cbind(portret$pl, portloss$pl))
#'
plop <-  function(x, wts = NULL, approxim = c(0, 1)) {
  if ((dim(wts)[1] == 1 || is.null(dim(wts)[1])) &&
      (!is.null(wts) && sum(wts) > 1)) {
    warning("Sum of weights in 'wts' exceeding 1.")
  }
  if ((dim(wts)[1] > 1 && !is.null(dim(wts)[1])) &&
      (!is.null(wts) && rowSums(wts) > 1)) {
    warning("Sum of weights in some row(s) of 'wts' exceeding 1.")
  }
  if (!length(approxim) %in% c(1, 2) || any(is.na(approxim)) ||
      !is.numeric(approxim) || !all(approxim %in% c(0, 1))) {
    stop("A single integer value must be passed to 'approxim'.",
         "Valid choices are 0 or 1.")
  }
  if (all(approxim == c(0, 1))) approxim = 1
  x <- as.matrix(x)
  ncol <- dim(x)[2]

  if (is.null(wts)) wts <- rep(1 / ncol, ncol)
  if (dim(wts)[1] == 1 && is.matrix(wts)) wts <- as.vector(wts)
  if (dim(wts)[1] == 1 || is.null(dim(wts)[1])) op <- '%*%'
  if (dim(wts)[1] > 1 && !is.null(dim(wts)[1])) op <- '*'

  switch(as.character(approxim),
         "0" =  {
           pl <- rowSums(get(op)(exp(x), wts) - 1)
         },
         "1" = {
           pl <- rowSums(get(op)(x, wts))
         })
  return(list(pl = pl,
              wts = wts))
}


