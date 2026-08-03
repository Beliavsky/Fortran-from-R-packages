#' Perform Portfolio Analysis
#'
#' This function performs portfolio analysis using various methods such as Mean-Variance (MV),
#' James-Stein (JM), LASSO, Ridge Regression, Equal Weighting (EW), among others.
#' It calculates weights, turnover, returns, Sharpe ratios, volatility, and maximum drawdown for each method.
#'
#' @param x A numeric matrix where each column represents asset returns and rows represent time periods.
#' @param mon A numeric vector representing the number of months since the start date for each time period.
#' @param count A numeric vector indicating the number of entries per month.
#' @param Date A vector of Date objects representing the dates of the time periods.
#' @param num_cores The number of cores to use for parallel processing. Default is 7.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{cumulative_return_plot}{A ggplot object representing the cumulative returns for each method.}
#'   \item{cumulative_turnover_plot}{A ggplot object representing the cumulative turnover for each method.}
#'   \item{turnover_mean}{A numeric vector of the mean turnover for each method.}
#'   \item{sharpe_ratio}{A numeric vector of the Sharpe ratio for each method.}
#'   \item{volatility}{A numeric vector of the annualized volatility for each method.}
#'   \item{max_drawdown}{A numeric vector of the maximum drawdown for each method.}
#'   \item{vw_to_mean}{The mean turnover for the volume-weighted (VW) portfolio.}
#'   \item{vw_sharpe_ratio}{The Sharpe ratio for the VW portfolio.}
#'   \item{vw_volatility}{The annualized volatility for the VW portfolio.}
#'   \item{vw_max_drawdown}{The maximum drawdown for the VW portfolio.}
#' }
#'
#' @details The function iterates through different time periods and calculates portfolio weights,
#' turnover, and returns for multiple methods including Mean-Variance (MV), James-Stein (JM), and various
#' regularization techniques. It also computes performance metrics like the Sharpe ratio, volatility,
#' maximum drawdown, and cumulative turnover for each method. Visualization of the cumulative returns
#' and turnover is generated using ggplot2.
#'
#' @examples
#' \donttest{
#' # Create a larger example dataset that aligns with the function's expectations
#' set.seed(123)
#' x <- matrix(runif(700), ncol = 10)  # 10 columns (assets), 70 rows (observations)
#' mon <- rep(1:10, each = 7)          # Example month identifiers, 7 observations per month
#' count <- rep(7, 10)                 # Example count per month (7 entries per month)
#' Date <- as.Date('2020-01-01') + 0:69 # Example date sequence (70 days)
#'
#' # Run the analysis with 2 cores
#' result <- perform_analysis(x, mon, count, Date, num_cores = 2)
#'
#' # Display results
#' print(result$cumulative_return_plot)
#' print(result$cumulative_turnover_plot)
#' }
#' @export
perform_analysis <- function(x, mon, count, Date, num_cores = 7) {
  p <- ncol(x)

  w.frame <- replicate(13, matrix(data = NA, nrow = max(mon) - 6, ncol = p))
  stat.frame <- matrix(data = NA, nrow = max(mon) - 6, ncol = 13)

  w <- w.frame
  turn_over <- stat.frame
  w0.tmp <- matrix(data = 1/p, nrow = p, ncol = 13)
  r <- matrix(NA, sum(count[7:length(count)]), 13)

  # Main analysis loop
  for (i in (7:max(mon))) {
    w.tmp <- matrix(data = NA, nrow = p, ncol = 13)
    x_tmp <- x[which((mon < i) & (mon >= i - 6)), ]
    n_tmp <- nrow(x_tmp)
    group <- buh.clust(x_tmp)

    # Perform calculations for different methods (ensure these functions are defined)
    w.tmp[,1] <- po.cols(rep(0, n_tmp), x_tmp)                                                 # MV
    w.tmp[,2] <- po.JM(x_tmp)                                                                 # JM
    w.tmp[,3] <- po.avg(rep(0, n_tmp), x_tmp, method = "LASSO")                               # TRP_min
    w.tmp[,4] <- po.avg(rep(0, n_tmp), x_tmp, method = "RIDGE")                               # APP-Ridge
    w.tmp[,5] <- po.grossExp(rep(0, n_tmp), x_tmp, method = "NOSHORT")                        # Fan.et.al_JM
    w.tmp[,6] <- po.grossExp(rep(0, n_tmp), x_tmp, method = "EQUAL")                          # FZY
    w.tmp[,7] <- rep(1/p, p)                                                                  # EW
    w.tmp[,9] <- po.bhu(rep(0, n_tmp), x_tmp, group, 100)                                     # TRP_clu
    w.tmp[,10] <- po.TZT(x_tmp, gamma = 3)                                                    # TZ
    w.tmp[,11] <- po.SW(x_tmp, b = round(p^0.7), sample = 1000)                               # SW
    w.tmp[,12] <- po.SW(x_tmp, b = round(length(group)), sample = 1000)                       # SW_clu
    w.tmp[,13] <- po.SW.lasso(rep(0, n_tmp), x_tmp, b = round(length(group)), sample = 1000)  # SW_shrink

    # Store weights and turnover
    for (j in 1:13) {
      w[i - 6,,j] <- w.tmp[,j]
      turn_over[i - 6,j] <- norm(as.matrix(w.tmp[,j]) - w0.tmp[,j])
    }

    # Define r_tmp before using it
    r_tmp <- (1 + (x[which(mon == i), ] / 100))

    # Calculate returns
    for (j in 1:13) {
      for (k in 1:nrow(r_tmp)) {
        idx <- which(mon == i)[k] - sum(count[1:6])
        if (idx > 0 && idx <= nrow(r)) {
          r[idx,j] <- r_tmp[k,] %*% w0.tmp[,j]
          w0.tmp[,j] <- (r_tmp[k,] * w0.tmp[,j]) / sum(r_tmp[k,] * w0.tmp[,j])
        }
      }
      w0.tmp[,j] <- w0.tmp[,j] / sum(w0.tmp[,j])
    }

    w0.tmp <- w.tmp
    message(sprintf("Processed month %d", i - 6))
  }

  # Performance Metrics
  colnames(turn_over) <- c("MV","JM","TRP_min","APP-Ridge","Fan.et.al_JM","FZY","EW","LW","TRP_clu","TZ","SW","SW_clu","SW_shrink")
  turnover_mean <- colMeans(turn_over, na.rm = TRUE) * 100

  sharpe_ratio <- sqrt(252) * (colMeans(r - 1, na.rm = TRUE) / apply(r - 1, 2, sd, na.rm = TRUE))
  volatility <- sqrt(252) * apply(r - 1, 2, sd, na.rm = TRUE) * 100

  v <- apply((1 + (r - 1)), 2, cumprod)
  v1 <- apply(v, 2, max, na.rm = TRUE)
  v1 <- matrix(rep(v1, each = nrow(v)), ncol = ncol(v))
  mdd <- (v1 - v) / v1 * 100
  max_drawdown <- apply(mdd, 2, max, na.rm = TRUE)

  # VW Calculation
  r_vw <- matrix(NA, sum(count[7:length(count)]), 1)
  w_vw <- matrix(data = NA, nrow = max(mon) - 6, ncol = p)

  vw0.tmp <- matrix(data = 1/p, nrow = p, ncol = 1)
  vw.to <- matrix(data = NA, nrow = max(mon) - 6, ncol = 1)

  for (i in 7:max(mon)) {
    x_tmp <- x[which((mon < i) & (mon >= i - 6)), ]
    n_tmp <- dim(x_tmp)[1]

    r_tmp <- (1 + (x[which(mon == i), ] / 100))

    vw.tmp <- vw0.tmp

    for (k in 1:nrow(r_tmp)) {
      idx <- which(mon == i)[k] - sum(count[1:6])
      if (idx > 0 && idx <= nrow(r_vw)) {
        r_vw[idx] <- r_tmp[k,] %*% vw0.tmp
        vw.tmp <- (r_tmp[k,] * vw0.tmp) / sum(r_tmp[k,] * vw0.tmp)
      }
    }
    w_vw[i - 6, ] <- vw.tmp
    vw.to[i - 6, ] <- norm(as.matrix(vw.tmp) - vw0.tmp)
    vw0.tmp <- vw.tmp
  }

  vw_to_mean <- colMeans(vw.to, na.rm = TRUE) * 100
  vw_sharpe_ratio <- sqrt(252) * (mean(r_vw - 1, na.rm = TRUE) / sd(r_vw - 1, na.rm = TRUE))
  vw_volatility <- sqrt(252) * sd(r_vw - 1, na.rm = TRUE) * 100

  vw_v <- cumprod(1 + (r_vw - 1))
  vw_v1 <- max(vw_v, na.rm = TRUE)
  vw_v1 <- rep(vw_v1, length(vw_v))  # Corrected line
  vw_mdd <- (vw_v1 - vw_v) / vw_v1 * 100
  vw_max_drawdown <- max(vw_mdd, na.rm = TRUE)

  ## Visualization
  # Combine all the returns for visualization
  r.final <- cbind(r, r_vw) # Combine all the returns of different methods
  r.final <- r.final[, colSums(is.na(r.final)) == 0] # Remove empty columns

  # Cumulative returns
  cum.r <- apply(r.final - 1, 2, cumsum)

  # Ensure cum.r and Date have matching lengths
  if (nrow(cum.r) != length(Date[-seq(1:sum(count[1:6]))])) {
    stop("Mismatch in rows between Date and cumulative returns (cum.r).")
  }

  df <- data.frame(time = Date[-seq(1:sum(count[1:6]))], cum.r)
  colnames(df) <- c("time", "MV", "JM", "FZY", "EW", "REN", "TZ", "VW")

  df <- melt(df, id.vars = "time")
  colnames(df) <- c("time", "Method", "value")

  p1.ff25 <- ggplot(df, aes(x = time, y = value, color = Method)) + theme_bw() + geom_line() +
    ggtitle("Cumulative Returns") +
    scale_color_manual(values = c('red','cyan','magenta','grey47','black','lightblue',"green","orange","brown")) +
    labs(y = "Cumulative Return", x = "Time") +
    theme(legend.justification = c(0, 1), legend.position = c(0, 1)) + theme(legend.position = "bottom")

  # Cumulative turnover
  to.final <- cbind(turn_over, vw.to) # Combine all the TOs of different methods
  to.final <- to.final[, colSums(is.na(to.final)) == 0] # Remove empty columns

  # Create time vector matching the turnover data
  time <- sapply(1:nrow(to.final), function(i) Date[which(mon == (i + 6))[1]])

  df2 <- data.frame(time = time, turnover = apply(abs(to.final), 2, cumsum))
  colnames(df2) <- c("time", "MV", "JM", "FZY", "EW", "REN", "TZ", "VW")

  df2 <- melt(df2, id.vars = "time")
  colnames(df2) <- c("time", "Method", "value")

  p2.ff25 <- ggplot(df2, aes(x = time, y = value, color = Method)) + theme_bw() + geom_line() +
    ggtitle("Cumulative Turnover") +
    scale_color_manual(values = c('red','cyan','magenta','grey47','black','lightblue',"green","orange","brown")) +
    labs(y = "Cumulative Turnover", x = "Time") +
    theme(legend.justification = c(0, 1), legend.position = c(0, 1)) + theme(legend.position = "bottom")

  return(list(
    cumulative_return_plot = p1.ff25,
    cumulative_turnover_plot = p2.ff25,
    turnover_mean = turnover_mean,
    sharpe_ratio = sharpe_ratio,
    volatility = volatility,
    max_drawdown = max_drawdown,
    vw_to_mean = vw_to_mean,
    vw_sharpe_ratio = vw_sharpe_ratio,
    vw_volatility = vw_volatility,
    vw_max_drawdown = vw_max_drawdown
  ))
}
