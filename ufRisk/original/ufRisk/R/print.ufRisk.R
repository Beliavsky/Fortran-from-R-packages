#' Print Method for Objects of Class 'ufRisk'
#'
#'This function regulates how objects of class \code{ufRisk} are printed.
#'
#' @param x an object of class \code{ufRisk}; for the current package version,
#' only the function \code{br_test} returns such an object.
#' @param ... implemented for compatibility with the generic function;
#' additional arguments, however, will not affect this print method.
#'
#' @export
#'
#' @return
#' None
#'
#' @author
#'\itemize{
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University) \cr
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'}
#'

print.ufRisk <- function(x, ...) {
    if (attr(x, "function") == "trafftest") {
        cat(" ", fill = TRUE)
        cat("###################################", fill = TRUE)
        cat("#       Backtesting results       #", fill = TRUE)
        cat("###################################", fill = TRUE)
        cat(" ", fill = TRUE)
        cat("# Traffic light zone boundaries #", fill = TRUE)
        df <- data.frame(Zone = c("Green zone:", "Yellow zone:", "Red zone:"),
            Probability = c("p < 95%", "95% <= p < 99.99%", "p >= 99.99%"))
        print.data.frame(df, row.names = FALSE, quote = FALSE, right = FALSE)
        cat(" ", fill = TRUE)
        pot.vals = c(x[["pot_VaR.v"]], x[["pot_VaR.e"]], round(x[["br.sum"]], 4))
        p.vals <- c(x[["p_VaR.v"]], x[["p_VaR.e"]], x[["p_ES"]])
        WAD <- round(x[["WAD"]], 4)
        result <- rep(NA, times = 3)
        result[p.vals < 0.95] <- "Green zone"
        result[p.vals >= 0.95 & p.vals < 0.9999] <- "Yellow zone"
        result[p.vals >= 0.9999] <- "Red zone"
        p.vals <- round(p.vals, 4)
        cat(paste0("# Test 1: ", (1 - x[["a.v"]]) * 100, "%-VaR #"), fill = TRUE)
        cat(" Number of violations:", pot.vals[1], fill = TRUE)
        cat(" p = ", p.vals[[1]], ": ", result[[1]], fill = TRUE, sep = "")
        cat(" ", fill = TRUE)
        cat(paste0("# Test 2: ", (1 - x[["a.e"]]) * 100, "%-VaR #"), fill = TRUE)
        cat(" Number of violations:", pot.vals[2], fill = TRUE)
        cat(" p = ", p.vals[[2]], ": ", result[[2]], fill = TRUE, sep = "")
        cat(" ", fill = TRUE)
        cat(paste0("# Test 3: ", (1 - x[["a.e"]]) * 100, "%-ES #"), fill = TRUE)
        cat(" Number of weighted violations:", pot.vals[3], fill = TRUE)
        cat(" p = ", p.vals[[3]], ": ", result[[3]], fill = TRUE, sep = "")
        cat(" ", fill = TRUE)
        cat("# Weighted Absolute Deviation #", fill = TRUE)
        cat(" WAD = ", WAD, fill = TRUE, sep = "")
    }
    if (attr(x, "function") == "covtest") {
        cat(" ", fill = TRUE)
        cat("##################################", fill = TRUE)
        cat("#          Test results          #", fill = TRUE)
        cat("##################################", fill = TRUE)
        cat(" ", fill = TRUE)
        cat("# Unconditional coverage test #", fill = TRUE)
        cat(" ", fill = TRUE)
        cat("H0: w = ", x$p, sep = "", fill = TRUE)
        cat("p_[uc] = ", round(x$p.uc, 4), sep = "", fill = TRUE)
        if (x$p.uc < (1 - x$conflvl)){
          cat("Decision: Reject H0", fill = TRUE)
        } else {
          cat("Decision: Fail to reject H0", fill = TRUE)
        }
        cat(" ", fill = TRUE)
        cat("# Independence test #", fill = TRUE)
        cat(" ", fill = TRUE)
        cat("H0: w_[00] = w_[10]", sep = "", fill = TRUE)
        cat("p_[ind] = ", round(x$p.ind, 4), sep = "", fill = TRUE)
        if (x$p.ind < (1 - x$conflvl)){
          cat("Decision: Reject H0", fill = TRUE)
        } else {
          cat("Decision: Fail to reject H0", fill = TRUE)
        }
        cat(" ", fill = TRUE)
        cat("# Conditional coverage test #", fill = TRUE)
        cat(" ", fill = TRUE)
        cat("H0: w_[00] = w_[10] = ", x$p, sep = "", fill = TRUE)
        cat("p_[cc] = ", round(x$p.cc, 4), sep = "", fill = TRUE)
        if (x$p.cc < (1 - x$conflvl)){
          cat("Decision: Reject H0", fill = TRUE)
        } else {
          cat("Decision: Fail to reject H0", fill = TRUE)
        }
        cat(" ", fill = TRUE)
    }
}
