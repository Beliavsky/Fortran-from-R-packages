#' Print Method for the Package 'quarks'
#'
#'This function regulates how objects created by the package \code{quarks} are
#'printed.
#'
#'@param x an input object of class \code{quarks}.
#'@param ... included for compatibility; additional arguments will however
#'not affect the output.
#'
#'@export
#'
#'@return
#'None

print.quarks <- function(x, ...) {
  if(attr(x, "function") == "rollcast") {
    if (x[["method"]] == "plain") {
      method = "Plain"
    }
    if (x[["method"]] == "age") {
      method = "Age Weighting"
    }
    if (x[["method"]] == "vwhs") {
      method = "Volatility Weighting"
    }
    if (x[["method"]] == "fhs") {
      method = "Filtered"
    }
    VaR <- x[["VaR"]]
    ES <- x[["ES"]]
    loss <- -x[["xout"]]
    p <- x[["p"]]
    It.var <- ifelse(VaR <= loss, 1, 0)
    br.var <- sum(It.var)
    It.es <- ifelse(ES <= loss, 1, 0)
    br.es <- sum(It.es)
    cat(" ", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("|              Specifications              |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat(" Out-of-sample size:   ", x[["nout"]], fill = TRUE)
    cat(" Rolling window size:  ", x[["nwin"]], fill = TRUE)
    cat(" Bootstrap sample size:", x[["nboot"]], fill = TRUE)
    cat(" Confidence level:     ", paste(p * 100, "%"), fill = TRUE)
    cat(" Method:               ", method, fill = TRUE)
    cat(" Model:                ", x[["model"]], fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("|           Forecast performance           |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat(" Out-of-sample losses exceeding VaR", fill = TRUE)
    cat(" ", fill = TRUE)
    cat(" Number of breaches:   ", br.var, fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat(" Out-of-sample losses exceeding ES", fill = TRUE)
    cat(" ", fill = TRUE)
    cat(" Number of breaches:   ", br.es, fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
  }

  if (attr(x, "function") == "trftest") {
    if (x[["method"]] == "plain") {
      method = "Plain"
    }
    if (x[["method"]] == "age") {
      method = "Age Weighting"
    }
    if (x[["method"]] == "vwhs") {
      method = "Volatility Weighting"
    }
    if (x[["method"]] == "fhs") {
      method = "Filtered"
    }
    cat(" ", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("|            Traffic Light Test            |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat(" Method:", method, fill = TRUE)
    if (method %in% c("Volatility Weighting", "Filtered")) {
    cat(" Model: ", x[["model"]], fill = TRUE)
    }
    cat("--------------------------------------------", fill = TRUE)
    cat("|       Traffic light zone boundaries      |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    df <- data.frame(Zone = c("Green zone:", "Yellow zone:", "Red zone:"),
                     Probability = c("p < 95%", "95% <= p < 99.99%", "p >= 99.99%"))
    print.data.frame(df, row.names = FALSE, quote = FALSE, right = FALSE)
    pot.vals = x[["pot_VaR"]]
    p.vals <- x[["p_VaR"]]
    result <- rep(NA, times = 3)
    result[p.vals < 0.95] <- "Green zone"
    result[p.vals >= 0.95 & p.vals < 0.9999] <- "Yellow zone"
    result[p.vals >= 0.9999] <- "Red zone"
    p.vals <- round(p.vals, 4)
    cat("--------------------------------------------", fill = TRUE)
    if (length(strsplit(paste((1 - x[["p"]]) * 100), split = ".", fixed = TRUE)[[1]]) > 1) {
      cat(paste0("|         Test result - ", (1 - x[["p"]]) * 100, "%-VaR          |"), fill = TRUE)
    } else {
      cat(paste0("|         Test result - ", (1 - x[["p"]]) * 100, "%-VaR            |"), fill = TRUE)
    }
    cat("--------------------------------------------", fill = TRUE)
    cat(" Number of violations:", pot.vals[1], fill = TRUE)
    cat(" ", fill = TRUE)
    cat(" p = ", p.vals[[1]], ": ", result[[1]], fill = TRUE, sep = "")
    cat("--------------------------------------------", fill = TRUE)
  }
  if (attr(x, "function") == "cvgtest") {
    if (x[["method"]] == "plain") {
      method = "Plain"
    }
    if (x[["method"]] == "age") {
      method = "Age Weighting"
    }
    if (x[["method"]] == "vwhs") {
      method = "Volatility Weighting"
    }
    if (x[["method"]] == "fhs") {
      method = "Filtered"
    }
    if (x[["method"]] == "NA") {
      method = "NA"
    }
    cat(" ", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("|               Test results               |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat(" Method:", method, fill = TRUE)
    cat(" Model: ", x[["model"]], fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("|        Unconditional coverage test       |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("H0: w = ", x[["p"]], sep = "", fill = TRUE)
    cat(" ", fill = TRUE)
    cat("p_[uc] = ", round(x[["p.uc"]], 4), sep = "", fill = TRUE)
    cat(" ", fill = TRUE)
    if (x[["p.uc"]] < (1 - x[["conflvl"]])) {
      cat("Decision: Reject H0", fill = TRUE)
    } else {
      cat("Decision: Fail to reject H0", fill = TRUE)
    }
    cat("--------------------------------------------", fill = TRUE)
    cat("|            Independence test             |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("H0: w_[00] = w_[10]", sep = "", fill = TRUE)
    cat(" ", fill = TRUE)
    cat("p_[ind] = ", round(x[["p.ind"]], 4), sep = "", fill = TRUE)
    cat(" ", fill = TRUE)
    if (x[["p.ind"]] < (1 - x[["conflvl"]])) {
      cat("Decision: Reject H0", fill = TRUE)
    } else {
      cat("Decision: Fail to reject H0", fill = TRUE)
    }
    cat("--------------------------------------------", fill = TRUE)
    cat("|         Conditional coverage test        |", fill = TRUE)
    cat("--------------------------------------------", fill = TRUE)
    cat("H0: w_[00] = w_[10] = ", x[["p"]], sep = "", fill = TRUE)
    cat(" ", fill = TRUE)
    cat("p_[cc] = ", round(x[["p.cc"]], 4), sep = "", fill = TRUE)
    cat(" ", fill = TRUE)
    if (x[["p.cc"]] < (1 - x[["conflvl"]])) {
      cat("Decision: Reject H0", fill = TRUE)
    } else {
      cat("Decision: Fail to reject H0", fill = TRUE)
    }
    cat("--------------------------------------------", fill = TRUE)
  }
}
