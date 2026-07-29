#' Print Method for the Package 'smoots'
#'
#'This function regulates how objects created by the package \code{smoots} are
#'printed.
#'
#'@param x an input object of class \code{smoots}.
#'@param ... included for compatibility; additional arguments will however
#'not affect the output.
#'
#'@export
#'
#'@return
#'None
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'}


# Print function for the R package 'smoots'-------------------------------
print.smoots <- function(x, ...) {


  if (attr(x, "function") == "rollCast") {
    cat("---------------------------", fill = TRUE)
    cat("| Results of the backtest |", fill = TRUE)
    cat("---------------------------", fill = TRUE)
    cat(" ", fill = TRUE)
    cat("Model: Semi-ARMA(", x[["model.par"]][["arma"]][[1]], ",",
      x[["model.par"]][["arma"]][[2]], "), bandwidth: ",
      sprintf("%.4f", x[["model.nonpar"]][["b0"]]), fill = TRUE,
      sep = "")
    if (x[["method"]] == "norm") {
      method <- "normal"
    } else {
      method <- "bootstrap"
    }
    if (x[["np.fcast"]] == "lin") {
      np.fcast <- "linear"
    } else {
      np.fcast <- "constant"
    }
    df1 <- data.frame(c(paste0(x[["alpha"]] * 100, "%"), method, np.fcast,
      sum(x[["breach"]]), sprintf("%.4f", x[["MASE"]]),
      sprintf("%.4f", x[["RMSSE"]])))
    colnames(df1) <- ""
    rownames(df1) <- c("Forecasting intervals:", "Method:", "Extrapolation:",
      "Breaches:", "MASE:", "RMSSE:")
    print.data.frame(df1, right = TRUE)
  }
  if (attr(x, "function") == "msmooth" | attr(x, "function") == "tsmooth") {
    if (attr(x, "method") == "lpr") {
      cat("-------------------------------------------------", fill = TRUE)
      cat("| Results of the nonparametric trend estimation |", fill = TRUE)
      cat("-------------------------------------------------", fill = TRUE)
      cat("Method: Local Polynomial Regression", fill = TRUE)
      result_vector <- c(as.character(x$n), x$niterations,
        sprintf("%.4f", x$b0))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Number of observations:",
        "Iterations until convergence:", "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart, x$p, x$mu, x$Mcf, x$bvc, x$InfR,
                   x$bb, x$cb)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value:", "Order of polynomial:",
                      "Smoothness parameter:",
                      "Variance factor estimation:",
                      "Enlarged bandwidth (var. factor):",
                      "Inflation rate:", "Boundary method:",
                      "Boundary cut-off:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "ws", "b0", "cf0")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
        "Weights:", "Optimal bandwidth:", "Estimated variance factor:")
      print.data.frame(abbr, right = FALSE)
      cat(" ", fill = TRUE)
      cat("Iterations:", fill = TRUE)
      cat("-----------", fill = TRUE)
      if (x$niterations < 10) {
        it.names <- paste0("i = ", 1:x$niterations)
      } else {
        it.names <- paste0("i = ", sprintf("%2.f", 1:x$niterations))
      }
      print.data.frame(data.frame(bandwidth = sprintf("%.4f", x$iterations),
        row.names = it.names))

    } else if (attr(x, "method") == "kr") {
      cat("-------------------------------------------------", fill = TRUE)
      cat("| Results of the nonparametric trend estimation |", fill = TRUE)
      cat("-------------------------------------------------", fill = TRUE)
      cat("Method: Kernel Regression", fill = TRUE)
      result_vector <- c(as.character(x$n), x$niterations,
        sprintf("%.4f", x$b0))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Number of observations:",
        "Iterations until convergence:", "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart, x$p, x$mu, x$Mcf, x$bvc, x$InfR,
                   x$bb, x$cb)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value:", "Order of polynomial:",
                      "Smoothness parameter:",
                      "Variance factor estimation:",
                      "Enlarged bandwidth (var. factor):",
                      "Inflation rate:", "Boundary method:",
                      "Boundary cut-off:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df, right = FALSE)

        cat("", fill = TRUE)
        cat("Components of the object ($):", fill = TRUE)
        cat("-----------------------------")
        abbreviations <- c("ye", "orig", "res", "b0", "cf0")
        abbr <- data.frame(abbreviations)
        colnames(abbr) <- ""
        rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                            "Optimal bandwidth:", "Estimated variance factor:")
        print.data.frame(abbr, right = FALSE)
        cat(" ", fill = TRUE)
        cat("Iterations:", fill = TRUE)
        cat("-----------", fill = TRUE)
        if (x$niterations < 10) {
          it.names <- paste0("i = ", 1:x$niterations)
        } else {
          it.names <- paste0("i = ", sprintf("%2.f", 1:x$niterations))
        }
        print.data.frame(data.frame(bandwidth = sprintf("%.4f", x$iterations),
          row.names = it.names))

    }
  } else if(attr(x, "function") == "dsmooth") {
      cat("------------------------------------------------------",
        fill = TRUE)
      cat("| Results of the nonparametric derivative estimation |",
        fill = TRUE)
      cat("------------------------------------------------------",
        fill = TRUE)
      cat("Method: Local Polynomial Regression", fill = TRUE)
      result_vector <- c(as.character(x[["v"]]), x$n, x$niterations,
        sprintf("%.4f", x$b0))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Order of derivative:", "Number of observations:",
        "Iterations until convergence:", "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart.p, x$pp, x$bStart, x$p, x$mu, x$InfR)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value (var. factor):",
                      "Order of polynomial (var. factor):",
                      "Bandwidth starting value (IPI):",
                      "Order of polynomial (IPI):",
                      "Smoothness parameter:",
                      "Inflation rate:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "ws", "b0", "cf0")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Weights:",
                          "Optimal bandwidth:", "Estimated variance factor:")
      print.data.frame(abbr, right = FALSE)
      cat(" ", fill = TRUE)
      cat("Iterations:", fill = TRUE)
      cat("-----------", fill = TRUE)
      if (x$niterations < 10) {
        it.names <- paste0("i = ", 1:x$niterations)
      } else {
        it.names <- paste0("i = ", sprintf("%2.f", 1:x$niterations))
      }
      print.data.frame(data.frame(bandwidth = sprintf("%.4f", x$iterations),
        row.names = it.names))

  } else if (attr(x, "function") == "gsmooth" |
             attr(x, "function") == "knsmooth") {
    if (x$bb == 0) {
      k_nearest <- "N"
    } else if (x$bb == 1) {
      k_nearest <- "Y"
    }
    cat("-------------------------------------------", fill = TRUE)
    cat("| Results of the nonparametric estimation |", fill = TRUE)
    cat("-------------------------------------------", fill = TRUE)
    if (attr(x, "function") == "gsmooth") {
      cat("Method: Local Polynomial Regression", fill = TRUE)
      result_vec <- c(x$n, x$v, x$p, x$mu, x$b, k_nearest)
      result_name <- c("Number of observations:", "Order of derivative:",
                       "Order of polynomial:", "Smoothness parameter:",
                       "Bandwidth:", "k-nearest:")
      result_df <- data.frame(result_vec)
      colnames(result_df) <- ""
      rownames(result_df) <- result_name
      print.data.frame(result_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "b", "ws")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                          "Bandwidth:", "Weighting system:")
      print.data.frame(abbr, right = FALSE)

    } else if (attr(x, "function") == "knsmooth") {
      cat("Method: Kernel Regression", fill = TRUE)
      result_vec <- c(x$n, 0, x$mu, x$b, k_nearest)
      result_name <- c("Number of observations:", "Order of derivative:",
                       "Smoothness parameter:", "Bandwidth:", "k-nearest:")
      result_df <- data.frame(result_vec)
      colnames(result_df) <- ""
      rownames(result_df) <- result_name
      print.data.frame(result_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "b")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                          "Bandwidth:")
      print.data.frame(abbr, right = FALSE)
    }
  } else if (attr(x, "function") == "confBounds") {
    cat("-----------------------------------------------", fill = TRUE)
    cat("| Results of the confidence bounds estimation |", fill = TRUE)
    cat("-----------------------------------------------", fill = TRUE)
    result <- c(as.character(x$n), x$v, sprintf("%.4f", x$b.ub))
    result_df <- data.frame(result)
    colnames(result_df) <- ""
    rownames(result_df) <- c("Number of observations:", "Order of derivative:",
      "Adjusted bandwidth:")
    print.data.frame(result_df)
  }


}

