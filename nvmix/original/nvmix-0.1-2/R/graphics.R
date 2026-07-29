### QQ-plot of Mahalanobis distances for visual GOF ############################



#' Perform univariate goodness-of-fit test
#'
#' @param maha2 data vector (typically, squared maha distance)
#' @param theodf a function of one argument; theoretical df against which test
#'        is performed
#' @param test character, one of c("KS.AD", "KS", "AD", "none")
#'        where "KS" = Kolm. Smirn. and "AD" = And. Darl.
#' @return either NULL, otherwise a list with one or two elements,
#'         each of which is the result returned by ks.test() and the other is
#'         ad.test(). The attribute gives a displayable string for the legend
#' @author Erik Hintz

gof_test <- function(maha2, theodf, test = c("KS.AD", "KS", "AD", "none")){
   test <- match.arg(test)
   res <- if(test == "none") NULL else if(test == "KS"){
      KS.out <- ks.test(maha2, y = theodf)
      KS.out$method <- "Kolm. Smirn. GoF test" # shorter
      list(KS.out = KS.out)
   } else {
      AD.out <- ad.test(maha2, distr.fun = theodf)
      AD.out$method <- "And. Darl. GoF test"
      if(test == "KS.AD"){
         KS.out <- ks.test(maha2, y = theodf)
         KS.out$method <- "Kolm. Smirn. GoF test" # shorter
         list(KS.out = KS.out, AD.out = AD.out)
      } else {
         list(AD.out = AD.out)
      }
   }
   res
}



##' @title QQ Plot of Mahalanobis distances versus their theoretical quantiles
##' @param x (n,d) data matrix
##' @param qmix see ?pnvmix()
##' @param loc see ?pnvmix()
##' @param scale see ?pnvmix()
##' @param fitnvmix_object object of class 'fitnvmix'; if provided, x, qmix,
##'        loc and scale are ignored.
##' @param trafo.to.normal logical, if TRUE the
##'        underlying Mahalanobis distances are mapped to normals by a probability-
##'        quantile-transform so that the resulting QQ plot is essentially a normal
##'        QQ plot.
##' @param test character; specifying if (and which) GoF test shall be performed.
##' @param boot.pars  list with elements 'B' (Bootstrap sample size for computing
##'        CIs) and 'level' (confidence level).
##' @param plot logical if the result should be plotted.
##' @param verbose see ?pnvmix()
##' @param control see ?pnvmix()
##' @param digits number of digits for the test statistic and the p-value
##' @param plot.pars list; see ?get_set_qqplot_param()
##' @param ... see ?pnvmix()
##' @return invisibly returns an object of class 'qqplot_maha'.
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
qqplot_maha <- function(x, qmix, loc, scale, fitnvmix_object,
                        trafo.to.normal = FALSE, test = c("KS.AD", "KS", "AD",
                                                          "none"),
                        boot.pars = list(B = 500, level = 0.95),
                        plot = TRUE, verbose = TRUE, control = list(),
                        digits = max(3, getOption("digits") - 4),
                        plot.pars = list(), ...){
   ## Initialize and check inputs
   control <- get_set_param(control)
   control$newton.df.reltol <- control$qqplot.df.reltol
   test <- match.arg(test)
   plot.pars <- get_set_qqplot_param(plot.pars)
   call <- match.call() # for return
   if(hasArg(fitnvmix_object)){
      stopifnot(inherits(fitnvmix_object, "fitnvmix"))
      ## Grab elements from the fitnvmix_object and then call recursively
      return(qqplot_maha(
         x = fitnvmix_object$data, qmix = fitnvmix_object$qmix,
         loc = fitnvmix_object$loc, scale = fitnvmix_object$scale,
         trafo.to.normal = trafo.to.normal, test = test,
         boot.pars = boot.pars, plot = plot, verbose = verbose,
         control = control, digits = digits, plot.pars = plot.pars,
         nu = fitnvmix_object[[1]]))
      ## Note: 'fitnvmix_object$qmix' is either a function(u, nu) or a string
      ## ("pareto", "inverse.gamma", "constant"). In all cases, 'nu' can
      ## be as the name of the argument supplied via the ellipsis argument
   } else {
      if(!is.matrix(x)) x <- rbind(x)
      notNA <- rowSums(is.na(x)) == 0
      x     <- x[notNA,, drop = FALSE] # non-missing data (rows)
      n     <- nrow(x)
      d     <- ncol(x)
      ## Obtain sorted Mahalanobis distances (X-loc)^T scale^{-1} (X-loc)
      maha2 <- sort(mahalanobis(x, center = loc, cov = scale))
   }
   pp <- ppoints(n)
   ## Compute theoretical quantiles
   theo_quant <- if(trafo.to.normal){
      maha2 <- qnorm(pgammamix(maha2, qmix = qmix, d = d, control = control,
                               verbose = verbose, ...)) # maha2 mapped to N(0,1)
      x_ <- qnorm(c(0.25, 0.75))
      theodf <- function(x) pnorm(x)
      list(q = qnorm(pp), log.density = dnorm(pp, log = TRUE))
   } else {
      x_ <- qgammamix(c(0.25, 0.75), qmix = qmix, d = d, control = control,
                     verbose = verbose, q.only = TRUE, ...)
      theodf <- function(x) pgammamix(x, qmix = qmix, d = d, control = control,
                                      ...) # hypothesized cdf for below
      qgammamix(pp, qmix = qmix, d = d, control = control,
                verbose = verbose, q.only = FALSE, stored.values = NULL, ...)
   }
   ## Already compute intercept and slope as in 'qqline()'
   y_ <- quantile(maha2, c(0.25, 0.75))
   slope <- diff(y_) / diff(x_)
   int <- x_[1L] - slope * y_[1L]
   ## Obtain asymptotic CI (see Fox (2008), pp 35-36)
   logSE <- (log(pp) + log(1-pp) - log(n))/2 -
      theo_quant$log.density # logarithmic SE
   asymptSE <- exp(logSE)
   ## Obtain bootstrap CI
   boot_CI <- if(boot.pars$B > 1){
      alpha.2 <- (1 - boot.pars$level) / 2
      bootsample <- sapply(1:boot.pars$B, function(i) maha2[sort(
         sample(1:n, size = n, replace = TRUE))])
      apply(bootsample, 1, quantile, probs = c(alpha.2, 1 - alpha.2),
            names = FALSE) # (2, n) matrix
   } else NULL
   ## Compute test statistic
   testout <- gof_test(maha2, theodf = theodf, test = test)
   ## Create S3-class object of type 'qqplot_maha'
   out <- class_qqplot_maha(
      maha2 = maha2, theo_quant = theo_quant$q, boot_CI = boot_CI,
      trafo.to.normal = trafo.to.normal, asymptSE = asymptSE, test = test,
      testout = testout, int = int, slope = slope, B = boot.pars$B,
      call = call)
   ## Create plot
   if(plot) {
      plot(out, digits = digits, plot.pars = plot.pars)
   } else {
      invisible(out)
   }
}

### S3 class functions and methods #############################################

#' Function to define S3 class 'qqplot_maha'
#'
#' @param maha2 observed Mahalanobis distances
#' @param theo_quant theoretical quantiles
#' @param method string, "bootstrap" or "normaltrafo"
#' @param asymptSE asymptotic standard errors
#' @param test string ("KS", "AD" or "none")
#' @param testout list returned by ks.test() or ad.test() or NULL if test = "none"
#' @param int intercept of the fitted line
#' @param slope slope of the fitted line
#' @param call function call
#' @return S3 object of class 'qqplot_maha'
#' @author Erik Hintz
class_qqplot_maha <- function(maha2, theo_quant, boot_CI, trafo.to.normal, asymptSE,
                              test, testout, int, slope, B, call){
   res <- list(maha2 = maha2, theo_quant = theo_quant, boot_CI = boot_CI,
               trafo.to.normal = trafo.to.normal, asymptSE = asymptSE, test = test,
               testout = testout, int = int, slope = slope, B = B, call = call)
   ## Return object of class "qqplot_maha"
   structure(res, class = "qqplot_maha")
}

## Method 'print' for S3 class 'qqplot_maha'
print.qqplot_maha <- function(x, ..., digits = max(3, getOption("digits") - 4)){
   ## Print function call to qqplot_maha()
   cat("Call: ", deparse(x$call), "\n", sep = "")
   cat("\n")
   ## Provide sample size
   cat(paste0("Input: ", length(x$maha2), " squared Mahalanobis distances."), "\n")
   ## Check 'trafo.to.normal'
   if(x$trafo.to.normal)
      cat("Squared Mahalanobis distances were transformed to N(0, 1).", "\n", sep = "")
   cat("\n")
   ## Check if a test was performed
   test_text <- testtext(x, breaktext = TRUE, digits = digits)
   cat(test_text, "\n")
   cat("\n")
   cat("Computed results stored in the object:", "\n")
   cat("- theoretical quantiles in $theo_quant;", "\n")
   cat("- sorted, squared Mahalanobis distances in $maha2;", "\n")
   cat("- estimated, asymptotic standard errors in $asymptSE;", "\n")
   if(!is.null(x$boot_CI))
      cat(paste0("- Bootstrap CIs (estimated from ", x$B, " resamples) in $boot_CI;"), "\n")
   if(!is.null(x$testout))
      cat("- GoF test results in $testout;", "\n")
   invisible(x)
}


## Method 'plot' for S3 class 'qqplot_maha'
plot.qqplot_maha <-
   function(x, ..., digits = max(3, getOption("digits") - 4), plot.pars = list())
{
   ## Set parameters
   plot.pars <- get_set_qqplot_param(plot.pars)
   ## Construct the text for the axis
   axistext <- testtext(x, digits = digits)
   ## Plot
   plot(x$theo_quant, x$maha2, xlab = plot.pars$xlab, ylab = plot.pars$ylab,
        xlim = plot.pars$xlim, ylim = plot.pars$ylim, main = plot.pars$main,
        sub = plot.pars$sub, log = plot.pars$log, pch = plot.pars$pch,
        col = plot.pars$col[1])
   ## Add line
   if(plot.pars$plot_line)
      abline(x$int, x$slope, lty = plot.pars$lty[1], col = plot.pars$col[2],
             untf = TRUE)
   ## Add asymptotic CI
   for(i in c(-1, 1))
      lines(x$theo_quant, x$maha2 + i * x$asymptSE, lty = plot.pars$lty[2],
            col = plot.pars$col[3])
   lgnd <- "Asymptotic CI"
   numlegend <- 1 # number of elements in the legend
   ## Add bootstrap CI
   if(!is.null(x$boot_CI)){
      lgnd <- c(lgnd, "Bootstrap CI")
      numlegend <- 2
      for(i in 1:2)
         lines(x$theo_quant, x$boot_CI[i, ], lty = plot.pars$lty[3],
               col = plot.pars$col[4])

   }
   if(plot.pars$plot_legend) legend("topleft", lgnd,
                          lty = plot.pars$lty[2:(1+numlegend)],
                          col = plot.pars$col[3:(2+numlegend)],
                          bty = "n")
   if(plot.pars$plot_test) mtext(axistext, side = 4) # empty if no test performed
   ## Invisibly return input object
   invisible(x)
}


## Method testtext() for "qqplot_maha"
testtext <- function(x, breaktext, digits) {
   UseMethod("testtext") # name of the generic function
}

testtext.qqplot_maha <- function(x, breaktext = FALSE,
                                 digits = max(3, getOption("digits") - 4)) {
   if(is.null(x$testout)) "No GoF test performed." else {
      if(length(x$testout) == 2){
         if(breaktext){
            paste0("KS test: D = ", round(x$testout$KS.out$statistic, digits),
                   ", p = ", format(x$testout$KS.out$p.value, scientific = TRUE,
                                    digits = digits), "\n", "AD test: D = ",
                   round(x$testout$AD.out$statistic, digits),
                   ", p = ", format(x$testout$AD.out$p.value, scientific = TRUE,
                                    digits = digits), ".")
         } else {
            paste0("KS test/AD test: D = ", round(x$testout$KS.out$statistic, digits),
                   " / ", round(x$testout$AD.out$statistic, digits), " ; p = ",
                   format(x$testout$KS.out$p.value, scientific = TRUE,
                          digits = digits), " / ",
                   format(x$testout$AD.out$p.value, scientific = TRUE,
                          digits = digits), ".")
         }
      } else {
         stopifnot(length(x$testout) == 1)
         paste0(x$testout[[1]]$method, ": D = ", round(x$testout[[1]]$statistic, digits),
                ", p = ",
                format(x$testout[[1]]$p.value, scientific = TRUE, digits = digits),
                ".")

      }
   }

}
