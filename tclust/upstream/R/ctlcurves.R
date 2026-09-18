## Fix drop[i,j], out[i,j]  and unrestr.fact[i,j]
##

##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#'
#' Classification Trimmed Likelihood Curves
#' 
#' @name ctlcurves
#' @aliases print.ctlcurves
#' @description The function applies \code{\link{tclust}} several times on a given dataset while parameters 
#'  \code{alpha} and \code{k} are altered. The resulting object gives an idea of the optimal 
#'  trimming level and number of clusters considering a particular dataset.
#' @details These curves show the values of the trimmed classification (log-)likelihoods 
#'  when altering the trimming proportion \code{alpha} and the number of clusters \code{k}. 
#'  The careful examination of these curves provides valuable information for choosing 
#'  these parameters in a clustering problem. For instance, an appropriate \code{k} 
#'  to be chosen is one that we do not observe a clear increase in the trimmed classification 
#'  likelihood curve for k with respect to the k+1 curve for almost all the range 
#'  of alpha values. Moreover, an appropriate choice of parameter alpha may be derived 
#'  by determining where an initial fast increase of the trimmed classification 
#'  likelihood curve stops for the final chosen k. A more detailed explanation can 
#'  be found in García-Escudero et al. (2011). 
#' @param x A matrix or data frame of dimension n x p, containing the observations (row-wise).
#' @param k A vector of cluster numbers to be checked. By default cluster numbers from 1 to 5 are examined.
#' @param alpha A vector containing the alpha levels to be checked. By default \code{alpha} 
#'  levels from 0 to 0.2 (continuously increased by 0.01), are checked.
#' @param restr.fact The restriction factor passed to \code{\link{tclust}}.
#' @param parallel A logical value, to be passed further to \code{tclust()}.
#' @param \ldots Further arguments (as e.g. \code{restr}), passed to \code{\link{tclust}} 
#' @param trace Defines the tracing level, which is set to \code{1} by default. 
#'  Tracing level \code{2} gives additional information on the current iteration.
#'  Tracing level \code{} suppresses all trace messages.
#' @return The function returns an S3 object of type \code{ctlcurves} containing the following components:
#'	\itemize{
#'  \item \code{par} A list containing all the parameters passed to this function 
#'	\item \code{obj}  An array containing the objective functions values of each computed cluster-solution 
#'	\item \code{min.weights} An array containing the minimum cluster weight of each computed cluster-solution 
#'  }
#' @references 
#'    \enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2011), 
#'    "Exploring the number of groups in robust model-based clustering." \emph{Statistics and Computing}, \bold{21} 
#'    pp. 585-599, <doi:10.1007/s11222-010-9194-z>
#'
#' @examples
#'
#' \donttest{
#'
#'  #--- EXAMPLE 1 ------------------------------------------
#'
#'  sig <- diag (2)
#'  cen <- rep (1, 2)
#'  x <- rbind(MASS::mvrnorm(108, cen * 0,   sig),
#'  	       MASS::mvrnorm(162, cen * 5,   sig * 6 - 2),
#'  	       MASS::mvrnorm(30, cen * 2.5, sig * 50))
#'
#'  ctl <- ctlcurves(x, k = 1:4)
#'  ctl
#'
#'    ##  ctl-curves 
#'  plot(ctl)  ##  --> selecting k = 2, alpha = 0.08
#'
#'    ##  the selected model 
#'  plot(tclust(x, k = 2, alpha = 0.08, restr.fact = 7))
#'
#'  #--- EXAMPLE 2 ------------------------------------------
#'
#'  data(geyser2)
#'  ctl <- ctlcurves(geyser2, k = 1:5)
#'  ctl
#'  
#'    ##  ctl-curves 
#'  plot(ctl)  ##  --> selecting k = 3, alpha = 0.08
#'
#'    ##  the selected model
#'  plot(tclust(geyser2, k = 3, alpha = 0.08, restr.fact = 5))
#'
#'
#'  #--- EXAMPLE 3 ------------------------------------------
#'  
#'  data(swissbank)
#'  ctl <- ctlcurves(swissbank, k = 1:5, alpha = seq (0, 0.3, by = 0.025))
#'  ctl
#'  
#'    ##  ctl-curves 
#'  plot(ctl)  ##  --> selecting k = 2, alpha = 0.1
#'  
#'    ##  the selected model
#'  plot(tclust(swissbank, k = 2, alpha = 0.1, restr.fact = 50))
#'  
#' }
#'
ctlcurves <- function(x, k=1:4, alpha=seq(0, 0.2, len=6), restr.fact=50, parallel=FALSE, trace=1, ...) {

## disabled parameter:
  mah.alpha <- 0.05


  stopifnot (length (restr.fact) == 1)

  if(trace >= 1)
  {
    cat ("Depending on arguments x, k and alpha, ",
         "this function needs some time to compute.\n",
         "(Remove this message by setting \"trace = 0\")\n", sep = "")
     flush.console()
  }

  rs <- array (NA, c (length (k), length (alpha)))
  dimnames (rs) <- list (k = k, alpha = round (alpha, 2))
  ndoubt <- drop <- out <- obj <- min.weights <- unrestr.fact <- rs

#  stab <-
#  stab[,] <- TRUE

  p <- ncol (x)
  n.k <- length (k)
  n.alpha <- length (alpha)
  n <- nrow (x)

  for (i in 1:n.k) {
    for (j in 1:n.alpha)   {
      cur.alpha <- alpha[j]

      if (trace >= 2)
        cat ("k =", k[i], "; alpha =", alpha[j])

      clus <- tclust(x, k=k[i], alpha=alpha[j], restr.fact=restr.fact, parallel=parallel, trace=0, ...)

      out[i, j] <- sum(clus$mah > qchisq(1-mah.alpha, p), na.rm = TRUE)

      drop[i, j] <- FALSE   # clus$warnings$drop | clus$warnings$size | clus$warnings$sizep

      dsc <- DiscrFact(clus)
      ndoubt[i, j] <- sum (dsc$assignfact > dsc$threshold)

      unrestr.fact[i, j] <- clus$unrestr.fact
      obj[i,j] <- clus$obj                      ##  the objective function criterion
      min.weights[i,j] <- min(clus$weights)     ##  weights of the smallest group

      if (trace >= 2)
        cat ("; obj = ", clus$obj, ";min (weights) =", min(clus$weights), "\n")
    }
  }

  par <- list(x=x, k=k, alpha=alpha, mah.alpha=mah.alpha, restr.fact=restr.fact)

  ret <- list (obj=obj, min.weights=min.weights, par=par,
    unrestr.fact=unrestr.fact, ndoubt=ndoubt, out=out, drop=drop)

  class (ret) <- "ctlcurves"
  ret
}

#' @name plot.ctlcurves
#' @title The \code{plot} method for objects of class \code{ctlcurves}
#' @rdname plot.ctlcurves
#' @description The \code{plot} method for class \code{ctlcurves}: This function implements 
#'  a series of plots, which display characteristic values 
#'  of the each model, computed with different values for \code{k} and \code{alpha}.
#' @param x The ctlcurves object to be shown
#' @param what A string indicating which type of plot shall be drawn. See the details section for more information.
#' @param main A character-string containing the title of the plot.
#' @param xlab,ylab,xlim,ylim Arguments passed to plot().
#' @param col A single value or vector of line colors passed to \code{\link[graphics]{lines}}.
#' @param lty A single value or vector of line colors passed to \code{\link[graphics]{lines}}.
#' @param \ldots Arguments to be passed to or from other methods.
#' @details These curves show the values of the trimmed classification (log-)likelihoods 
#'  when altering the trimming proportion \code{alpha} and the number of clusters \code{k}.
#'  The careful examination of these curves provides valuable information for choosing these 
#'  parameters in a clustering problem. For instance, an appropriate \code{k} to be chosen 
#'  is one that we do not observe a clear increase in the trimmed classification likelihood
#'  curve for \code{k} with respect to the \code{k+1} curve for almost all the range of
#'  \code{alpha} values. Moreover, an appropriate choice of parameter \code{alpha} may
#'  be derived by determining where an initial fast increase of the trimmed classification 
#'  likelihood curve stops for the final chosen \code{k}. A more detailed explanation 
#'  can be found in García-Escudero et al. (2011).
#'	  
#'  This function implements a series of plots, which display characteristic values 
#'  of the each model, computed with different values for \code{k} and \code{alpha}.
#   The plot type is selected by setting argument \code{what} to one of the following values:
#'  \describe{
#'    \item{\code{"obj"}}{Objective function values.}
#'    \item{\code{"min.weights"}}{The minimum cluster weight found for each computed model. 
#'      This plot is intended to spot spurious clusters, which in 
#'	     general yield quite small weights.}
#'    \item{\code{"doubtful"}}{The number of "doubtful" decisions identified by \code{\link{DiscrFact}}.}
#'  }
#'
#' @references 
#'    \enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2011), 
#'    "Exploring the number of groups in robust model-based clustering." \emph{Statistics and Computing}, \bold{21} 
#'    pp. 585-599, <doi:10.1007/s11222-010-9194-z>
#'
#' @examples
#'
#'  #--- EXAMPLE 1 ------------------------------------------
#'
#'  \donttest{
#'  sig <- diag (2)
#'  cen <- rep (1, 2)
#'  x <- rbind(MASS::mvrnorm(108, cen * 0,   sig),
#'  	       MASS::mvrnorm(162, cen * 5,   sig * 6 - 2),
#'  	       MASS::mvrnorm(30, cen * 2.5, sig * 50))
#'
#'  (ctl <- ctlcurves(x, k = 1:4))
#'
#'  plot(ctl)
#'  }
#'
plot.ctlcurves <- function(x, what=c("obj", "min.weights", "doubtful"),
         main, xlab, ylab, xlim, ylim, col, lty=1, ...)
{
  mark.art = FALSE

  what <- match.arg(what[1], eval(formals ()$what))

  n.k <- length (x$par$k)
  idxb.use <- rep (TRUE, n.k)
  if(what == "obj") {
    dat <- x$obj
    set.ylab <- "Objective Function Value"
    set.main <- "CTL-Curves"
  } else if(what == "out") {
    dat <- x$out
    set.ylab <- "Number of Outliers among Clusters"
    set.main <- "Outlier-Curves"
  } else if(what == "doubtful") {
    dat <- x$ndoubt
    set.ylab <- "Number of Doubtful Decision"
    set.main <- "Doubtfull Decision-Curves"
  } else if(what == "min.weights") {
    idxb.use <- x$par$k != 1       # k == 1 -> min.weights == 1 -> we're not interested in that.
    dat <- x$min.weights
    set.ylab <- "Minimum Weigths"
    set.main <- "Minimum Weight-Curves"
  }

  if (missing (ylim)) ylim <- range (dat[idxb.use,])
  if (missing (xlim)) xlim <- range (x$par$alpha)
  if (missing (xlab)) xlab <- expression (alpha)
  if (missing (ylab)) ylab <- set.ylab
  if (missing (main)) main <- set.main

  n.alpha <- length (x$par$alpha)
  lty <- rep (lty, length (x$par$k))

  n.use.k <- sum (idxb.use)

  if (missing (col))
    col <- 1 + (1:n.k)
  else
    col <- rep (col, len = n.k)

  plot (0, type = "n", ylim = ylim, xlim = xlim ,
        main = main, xlab = xlab, ylab = ylab)
  mtext (paste ("Restriction Factor =", x$par$restr.fact), line = 0.25)

  pch <- rownames (dat)

#  idx.grey <- !x$stable
  if (mark.art)
    idx.grey <- x$par$restr.fact < x$unrestr.fact
  else
    idx.grey <- FALSE

  for (j in n.k:1)
  {
    if (!idxb.use[j])
    next

  cur.col <- rep (col[j], len = n.alpha)
  if (any (idx.grey))
      cur.col [idx.grey[j,]] <- 8

#    lines (x$par$alpha, dat[j,], type="b", col = cur.col, lty = lty[j],
#           pch = as.character (pch[j]))
    lines (x$par$alpha, dat[j,], type = "b", col = col[j], lty = lty[j],
           pch = "  ")

    text (x = x$par$alpha, y = dat[j,], col = cur.col, labels = pch[j])

  }

  if (what == "out")
    lines (x$par$alpha, nrow (x$par$x) * (1-x$par$alpha) * x$par$mah.alpha, lty = 2)
}

#' @export

print.ctlcurves <- function (x, ...) {

  cat ("Computed ", length (x$par$k) * length (x$par$alpha),
       " solutions (chosen restr.fact = ", x$par$restr.fact, ").\n\n",
     sep = "")

  idx.ar <- x$par$restr.fact < x$unrestr.fact

  if (any (idx.ar | x$drop))
  {
    uf <- array ("", dim = dim (x$unrestr.fact))
    uf[idx.ar] <- "*"
    uf[x$drop] <- paste (uf[x$drop], "k", sep = "")

    attributes (uf) <- attributes (x$unrestr.fact)

    print (uf, justify = "right", quote = FALSE)

    if (any (idx.ar))
        cat ("\n(*) Identified ", sum (idx.ar), " artificially restricted solutions.", sep = "")
    if (any (x$drop))
        cat ("\n(k) Identified ", sum (x$drop), " solutions with very small/dropped clusters.", sep = "")
  }
  else
    cat ("\nNo artificially restricted solutions or dropped clusters found.")

 cat ("\n")

 invisible(x)
}
