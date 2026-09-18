## The parameters from a tclust object that we need:
##  - x, n, p, alpha, k, size, centers, cov
##  TO DO:
##  - x can be also object different from tclust: tkmeans or rlg?
##  - The optimization method used for computing x can be either HARD or MIXT
##
##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/rrdev/robClus", load_code=roxygen2:::load_installed)

#'
#' Discriminant Factor analysis for \code{tclust} objects
#' 
#' @name DiscrFact
#' @aliases print.DiscrFact
#' @description Analyzes a \code{tclust}-object by calculating discriminant factors 
#'  and comparing the quality of the actual cluster assignments to that of the second best 
#'  possible assignment for each observation. Cluster assignments of observations 
#'  with large discriminant factors are considered "doubtful" decisions. Silhouette 
#'  plots give a graphical overview of the discriminant factors distribution 
#'  (see \code{\link{plot.DiscrFact}}). More details can be found in García-Escudero et al. (2011).
#' @param x A \code{tclust} object.
#' @param threshold  A cluster assignment or a trimming decision for an observation with a 
#'  discriminant factor larger than \code{log(threshold)} is considered a "doubtful" decision. 
#' @return The function returns an S3 object of type \code{DiscrFact} containing the following components:
#'	\itemize{
#'      \item \code{x} A \code{tclust} object. 
#'      \item \code{ylimmin} A minimum y-limit calculated for plotting purposes. 
#'	    \item \code{ind} The actual cluster assignment. 
#'	    \item \code{ind2} The second most likely cluster assignment for each observation. 
#'	    \item \code{lik} The (weighted) likelihood of the actual cluster assignment of each observation. 
#'	    \item \code{lik2} The (weighted) likelihood of the second best cluster assignment of each observation. 
#'	    \item \code{assignfact} The factor \code{log(disc/disc2)}. 
#'	    \item \code{threshold} The threshold used for deciding whether \code{assignfact} indicates a "doubtful" assignment. 
#'	    \item \code{mean.DiscrFact} A vector of length \code{k + 1} containing the mean discriminant 
#'      factors for each cluster (including the outliers). 
#'  }
#' @references 
#'    \enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2011), 
#'    "Exploring the number of groups in robust model-based clustering." \emph{Statistics and Computing}, \bold{21} 
#'    pp. 585-599, <doi:10.1007/s11222-010-9194-z>
#'
#' 
DiscrFact <- function(x, threshold=1/10)
{
    if (!any (class (x) == "tclust"))
      stop ("parameter x: expected object of type \x22tclust\x22")

    p <- x$int$dim[2]
    alpha <- x$par$alpha
    n <- nrow(x$par$x)
    
    ll <- matrix(nrow=n, ncol=x$k)
    lik <- lik2 <- ind <- array (NA, n)

    no.trim = floor(n*(1-alpha))

    for(k in 1:x$k)
      ll[,k] <- (x$size[k] / no.trim) * dmvnrm(x$par$x, x$centers[,k], as.matrix(x$cov[,,k]))

    llo <- apply(-ll, 1, order)         # the row-wise order of matrix ll 
    if(!is.matrix(llo))
      llo <- matrix(llo, ncol=nrow(ll)) # -> transposed..

    lik <- ll[cbind (1:n, llo[1, ])]
    ind <- llo[1,]

    if(nrow(llo) >= 2) {
      lik2 <- ll[cbind(1:n, llo[2, ])]
      ind2 <- llo[2,]
    } else {
      ind2 <- rep (0, n)
      lik2 <- lik * threshold
     }

    mropt <- sort(lik)[n - floor (n * (1 - alpha)) + 1] 
    idx.out <- lik < mropt

    ind2[idx.out] <- ind[idx.out]
    lik2[idx.out] <- lik[idx.out]
    lik[idx.out] <- mropt
    ind[idx.out] <- 0

    assignfact <- log(lik2 / lik)
    idx.fin <- is.finite(assignfact)
    ylimmin <- min(assignfact[idx.fin]) * 1.5
    assignfact[!idx.fin] <- ylimmin * 2

    mean.DiscrFact <- array(,x$k)
    for (i in 0:x$k)
      mean.DiscrFact[i + 1] = mean(assignfact[ind == i])
      
    names(mean.DiscrFact) <- c ("O", 1:x$k)

    ret <- list(x=x, ylimmin=ylimmin, ind=ind, ind2=ind2, assignfact=assignfact, 
            lik=lik, lik2=lik2, threshold=log(threshold), 
            mean.DiscrFact=mean.DiscrFact)
    class (ret) <- "DiscrFact"
    ret
  }

#' @export

print.DiscrFact <- function (x, ...) {
  cat ("Mean overall discriminant factor:", mean (x$assignfact), "\n")
  cat ("Mean discriminant factor per cluster:\n")
  print (x$mean.DiscrFact)

  idx <- x$assignfact > x$threshold

  if(!sum(idx))
    cat("No decision is considered as doubtful\n")
  else
    cat(sum(idx), "decisions are considered as doubtful\n")

  invisible(x)
}

#' @name plot.DiscrFact
#' @title The \code{plot} method for objects of class \code{DiscrFact}
#' @rdname plot.DiscrFact
#' @description The \code{plot} method for class \code{DiscrFact}: Next to a plot of the \code{tclust} 
#'  object which has been used for creating the \code{DiscrFact} object, a silhouette plot 
#'  indicates the presence of groups with a large amount of doubtfully assigned 
#'  observations. A third plot similar to the standard \code{tclust} plot serves 
#'  to highlight the identified doubtful observations.
#'
#' @param x An object of class \code{DiscrFact} as returned from DiscrFact()
#' @param enum.plots A logical value indicating whether the plots shall be enumerated 
#'  in their title ("(a)", "(b)", "(c)").
#' @param xlab,ylab,xlim Arguments passed to funcion \code{plot.tclust()}
#' @param print.DiscrFact A logical value indicating whether each clusters mean discriminant factor shall be plotted
#' @param col.nodoubt Color of all observations not considered as to be assigned doubtfully.
#' @param by.cluster Logical value indicating whether optional parameters pch and col 
#'  (if present) refer to observations (FALSE) or clusters (TRUE)
#' @param \ldots Arguments to be passed to or from other methods
#'
#' @details \code{plot_DiscrFact_p2} displays a silhouette plot based on the discriminant 
#'  factors of the observations. A solution with many large discriminant factors is 
#'  not reliable. Such clusters can be identified with this silhouette plot. 
#'  Thus \code{plot_DiscrFact_p3} displays the dataset, highlighting observations with 
#'  discriminant factors greater than the given threshold. The function \code{plot.DiscrFact()} 
#'  combines the standard plot of a \code{tclust} object, and the two plots introduced here.
#'
#' @references 
#'    \enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2011), 
#'    "Exploring the number of groups in robust model-based clustering." \emph{Statistics and Computing}, \bold{21} 
#'    pp. 585-599, <doi:10.1007/s11222-010-9194-z>
#'
#' @examples
#'  sig <- diag (2)
#'  cen <- rep (1, 2)
#'  x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
#'  	       MASS::mvrnorm(540, cen * 5,   sig * 6 - 2),
#'  	       MASS::mvrnorm(100, cen * 2.5, sig * 50))
#'
#'  clus.1 <- tclust(x, k = 2, alpha=0.1, restr.fact=12)
#'  clus.2 <- tclust(x, k = 3, alpha=0.1, restr.fact=1)
#'
#'  dsc.1 <- DiscrFact(clus.1)
#'  plot(dsc.1)
#'
#'  dsc.2 <- DiscrFact(clus.2)
#'  plot(dsc.2)
#'
plot.DiscrFact <- function (x, enum.plots = FALSE, 
    xlab="Discriminant Factor", ylab="Clusters", print.DiscrFact=TRUE, xlim, 
    col.nodoubt= grey(0.8), by.cluster=FALSE, ...)
{
  if(enum.plots)
    main.pre <-  c ("(a)", "(b)", "(c)")
  else
    main.pre <- NULL

  old.par <- par (mfrow = c (1,3))

  plot(x$x, main.pre = main.pre[1], ...)
  plot_DiscrFact_p2 (x, main.pre = main.pre[2], xlab=xlab, ylab=ylab, print.DiscrFact=print.DiscrFact, xlim=xlim, ...)
  plot_DiscrFact_p3 (x, main.pre = main.pre[3], col.nodoubt=col.nodoubt, by.cluster=by.cluster, ...)  

  par (old.par)
}

plot_DiscrFact_p2 <- function(x, xlab="Discriminant Factor", ylab="Clusters", main, xlim,
          print.DiscrFact=TRUE, main.pre, ...)
{
  n <- x$x$int$dim[1]

  if(missing (main))
    main <- "Silhouette Plot"  
    
  if(!missing(main.pre) && !is.null(main.pre))
    main <- paste(main.pre, main)

  if(missing (xlim))
      xlim = c(x$ylimmin,0)
      
  plot(0, 0, xlim=xlim, ylim=c(1,n), type="n", xlab=xlab, ylab=ylab,
        main=main, axes=FALSE, ...)
  

  axis(side=1)

  cs <- c(0, cumsum(c(x$x$int$dim[1] - sum(x$x$size), x$x$size)))

  {
    ylines <- cs[-1]
    ylines <- ylines[-length (ylines)]
    abline (h = ylines, lty = 2)
  }

  cs = (cs[-1] + cs[-length(cs)]) / 2
  axis (side = 2, at = cs, labels = c ("O", 1:x$x$k))
  box ()

  cury = 0
  for (k in 0:x$x$k)
  {
    grupo.k <- sort (x$assignfact[x$ind == k])
    gs = length (grupo.k)
    if (gs > 0)
      polygon (c (0, grupo.k, 0), c (1, 1:gs, gs) + cury, border = 0, col =
              k + 1)
    
    {  
      ll <- cury
      ul <- cury + gs

      if (k == 0)
        ll <- par ("usr")[3]
      if (k == x$x$k)
        ul <- par ("usr")[4]

    }
    cury = cury + gs
  }

  xpos = sum (par ("usr")[1:2] * c (1, 4)) / 5
  
  if(print.DiscrFact)
    legend ("topleft", legend = format (x$mean.DiscrFact[(x$x$k + 1):1],
            digits = 4), inset = 0.04, col = 1 + (x$x$k:0), pch = 15,
            title = "Mean Discriminant Factors", box.lwd = 0, bty = "n")
  abline (v = x$threshold + 1, lty = 2)

}

plot_DiscrFact_p3 <- function(x, main="Doubtful Assignments", col, pch, 
          col.nodoubt= grey(0.8), by.cluster=FALSE, ...)
{
  idxplot <- x$assignfact > x$threshold
  n <- length(x$x$cluster)

  if(missing(col)) {
    col <- rep(col.nodoubt, n)
    col[idxplot] <- x$ind[idxplot] + 1 
  }

  if(missing (pch))
    pch <- x$x$cluster + 1
  
  plot(x$x, by.cluster=FALSE, col=col, pch=pch, main=main, sub="",sub1="", ...)
}


#' @name summary.DiscrFact
#' @title The \code{summary} method for objects of class \code{DiscrFact}
#' @rdname summary.DiscrFact
#' @description The summary method for class \code{DiscrFact}.
#'
#' @param object An object of class \code{DiscrFact} as returned from \code{DiscrFact()}.
#' @param hide.emtpy A logical value specifying whether clusters without doubtful 
#'  assignment shall be hidden.
#' @param show.clust A logical value specifying whether the number of doubtful 
#'  assignments per cluster shall be displayed.
#' @param show.alt A logical value specifying whether the alternative cluster 
#'  assignment shall be displayed.
#' @param \ldots Arguments passed to or from other methods.
#' @references 
#'    \enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2011), 
#'    "Exploring the number of groups in robust model-based clustering." \emph{Statistics and Computing}, \bold{21} 
#'    pp. 585-599, <doi:10.1007/s11222-010-9194-z>
#'
#' @examples
#'  sig <- diag (2)
#'  cen <- rep (1, 2)
#'  x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
#'  	       MASS::mvrnorm(540, cen * 5,   sig * 6 - 2),
#'  	       MASS::mvrnorm(100, cen * 2.5, sig * 50)
#'  )
#'
#'  clus.1 <- tclust(x, k = 2, alpha=0.1, restr.fact=12)
#'  clus.2 <- tclust(x, k = 3, alpha=0.1, restr.fact=1)
#'
#'  dsc.1 <- DiscrFact(clus.1)
#'  summary(dsc.1)
#'
#'  dsc.2 <- DiscrFact(clus.2)
#'  summary(dsc.2)
#'
summary.DiscrFact <- function(object, hide.emtpy=TRUE, show.clust, show.alt, ...) {
  idx <- object$assignfact > object$threshold
  k <- object$x$k
  
  if(missing(show.clust))
    show.clust <- k <= 20
  if(missing(show.alt))
    show.alt <- k <= 5

  k1 <- k + 1

  dn <- c ("O", 1:k)

  if(show.clust) {
    cat ("\nNumber of doubtful assignments in clusters:\n")
    facts1 <- object$ind[idx] + 1
    cha <- as.array(tabulate(facts1, nbins = k1))
    dimnames(cha) <- list(dn)
    if(hide.emtpy)
      cha <- cha[cha != 0]
    print(cha)
  }
  
  if(show.alt)  {
    facts2 <- object$ind[idx] * k1 + object$ind2[idx]  + 1

    cat("\nObservations with doubtful decision are alternatively assigned to clusters:\n")  
    chma <- matrix(tabulate(facts2, nbins = k1 * k1), ncol=k1, nrow=k1)
    dimnames(chma) <- list(dn, dn)
    if(hide.emtpy) {
      chma <- chma[apply(chma, 1, sum) != 0, , drop = FALSE]
      chma <- chma[, apply(chma, 2, sum) != 0, drop = FALSE]
    }
    print(chma)
  }  
}


