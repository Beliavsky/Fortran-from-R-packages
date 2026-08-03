### From Master Thesis of Emmanuel Profumo (w/ M.Maechler) Autumn 2016--March 2017
### "Generalized clusGap()" :  We cannot be 100% compatible to clusGap()

#' @param x the data, can be a data frame or a matrix
#' @param algo, a clustering algorithm function taking the prepared data and
#' a number of clusters as arguments
#' @param index, a function taking a clustering vector and the prepared data
#' which returns the value of a validity index. Index can also be a list of such
#' functions to obtain results for different indices.
#' For coherence with the originally proposed Gap Statistic in Tibshirani et al's
#' a LOWER value of the validity index implies a better clustering quality, so
#' indices such as average silhouette width should be added a minus sign.
#' This can be changed by setting the argument low=FALSE.
#' @param modelH0, a function which takes as argument at least the data x,
#' parameters estimated from the data, and further arguments in ...
#' @param K.max, number of different clusters for which the index should be
#' evaluated.
#' @param B, the number of bootstraps sample.
#' @param transformData, a function which takes the data x as argument and
#' processed it for clustering
#' @param modelH0Param, a function which takes as argument the data x and returns
#' a list of modelH0 parameters with matching names.
#' @param low, logical, if FALSE a HIGHER value of the index or the indices in the
#' user provided list implies a better clustering quality
#'
#' @return if index is just one function a list with components:
#' Indks, the bootstrap validity plots
#' Ind, the validity plot corresponding to the data
#' E.Ind, the sample mean of the bootstrap validity plots
#' gap, the calibrated validity plot, difference between E.Ind and Ind
#' gapHen, the gap divided by the standard deviation of the bootstraps validity plots
#' SE.sim, the sample standard error of the bootstrap validity plots with a correction
#' term for bootstrap estimation
#' SE, the sample standard error
#' if index is a list of index then each of the components above are lists with values
#' for each index
clusGapGen <- function(x, algo, index, modelH0, K.max, B = 100,
                       transformData = identity,
                       modelH0Param = function(y) list(),
                       low=TRUE, verbose = interactive(), ...)
{
  ind.isList <- is.list(index)
  if (is.function(index))
    index <- list(index)
  else if (!ind.isList || !all(vapply(index, is.function, NA)))
    stop("index has to be a function or a list of function")

  Ind <- E.Ind <- SE.sim <- SE <- index
  for (i in seq_along(index)) Ind[[i]] <- E.Ind[[i]] <- SE.sim[[i]] <- numeric(K.max)

  if(verbose) cat("Clustering k = 1,2,..., K.max (= ",K.max,"): .. ", sep='')
  xt <- transformData(x)
  for(k in 1:K.max){
    cls <- algo(xt,k)
    for (i in seq_along(index))
      Ind[[i]][k] <- index[[i]](cls,xt)
  }
  if(verbose) cat("done\n")
  Indks <- index
  for (i in seq_along(index)) Indks[[i]] <- matrix(0, B, K.max)

  param <- modelH0Param(x)
  if(verbose) cat("Bootstrapping, b = 1,2,..., B (= ", B,
                  ")  [one \".\" per sample]:\n", sep="")
  for (b in 1:B) {
    z <- do.call(modelH0,c(list(x=x),param,list(...)))
    zt <- transformData(z)
    for(k in 1:K.max) {
      cls <- algo(zt,k)
      for (i in seq_along(index))
        Indks[[i]][b,k] <- index[[i]](cls,zt)
    }
    if(verbose) cat(".", if(b %% 50 == 0) paste(b,"\n"))
  }
  if(verbose && (B %% 50 != 0)) cat("",B,"\n")
  gap <- gapHen <- index
  for (i in seq_along(index)){
    E.Ind[[i]] <- colMeans(Indks[[i]])
    var.i <- apply(Indks[[i]], 2, var)
    SE[[i]] <- sqrt(var.i)
    SE.sim[[i]] <- sqrt((1 + 1/B) * var.i)
    gap[[i]] <- gap.i <- E.Ind[[i]] - Ind[[i]]
    gapHen[[i]] <- gap.i/SE[[i]]
    if (!low) {
      gap[[i]] <- -gap[[i]]
      gapHen[[i]] <- -gapHen[[i]]
    }
  }
  ## TODO: really? make distinction of *list* of indices vs 1 index?
  ## ---   well maybe, keep it: *The* usual case = _one_ index (or not?)
  if (ind.isList) {
    list(Indks=Indks,Ind=Ind, E.Ind=E.Ind, gap = gap,
         gapHen = gapHen ,SE.sim=SE.sim,SE=SE)
  }
  else {list(Indks=Indks[[1]],Ind=Ind[[1]], E.Ind=E.Ind[[1]], gap = gap[[1]],
             gapHen = gapHen[[1]] ,SE.sim=SE.sim[[1]],SE=SE[[1]])
  }
}


#' @param clusGapRes, a list returned by a call of function clusGapGen
#' @param main, the main title to the plots
#' @param divBySd, logical, if TRUE plot for the standardize version of the gap

clusGapGen.plot <- function(clusGapRes,divBySd=FALSE,main=""){


  if (!is.list(clusGapRes$Ind)) clusGapRes <- lapply(clusGapRes,
                                                     function(el) list(el))
  for (i in seq_along(clusGapRes$Ind)){
    B <- nrow(clusGapRes$Indks[[i]])
    K.max <- ncol(clusGapRes$Indks[[i]])
    std <- t(replicate(B,rep(1,K.max)))
    ylm <- range(rbind(clusGapRes$Indks[[i]],clusGapRes$Ind[[i]]),na.rm=TRUE)
    ylb <- names(clusGapRes$Ind[i])
    if (is.null(ylb)) ylb <- paste("Index",as.character(i))
    gp <- "gap"
    namegap <- paste(gp,ylb)
    if (divBySd) {
      std <- t(replicate(B,clusGapRes$SE.sim[[i]]))
      gp <- "gapHen"
      namegap <- paste(gp,ylb)
    }
    matplot(replicate(B,1:K.max),t(clusGapRes$Indks[[i]]),
            pch = "-",xlab = "k",ylab = ylb,
            type="l",ylim=ylm,main=main
    )

    lines(1:K.max,clusGapRes$E.Ind[[i]],type="l",col="white",lwd=2)
    lines(1:K.max,clusGapRes$Ind[[i]],lwd=2)

    boxplot((clusGapRes$Indks[[i]]-t(replicate(B,clusGapRes$Ind[[i]])))/std,
            pch = "*",xlab = "k", ylab = namegap,type="l",col=c("light blue"),
            notch=TRUE, border="grey",main=main
    )

    lines(1:K.max,clusGapRes[[gp]][[i]],type="l",xlab="k",ylab="",
          col="orangered",lwd=1.5)

  }
}
