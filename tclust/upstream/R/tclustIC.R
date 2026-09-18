##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#'
#' Performs cluster analysis by calling \code{\link{tclust}} for different
#'  number of groups \code{k} and restriction factors \code{c}
#' 
#' @name tclustIC
#' @aliases print.tclustIC
#' @description Computes the values of BIC (MIXMIX), ICL (MIXCLA) or CLA (CLACLA),
#'  for different values of \code{k} (number of groups) and different values of \code{c}
#'  (restriction factor), for a prespecified level of trimming (the last two letters in the name
#'  stand for 'Information Criterion'). 
#'
#' @param x A matrix or data frame of dimension n x p, containing the observations (row-wise).
#' @param kk an integer vector specifying the number of mixture components (clusters) 
#'  for which the information criteria are be calculated. By default \code{kk=1:5}.
#' @param cc an  vector specifying the values of the restriction factor which have to 
#'  be considered. By default \code{cc=c(1, 2, 4, 8, 16, 32, 64, 128)}.
#' @param whichIC A character value which specifies which information criteria must be computed
#'  for each \code{k} (number of groups) and each value of the restriction factor 
#'  \code{c}. Possible values for \code{whichIC} are:
#'  \itemize{
#'   \item "MIXMIX": a mixture model is fitted and for computing the information criterion
#'      the mixture likelihood is used. This option corresponds to the use of the Bayesian
#'      Information criterion (BIC). In output just the matrix \code{MIXMIX} is given.
#'  \item "MIXCLA": a mixture model is fitted but to compute the information criterion
#'      the classification likelihood is used. This option corresponds to the use of the
#'      Integrated Complete Likelihood (ICL). In the output just the matrix \code{MIXCLA} is given.
#'  \item "CLACLA": everything is based on the classification likelihood. This information
#'      criterion will be called CLA. In the output just the matrix \code{CLACLA} is given.
#'  \item "ALL": both classification and mixture likelihood are used. In this case all
#'      three information criteria CLA, ICL and BIC are computed. In the output all
#'      three matrices \code{MIXMIX}, \code{MIXCLA} and \code{CLACLA} are given.
#'  }
#'
#' @param alpha The proportion of observations to be trimmed.
#' @param parallel A logical value, specifying whether the calls to \code{\link{tclust}} should be done in parallel.
#' @param n.cores The number of cores to use when paralellizing, only taken into account if \code{parallel=TRUE}.
#' @param trace Whether to print intermediate results. Default is \code{trace=FALSE}.
#' @param \ldots Further arguments (as e.g. \code{restr}), passed to \code{\link{tclust}} 
#' @return The functions \code{print()} and \code{summary()} are used to obtain and print a
#'  summary of the results. The function returns an S3 object of type 
#'  \code{tclustIC} containing the following components:
#'	\itemize{
#'  \item call the matched call
#'  \item kk a vector containing the values of \code{k} (number of components) which have been considered.
#'      This vector is identical to the optional argument \code{kk} (default is \code{kk=1:5}.
#'  \item cc a vector containing the values of \code{c} (values of the restriction factor) which
#'      have been considered. This vector is identical to the optional argument \code{cc} (defalt is \code{cc=c(1, 2, 4, 8, 16, 32, 64, 128)}.
#'  \item alpha trimming level
#'  \item whichIC Information criteria used 
#'  \item CLACLA a matrix of size \code{length(kk)-times-length(cc)} containinig the value of
#'      the penalized classification likelihood. This output is present only if \code{whichIC="CLACLA"} or  \code{whichIC="ALL"}.
#'  \item IDXCLA a matrix of lists of size \code{length(kk)-times-length(cc)} containinig the assignment of each unit
#'      using the classification model. This output is present only if \code{whichIC="CLACLA"} or  \code{whichIC="ALL"}.
#'  \item MIXMIX a matrix of size \code{length(kk)-times-length(cc)} containinig the value of
#'      the penalized mixtrue likelihood. This output is present only if \code{whichIC="MIXMIX"} or  \code{whichIC="ALL"}.
#'  \item IDXMIX a matrix of lists of size \code{length(kk)-times-length(cc)} containinig the assignment of each unit
#'      using the classification model. This output is present only if \code{whichIC="MIXMIX"} or  \code{whichIC="ALL"}.
#'  \item MIXCLA a matrix of size \code{length(kk)-times-length(cc)} containinig the value of
#'      the ICL criterion. This output is present only if \code{whichIC="MIXCLA"} or  \code{whichIC="ALL"}.
#'  \item x the input data matrix of size \code{length(n)-times-length(p)} with which 
#'      the Information Criteria were computed.
#'  }
#' @references
#'      Cerioli, A., Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2017).
#'      Finding the Number of Groups in Model-Based Clustering via Constrained Likelihoods,
#'      \emph{Journal of Computational and Graphical Statistics}, pp. 404-416,
#'      https://doi.org/10.1080/10618600.2017.1390469.
#'
#' @seealso \code{\link{tclust}}
#' @examples
#'
#'  #--- EXAMPLE 1 ------------------------------------------
#'  \donttest{
#'  data(geyser2)
#'  (out <- tclustIC(geyser2, whichIC="MIXMIX", alpha=0.1))
#'  summary(out)
#'  ## Find the smallest value inside the table and write the corresponding
#'  ## values of k (number of groups) and c (restriction factor)
#'  inds <- which(out$MIXMIX == min(out$MIXMIX), arr.ind=TRUE)
#'  vals <- out$MIXMIX[inds]
#'  cat("\nThe smallest value of the IC is ", vals, 
#'      " and takes place for k=", out$kk[inds[1]], " and c=",   
#'      out$cc[inds[2]], "\n")
#'  }
#'
#'  #--- EXAMPLE 2 ------------------------------------------
#'  \donttest{
#'  data(flea)
#'  Y <- as.matrix(flea[, 1:(ncol(flea)-1)])    # select only the numeric variables
#'  rownames(Y) <- 1:nrow(Y)
#'  head(Y)
#'
#'  (out <- tclustIC(Y, whichIC="CLACLA", alpha=0.1))
#'  summary(out)
#'  ## Find the smallest value inside the table and write the corresponding
#'  ## values of k (number of groups) and c (restriction factor)
#'  inds <- which(out$CLACLA == min(out$CLACLA), arr.ind=TRUE)
#'  vals <- out$CLACLA[inds]
#'  cat("\nThe Smallest value of the IC is ", vals, 
#'      " and takes place for k=", out$kk[inds[1]], " and c=",   
#'      out$cc[inds[2]], "\n")
#'  }
#'
#'  #--- EXAMPLE 3 ------------------------------------------
#'  \donttest{
#'  data(swissbank)
#'  (out <- tclustIC(swissbank, whichIC="ALL"))
#'  
#'  plot(out)  ##  --> selecting k=3, c=128
#'  
#'  ##  the selected model
#'  plot(tclust(swissbank, k = 3, alpha = 0.1, restr.fact = 128))
#'  
#'  }
#'
tclustIC <- function(x, kk=1:5, cc=c(1, 2, 4, 8, 16, 32, 64, 128), alpha=0.05, 
    whichIC=c("ALL", "MIXMIX", "MIXCLA", "CLACLA"), parallel=FALSE, n.cores=-1, trace=FALSE, ...) {

    ## Check the dots ...
    args <- list(...)
    if("opt" %in% names(args))
        stop("Argument 'opt' cannot be used with tclustIC()!")
    if("restr.fact" %in% names(args))
        stop("Argument 'restr.fact' cannot be used with tclustIC()!")
    
    rlang::check_dots_used()
    
    whichIC <- match.arg(whichIC)
    if(trace) {
        cat ("Depending on arguments x, kk and cc, ",
             "this function needs some time to compute.\n",
             "(Remove this message by setting \"trace=FALSE\")\n", sep = "")
         flush.console()
    }

    n <- nrow (x)
    p <- ncol (x)
    nkk <- length(kk)
    ncc <- length(cc)
    aa <- cbind(rep(1:nkk, each=ncc), rep(1:ncc, nkk))

    if(trace){
        if(parallel) {
            cat("\n Running in parallel.", 
                "\nProgress bar will not display accurate information.\n")
        }
        pb <- txtProgressBar(min = 0, max = nkk*ncc, style = 3)
    }
    
    if(!parallel) {
        IClist <- vector("list", length=nrow(aa))
        for(ii in 1:nrow(aa)) {
            ik <- aa[ii, 1]; ic <- aa[ii, 2]
            
            ##  if(trace)
            ##      cat("k =", kk[ik], "; c=", cc[ic], "\n")
            if(trace) {
                setTxtProgressBar(pb, ii)
                flush.console()
            }
            
            ax <- list(ik=ik, ic=ic, CLACLA=NULL, MIXMIX=NULL, MIXCLA=NULL, IDXCLA=NULL, IDXMIX=NULL)
            if(whichIC == "CLACLA" || whichIC == "ALL") {        
                out <- tclust(x, k=kk[ik], alpha=alpha, restr.fact=cc[ic], parallel=FALSE, trace=0, opt="HARD", ...)
                ax$CLACLA <- out$CLACLA
                ax$IDXCLA <- out$cluster
            } 
            if(whichIC == "MIXMIX" || whichIC == "MIXCLA" || whichIC == "ALL") {
                out <- tclust(x, k=kk[ik], alpha=alpha, restr.fact=cc[ic], parallel=FALSE, trace=0, opt="MIXT", ...)
                
                if(whichIC == "MIXMIX" || whichIC == "ALL") {
                    ax$MIXMIX <- out$MIXMIX
                    ax$IDXMIX <- out$cluster
                }
                if(whichIC == "MIXCLA" || whichIC == "ALL") {
                    ax$MIXCLA <- out$MIXCLA
                }
            }
            IClist[[ii]] <- ax
            
        }
    } else {
                    
        ## Setup parallel cluster
        if(n.cores == -1) {
            n.cores <- detectCores()
        } else if(n.cores == -2){
            n.cores <- detectCores() - 1
        }
        
        parclus <- makeCluster(n.cores)
        registerDoParallel(parclus)
        

        comb <- function(results, x) { 
            if(trace) {
                setTxtProgressBar(pb, length(results))
                flush.console()
            }
            
            return(append(results, list(x)))
        }
  
        IClist <- foreach(ii = 1:nrow(aa),
                                .packages = "tclust",
                                .combine = "comb",
                                .init = list(),
                                ##.multicombine = TRUE,
                                .inorder = FALSE) %dopar% {

            ik <- aa[ii, 1]; ic <- aa[ii, 2]
          
            ax <- list(ik=ik, ic=ic, CLACLA=NULL, MIXMIX=NULL, MIXCLA=NULL, IDXCLA=NULL, IDXMIX=NULL)
            if(whichIC == "CLACLA" || whichIC == "ALL") {        
                out <- tclust(x, k=kk[ik], alpha=alpha, restr.fact=cc[ic], parallel=FALSE, trace=0, opt="HARD", ...)
                ax$CLACLA <- out$CLACLA
                ax$IDXCLA <- out$cluster
            } 
            if(whichIC == "MIXMIX" || whichIC == "MIXCLA" || whichIC == "ALL") {
                out <- tclust(x, k=kk[ik], alpha=alpha, restr.fact=cc[ic], parallel=FALSE, trace=0, opt="MIXT", ...)
                
                if(whichIC == "MIXMIX" || whichIC == "ALL") {
                    ax$MIXMIX <- out$MIXMIX
                    ax$IDXMIX <- out$cluster
                }
                if(whichIC == "MIXCLA" || whichIC == "ALL") {
                    ax$MIXCLA <- out$MIXCLA
                }
            }
            ax
        }
        stopCluster(parclus)
    }

    if(trace)
        close(pb)
    
    CLACLA <- matrix(NA, nrow=length(kk), ncol=length(cc))
    MIXMIX <- matrix(NA, nrow=length(kk), ncol=length(cc))
    MIXCLA <- matrix(NA, nrow=length(kk), ncol=length(cc))
    cla <- list()
    mix <- list()
    for(ii in 1:(nkk*ncc)){
        ik <- IClist[[ii]]$ik
        ic <- IClist[[ii]]$ic
        
        if(whichIC == "CLACLA" || whichIC == "ALL") {        
            CLACLA[ik, ic] <- IClist[[ii]]$CLACLA
            cla[[(ik-1)*ncc + ic]] <- IClist[[ii]]$IDXCLA
        }
        if(whichIC == "MIXMIX" || whichIC == "MIXCLA" || whichIC == "ALL") {
            if(whichIC == "MIXMIX" || whichIC == "ALL") {
                MIXMIX[ik, ic] <- IClist[[ii]]$MIXMIX
                mix[[(ik-1)*ncc + ic]] <- IClist[[ii]]$IDXMIX
            } 
            if(whichIC == "MIXCLA" || whichIC == "ALL") {
                MIXCLA[ik, ic] <- IClist[[ii]]$MIXCLA
            }
        }
    }   

    ret <- list(call=match.call(), kk=kk, cc=cc, alpha=alpha, whichIC=whichIC, x=x)
    xkk <- paste0("k=", kk)
    xcc <- paste0("c=", cc)

    if(whichIC == "CLACLA" || whichIC == "ALL") {
        ret$CLACLA <- CLACLA
        dimnames(ret$CLACLA) <- list(xkk, xcc)
        ret$IDXCLA <- matrix(cla, nrow=length(kk), byrow=TRUE)
        dimnames(ret$IDXCLA) <- list(xkk, xcc)
    }
    if(whichIC == "MIXMIX" || whichIC == "ALL") {
        ret$MIXMIX <- MIXMIX
        dimnames(ret$MIXMIX) <- list(xkk, xcc)
        ret$IDXMIX <- matrix(mix, nrow=length(kk), byrow=TRUE)
        dimnames(ret$IDXMIX) <- list(xkk, xcc)
    }
    if(whichIC == "MIXCLA" || whichIC == "ALL") {
        ret$MIXCLA <- MIXCLA
        dimnames(ret$MIXCLA) <- list(xkk, xcc)
    }
    
    class(ret) <- "tclustIC"
    
    ret
}

#' @export
print.tclustIC <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nInformation criteria for TCLUST:", x$whichIC, "\n", paste0("Trimming = ", x$alpha))
    cat("\nNumber of mixture components (clusters):", x$kk)
    cat("\nvalues of the restriction factor:", x$cc, "\n")
    if(!is.null(x$MIXMIX))
    {
        cat("\n\nPenalized mixture likelihood:\n")
        print.default(format(x$MIXMIX, digits = digits), print.gap = 2, quote = FALSE)
    }
    if(!is.null(x$CLACLA))
    {
        cat("\n\nPenalized classification likelihood:\n")
        print.default(format(x$CLACLA, digits = digits), print.gap = 2, quote = FALSE)
    }
    if(!is.null(x$MIXCLA))
    {
        cat("\n\nICL criterion:\n")
        print.default(format(x$MIXCLA, digits = digits), print.gap = 2, quote = FALSE)
    }

     invisible(x)
}

#' @export
summary.tclustIC <- function (object, ...) {
    ans <- list(tclustobj=object)
    class(ans) <- "summary.tclustIC"
    ans
}

#' @export
print.summary.tclustIC <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",
	paste(deparse(x$tclustobj$call), sep = "\n", collapse = "\n"), "\n", sep = "")

    cat("\nInformation criteria for TCLUST:", x$tclustobj$whichIC, "\n", paste0("Trimming = ", x$tclustobj$alpha))
    cat("\nNumber of mixture components (clusters):", x$tclustobj$kk)
    cat("\nValues of the restriction factor:", x$tclustobj$cc, "\n")

    invisible(x)
}

#' @name plot.tclustIC
#' @title The \code{plot} method for objects of class \code{tclustIC}
#' @rdname plot.tclustIC
#' @description The \code{plot} method for class \code{tclustIC}: This function implements 
#'  a series of plots, which display characteristic values 
#'  of each model, computed with different values for \code{k} and \code{c} for a fixed \code{alpha}.
#' @param x The \code{tclustIC} object to be shown
#' @param whichIC A string indicating which information criterion will be used. See the details section for more information.
#' @param cc choose which curves to plot (for which restriction factors \code{c}).
#'  If missing, by default all curves will be printed.
#' @param main A character-string containing the title of the plot.
#' @param xlab,ylab,xlim,ylim Arguments passed to plot().
#' @param col A single value or vector of line colors passed to \code{\link[graphics]{lines}}.
#' @param lty A single value or vector of line types passed to \code{\link[graphics]{lines}}.
#' @param \ldots Arguments to be passed to or from other methods.
#'
#' @references
#'      Cerioli, A., Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2017).
#'      Finding the Number of Groups in Model-Based Clustering via Constrained Likelihoods,
#'      \emph{Journal of Computational and Graphical Statistics}, pp. 404-416,
#'      https://doi.org/10.1080/10618600.2017.1390469.
#'
#' @examples
#'
#'  \donttest{
#'  sig <- diag (2)
#'  cen <- rep (1, 2)
#'  x <- rbind(MASS::mvrnorm(108, cen * 0,   sig),
#'  	       MASS::mvrnorm(162, cen * 5,   sig * 6 - 2),
#'  	       MASS::mvrnorm(30, cen * 2.5, sig * 50))
#'
#'  (out <- tclustIC(x, whichIC="ALL"))
#'
#'  plot(out)
#'  }
#'
plot.tclustIC <- function(x, whichIC, cc, main, xlab, ylab, xlim, ylim, col, lty, ...)
{
   
    nkk <- length(x$kk)

    if(missing(whichIC)) {
        whichIC <- x$whichIC
        if(whichIC == "ALL")
            whichIC <- "MIXMIX"  
    } else {
        if(x$whichIC == "ALL")
            choices <- c("MIXMIX", "MIXCLA", "CLACLA")
        else if(x$whichIC == "MIXMIX")
            choices <- c("MIXMIX")
        else if(x$whichIC == "MIXCLA")
            choices <- c("MIXCLA")
        else if(x$whichIC == "CLACLA")
            choices <- c("CLACLA")
            
        whichIC <- match.arg(whichIC, choices)
    }
    
    if(whichIC == "CLACLA") {
        dat <- x$CLACLA
        set.ylab <- "CLACLA"
        set.main <- "Penalized mixture likelihood"
    } else if(whichIC == "MIXMIX") {
        dat <- x$MIXMIX
        set.ylab <- "MIXMIX"
        set.main <- "Penalized mixture likelihood"
    } else if(whichIC == "MIXCLA") {
        dat <- x$MIXCLA
        set.ylab <- "MIXCLA"
        set.main <- "ICL criterion"
    } 
  

    if(missing(ylim)) ylim <- range(dat)
    if(missing(xlim)) xlim <- range(x$kk)
    if(missing(xlab)) xlab <- "Number of groups"
    if(missing(ylab)) ylab <- set.ylab
    if(missing(main)) main <- set.main
    
    ncc <- length(x$cc)
    
    if(missing(lty))
        lty <- 1 + (1:ncc)
    else
        lty <- rep(lty, len=ncc)
    
    if(missing(col))
        col <- 1 + (1:ncc)
    else
        col <- rep(col, len=ncc)
    
    pch <- 1 + (1:ncc)
    
    if(missing(cc))
    {
        plot(0, type="n", ylim=ylim, xlim=xlim, main=main, xlab=xlab, ylab=ylab)
        mtext(paste("Trimming =", x$alpha), line = 0.25)
        for(j in ncc:1) {    
            lines(x$kk, dat[,j], type = "b", col=col[j], lty=lty[j], pch=pch[j])   
        }
        
        legend("topright", legend=colnames(dat), lty=lty, pch=pch, col=col, ...)
    } else {    

        if(length(which(x$cc %in% cc)) != length(cc))
            stop(paste("All elements of cc must be in ", paste(x$cc, collapse=","), ".", sep=""))

        xcc <- if(length(cc) == 1) paste("c=", cc) else paste("c=[", paste(cc, collapse=","), "]", sep="")
        plot(0, type="n", ylim=ylim, xlim=xlim, main=xcc, xlab=xlab, ylab=ylab)
        for(j in which(x$cc %in% cc)) {    
            lines(x$kk, dat[,j], type = "b", col=col[j], lty=lty[j], pch=pch[j])   
        }
    }    
    
    invisible(x)
}

