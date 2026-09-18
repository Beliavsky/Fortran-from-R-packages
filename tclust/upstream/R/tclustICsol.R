##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#'
#' Extracts a set of best relevant solutions obtained by \code{\link{tclustIC}}
#' 
#' @name tclustICsol
#' @aliases print.tclustICsol
##' @description  The function \code{tclustICsol()} takes as input an object of class
#'  \code{\link{tclustIC}}, the output
#'  of function \code{\link{tclustIC}} (that is a series of matrices which contain
#'  the values of the information criteria BIC/ICL/CLA for different values of \code{k}
#'  and \code{c}) and extracts the first best solutions. Two solutions are considered
#'  equivalent if the value of the adjusted Rand index (or the adjusted Fowlkes and
#'  Mallows index) is above a certain threshold. For each tentative solution the program
#'  checks the adjacent values of \code{c} for which the solution is stable.
#'  A matrix with adjusted Rand indexes is given for the extracted solutions.
#'
#' @param obj An S3 object of class \code{\link{tclustIC}}
#'  (output of \code{\link{tclustIC}}) containing the values
#'  of the information criteria BIC (MIXMIX), ICL (MIXCLA) or CLA (CLACLA),
#'  for different values of \code{k} (number of groups) and different
#'  values of \code{c} (restriction factor), for a prespecified level of trimming.
#' @param whichIC A character value which Specifies the information criterion to 
#'  use to extract best solutions. Possible values for \code{whichIC} are:
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
#' @param nsol Number of best solutions to extract from BIC/ICL matrix.
#'  The default value of NumberOfBestSolutions is 5
#' @param index Index to use to compare partitions. If \code{index=Rand} (default) 
#'  the adjusted Rand index is used, else, \code{index="FM"}, the adjusted 
#'  Fowlkes and Mallows index is used
#' @param thresholdRI Threshold to identify spurious solutions - the threshold
#'  of the adjusted Rand index to use to consider two solutions as equivalent.
#' The default value of ThreshRandIndex is 0.7
#' @param trace Whether to print intermediate results. Default is \code{trace=FALSE}.
#' @return The function returns an S3 object of type \code{tclustICsol} containing the following components:
#'  \item{call}{the matched call}
#'  \item{kk}{a vector containing the values of \code{k} (number of components) which have been considered.
#'      This vector is identical to the optional argument \code{kk} (default is \code{kk=1:5}.}
#'  \item{cc}{a vector containing the values of \code{c} (values of the restriction factor) which
#'      have been considered. This vector is identical to the optional argument \code{cc}
#'   (defalt is \code{cc=c(1, 2, 4, 8, 16, 32, 64, 128)}.}
#'  \item{alpha}{trimming level}
#'  \item{whichIC}{Information criteria used}
#'  \item{MIXMIXbs}{a matrix of lists of size \code{NumberOfBestSolutions-times-5} which
#'  contains the details of the best solutions for MIXMIX (BIC). Each row refers to a
#'  solution. The information which is stored in the columns is as follows.
#'  \itemize{
#'  \item 1st col = value of k for which solution takes place
#'  \item 2nd col = value of c for which solution takes place;
#'  \item 3rd col = a vector of length \code{d} which contains the values of \code{c}
#'      for which the solution is uniformly better.
#'  \item 4th col = a vector of length \code{d + r} which contains the values of \code{c}
#'      for which the solution is considered stable (i.e. for which the value
#'      of the adjusted Rand index (or the adjusted Fowlkes and Mallows index)
#'      does not go below the threshold defined in input option \code{ThreshRandIndex}).
#'  \item 5th col = string which contains 'true' or 'spurious'. The solution is labelled
#'      spurious if the value of the adjusted Rand index with the previous solutions
#'      is greater than ThreshRandIndex.
#'  }
#'
#'  Remark: the field MIXMIXbs is present only if \code{whichIC=ALL} or \code{whichIC="MIXMIX"}.
#'  }
#'  \item{MIXMIXbsari}{a matrix of adjusted Rand indexes (or Fowlkes and Mallows indexes)
#'      associated with the best solutions for MIXMIX. A matrix of size \code{NumberOfBestSolutions-times-NumberOfBestSolutions}
#'      whose \code{i,j}-th entry contains the adjusted Rand index between classification produced by solution
#'      \code{i} and solution \code{j}, \code{i,j=1,2, ...,NumberOfBestSolutions}.
#'
#'  Remark: the field \code{MIXMIXbsari} is present only if \code{whichIC=ALL} or \code{whichIC="MIXMIX"}.
#'  }
#'  \item{ARIMIX}{a matrix of adjusted Rand indexes between two consecutive value of \code{c}.
#'      Matrix of size \code{k-by-length(cc)-1}. The first column contains the ARI indexes
#'      between \code{cc[2]} and \code{cc[1]} given \code{k}.
#'      The second column contains the the ARI indexes between \code{cc[3]} and \code{cc[2]} given \code{k}.
#'
#'  Remark: the field \code{ARIMIX} is present only if \code{whichIC=ALL} or \code{whichIC="MIXMIX"} or \code{whichIC="MIXCLA"}.
#'  }
#'  \item{MIXCLAbs}{has the same structure as \code{MIXMIXbs} but referres to MIXCLA.
#'
#'  Remark: the field MIXCLAbs is present only if \code{whichIC=ALL} or \code{whichIC="MIXCLA"}.
#'  }
#'  \item{MIXCLAbsari}{has the same structure as \code{MIXMIXbsari} but referres to MIXCLA.
#'
#'  Remark: the field \code{MIXMIXbsari} is present only if \code{whichIC=ALL} or \code{whichIC="MIXCLA"}.
#'  }
#'  \item{CLACLAbs}{has the same structure as \code{MIXMIXbs} but referres to CLACLA.
#'
#'  Remark: the field CLACLAbs is present only if \code{whichIC=ALL} or \code{whichIC="CLACLA"}.
#'  }
#'  \item{CLACLAbsari}{has the same structure as \code{MIXMIXbsari} but referres to CLACLA.
#'
#'  Remark: the field \code{CLACLAbsari} is present only if \code{whichIC=ALL} or \code{whichIC="CLACLA"}.
#'  }
#'  \item{ARICLA}{a matrix of adjusted Rand indexes between two consecutive value of \code{c}.
#'      Matrix of size \code{k-by-length(cc)-1}. The first column contains the ARI indexes
#'      between \code{cc[2]} and \code{cc[1]} given \code{k}.
#'      The second column contains the the ARI indexes between \code{cc[3]} and \code{cc[2]} given \code{k}.
#'
#'  Remark: the field \code{ARICLA} is present only if \code{whichIC=ALL} or \code{whichIC="CLACLA"}.
#'  }
#'  \item{x}{the input data matrix of size \code{n-times-p} with which
#'      the Information Criteria were computed.}
#'
#' @references
#'      Cerioli, A., Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2017).
#'      Finding the Number of Groups in Model-Based Clustering via Constrained Likelihoods,
#'      \emph{Journal of Computational and Graphical Statistics}, pp. 404--416,
#'      https://doi.org/10.1080/10618600.2017.1390469.
#'
#'  Hubert L. and Arabie P. (1985), Comparing Partitions, \emph{Journal of Classification}, 
#'      Vol. 2, pp. 193--218.
#'
#' @seealso \code{\link{tclust}}, \code{\link{tclustIC}}
#' @examples
#'
#'  #--- EXAMPLE 1 ------------------------------------------
#'  \donttest{
#'  data(geyser2)
#'  (out <- tclustIC(geyser2, whichIC="MIXMIX", alpha=0.1))
#'
#'  ##  Show the first two best solutions using as Information criterion MIXMIX
#'  cat("\nBest solutions using MIXMIX\n")
#'  outsol <- tclust::tclustICsol(out, whichIC="MIXMIX", nsol=2)
#'  print(outsol$MIXMIXbs)
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
#'  ## Find the smallest value inside the table and write the corresponding
#'  ## values of k (number of groups) and c (restriction factor)
#'  inds <- which(out$CLACLA == min(out$CLACLA), arr.ind=TRUE)
#'  vals <- out$CLACLA[inds]
#'  cat("\nThe Smallest value of the IC is ", vals, 
#'      " and takes place for k=", out$kk[inds[1]], " and c=",   
#'      out$cc[inds[2]], "\n")
#'
#'  ##  Show the first two best solutions using as Information criterion CLACLA
#'  cat("\nBest solutions using CLACLA\n")
#'  outsol <- tclust::tclustICsol(out, whichIC="CLACLA", nsol=2)
#'  print(outsol$CLACLAbs)
#'  
#'  }
#'
#'  #--- EXAMPLE 3 ------------------------------------------
#'  \donttest{
#'  data(swissbank)
#'  (out <- tclustIC(swissbank, whichIC="ALL"))
#'  
#'  outsol <- tclust::tclustICsol(out, whichIC="ALL", nsol=2)
#'  print(outsol$CLACLAbs)
#'  
#'  }
#'

tclustICsol <- function(obj, whichIC=c("ALL", "MIXMIX", "MIXCLA", "CLACLA"), 
        nsol=5, index=c("Rand", "FM"), thresholdRI=0.7, trace=FALSE) {

    if(!inherits(obj, "tclustIC"))
      stop ("parameter obj: expected object of type \x22tclustIC\x22")

    whichIC <- match.arg(whichIC)
    index <- match.arg(index)

    if(whichIC != obj$whichIC && obj$whichIC != "ALL") {        # Not OK combinations
        if(whichIC == "ALL")
            stop(paste0("The \x22tclustIC\x22 object \x22obj\x22 was computed with IC=\x22", obj$whichIC, "\x22!"))
        if(whichIC == "CLACLA" && obj$whichIC != "CLACLA")
            stop("The \x22tclustIC\x22 object \x22obj\x22 was not computed with IC=\x22CLACLA\x22!")
        if(whichIC != "CLACLA" && obj$whichIC == "CLACLA")
            stop("The \x22tclustIC\x22 object \x22obj\x22 was computed with IC=\x22CLACLA\x22!")
    }
    
    kk <- obj$kk
    cc <- obj$cc
    alpha <- obj$alpha
    nkk <- length(kk)
    ncc <- length(cc)

    ARIMIX <- ARICLA <- NULL
    if(whichIC != "CLACLA")                     
        ARIMIX <- matrix(1, nrow=nkk, ncol=ncc)
    if(whichIC == "CLACLA" || whichIC == "ALL")  
        ARICLA <- matrix(1, nrow=nkk, ncol=ncc)
    
    ## loop for different values of k (number of groups)    
    for(k in 1:nkk) {
        if(kk[k] == 1)   # do the computations for ARI or FM indexes when k is different from 1
            next
        
        ## loop through the different values of c
        for(j in 2:ncc) {
            if(whichIC != "CLACLA") {
                idjm1 <- obj$IDXMIX[[k, j-1]]
                idj <- obj$IDXMIX[[k, j]]
                ARIMIX[k, j] <- if(index == "Rand") randIndex(idjm1, idj)$AR else 
                                                    FowlkesMallowsIndex(idjm1, idj)$ABk
            }
            if(whichIC == "CLACLA" || whichIC == "ALL") {
                idjm1 <- obj$IDXCLA[[k, j-1]]
                idj <- obj$IDXCLA[[k, j]]
                ARICLA[k, j] <- if(index == "Rand") randIndex(idjm1, idj)$AR else 
                                                    FowlkesMallowsIndex(idjm1, idj)$ABk
            }
        }
    }

    ret <- list(call=match.call(), kk=kk, cc=cc, alpha=alpha, whichIC=whichIC, x=obj$x, ARIMIX=ARIMIX, ARICLA=ARICLA)

    if(whichIC == "MIXMIX" || whichIC == "ALL") {
    
        objbs <- findBestSolutions(pll=obj$MIXMIX, ARI=ARIMIX, IDX=obj$IDXMIX, 
            kk=kk, cc=cc, nsol=nsol, thresholdRI=thresholdRI, trace=trace)
        
        ret$MIXMIXbs <- objbs$Bestsols
        ret$MIXMIXbsari <- objbs$ARIbest
        
        ret$ARIMIX <- ARIMIX[, 2:ncol(ARIMIX)]

        ## Store matrix which contains in the columns the details of the
        ## classification             
        ##  - nrow(ret$MIXMIXbs) is the number of existing solutions - could be less than nsol             
        for(i in 1:nrow(ret$MIXMIXbs)) {
            kx <- which(kk == ret$MIXMIXbs[[i, 1]])
            cx <- which(cc == ret$MIXMIXbs[[i, 2]])
            IDX <- obj$IDXMIX[[kx, cx]]
            if(i == 1)
                ret$MIXMIXbsIDX <- matrix(NA, nrow=length(IDX), ncol=nsol)
            ret$MIXMIXbsIDX[, i] <- IDX    
        }
    }
   
    if(whichIC == "MIXCLA" || whichIC == "ALL") {
    
        objbs <- findBestSolutions(pll=obj$MIXCLA, ARI=ARIMIX, IDX=obj$IDXMIX, 
            kk=kk, cc=cc, nsol=nsol, thresholdRI=thresholdRI, trace=trace)
        
        ret$MIXCLAbs <- objbs$Bestsols
        ret$MIXCLAbsari <- objbs$ARIbest
        
        ret$ARIMIX <- ARIMIX[, 2:ncol(ARIMIX)]

        ## Store matrix which contains in the columns the details of the
        ## classification             
        ##  - nrow(ret$MIXCLAbs) is the number of existing solutions - could be less than nsol             
        for(i in 1:nrow(ret$MIXCLAbs)) {
            kx <- which(kk == ret$MIXCLAbs[[i, 1]])
            cx <- which(cc == ret$MIXCLAbs[[i, 2]])
            IDX <- obj$IDXMIX[[kx, cx]]
            if(i == 1)
                ret$MIXCLAbsIDX <- matrix(NA, nrow=length(IDX), ncol=nsol)
            ret$MIXCLAbsIDX[, i] <- IDX    
        }
    }
   
    if(whichIC == "CLACLA" || whichIC == "ALL") {
    
        objbs <- findBestSolutions(pll=obj$CLACLA, ARI=ARICLA, IDX=obj$IDXCLA, 
            kk=kk, cc=cc, nsol=nsol, thresholdRI=thresholdRI, trace=trace)

        ret$CLACLAbs <- objbs$Bestsols
        ret$CLACLAbsari <- objbs$ARIbest
        
        ret$ARICLA <- ARICLA[, 2:ncol(ARICLA)]

        ## Store matrix which contains in the columns the details of the
        ## classification.
        ##  - nrow(ret$CLACLAbs) is the number of existing solutions - could be less than nsol             
        for(i in 1:nrow(ret$CLACLAbs)) {
            kx <- which(kk == ret$CLACLAbs[[i, 1]])
            cx <- which(cc == ret$CLACLAbs[[i, 2]])
            IDX <- obj$IDXCLA[[kx, cx]]
            if(i == 1)
                ret$CLACLAbsIDX <- matrix(NA, nrow=length(IDX), ncol=nsol)
            ret$CLACLAbsIDX[, i] <- IDX    
        }
    }
   
    xkk <- paste0("k=", kk)
    xcc <- paste0("c", cc[2:length(cc)], "_vs_c", cc[1:(length(cc)-1)])

    if(!is.null(ret$ARIMIX))
        dimnames(ret$ARIMIX) <- list(xkk, xcc)
    if(!is.null(ret$ARICLA))
        dimnames(ret$ARICLA) <- list(xkk, xcc)

    if(!is.null(ret$MIXMIXbs))
        dimnames(ret$MIXMIXbs) <- list(1:nrow(ret$MIXMIXbs), c("k", "c", "c (uniformly better)", "c (stable)", "Solution"))
    if(!is.null(ret$MIXCLAbs))
        dimnames(ret$MIXCLAbs) <- list(1:nrow(ret$MIXCLAbs), c("k", "c", "c (uniformly better)", "c (stable)", "Solution"))
    if(!is.null(ret$CLACLAbs))
        dimnames(ret$CLACLAbs) <- list(1:nrow(ret$CLACLAbs), c("k", "c", "c (uniformly better)", "c (stable)", "Solution"))

    if(!is.null(ret$MIXMIXbsari))
        dimnames(ret$MIXMIXbsari) <- list(1:nrow(ret$MIXMIXbsari), 1:ncol(ret$MIXMIXbsari))
    if(!is.null(ret$MIXCLAbsari))
        dimnames(ret$MIXCLAbsari) <- list(1:nrow(ret$MIXCLAbsari), 1:ncol(ret$MIXCLAbsari))
    if(!is.null(ret$CLACLAbsari))
        dimnames(ret$CLACLAbsari) <- list(1:nrow(ret$CLACLAbsari), 1:ncol(ret$CLACLAbsari))

    if(!is.null(ret$MIXMIXbsIDX))
        dimnames(ret$MIXMIXbsIDX) <- list(1:nrow(ret$MIXMIXbsIDX), 1:ncol(ret$MIXMIXbsIDX))
    if(!is.null(ret$MIXCLAbsIDX))
        dimnames(ret$MIXCLAbsIDX) <- list(1:nrow(ret$MIXCLAbsIDX), 1:ncol(ret$MIXCLAbsIDX))
    if(!is.null(ret$CLACLAbsIDX))
        dimnames(ret$CLACLAbsIDX) <- list(1:nrow(ret$CLACLAbsIDX), 1:ncol(ret$CLACLAbsIDX))

    class(ret) <- "tclustICsol"
    ret
}

##  pll - matrix of size length(kk) x length(cc) containing penalized log likelihood
##      (MIXMIX, MIXCLA or CLACLA returned by tclustIC()).
##  ARI - matrix of size length(kk) x length(cc) containing Adjusted Rand
##      indexes between two consecutive values of c for fixed k
##  IDX - matrix of lists of size length(kk) x length(cc) containing the 
##      classification vector for each k (rows) and c (columns)
##  kk  - vector containing values of k which have to be considered
##  cc  - vector containing values of c which have to be considered
findBestSolutions <- function(pll, ARI, IDX, kk, cc, nsol, thresholdRI, trace=FALSE) {

    ARI <- t(ARI)
    Xcmod <- t(pll)
    lcc <- length(cc)
    seqcc <- 1:lcc
    seqkk <- 1:length(kk)

    ## Initialize Bestsols as a matrix of lists
    Bestsols <- matrix(vector("list", nsol * 5), nrow = nsol, ncol = 5)
    Bestsols[[1, 5]] <- 'true'

    endofloop <- FALSE
    lab <- "c="
    NumberOfExistingSolutions <- nsol

    for (z in 1:nsol) {
        valmin <- apply(Xcmod, 2, min)
        indmin <- apply(Xcmod, 2, which.min)
        minBICk <- which.min(valmin)
        minBICc <- indmin[minBICk]

        if(min(valmin) < Inf) {
            if(trace) {
                cat("Best solution number:", z, "\n")
                cat("k=", kk[minBICk], "\n")
                cat(lab, cc[minBICc], "\n")
            }

            Xcmod[indmin[minBICk], minBICk] <- Inf
            cwithbestsol <- rep(NA, length(cc))
            cwithbestsol[minBICc] <- 1
            Bestsols[[z, 1]] <- kk[minBICk]
            Bestsols[[z, 2]] <- cc[minBICc]
            
            if(trace) {
                cat("Checking adjacent values for extension\n")
            }

            XcmodWithoutBestk <- Xcmod
            XcmodWithoutBestk[, minBICk] <- Inf
            minICconstr <- min(XcmodWithoutBestk)
            cctoadd <- rep(0, length(cc))

            candcabove <- seqcc[seqcc > minBICc]
            if (length(candcabove) > 0) {
                for (r in candcabove) {
                    if (ARI[r, minBICk] > thresholdRI && Xcmod[r, minBICk] < minICconstr) {
                        Xcmod[r, minBICk] <- Inf
                        cctoadd[r] <- 1
                    } else {
                        break
                    }
                }
            }

            candcbelow <- seqcc[seqcc < minBICc]
            if (length(candcbelow) > 0) {
                for (r in rev(candcbelow)) {
                    if (ARI[r + 1, minBICk] > thresholdRI && Xcmod[r, minBICk] < minICconstr) {
                        Xcmod[r, minBICk] <- Inf
                        cctoadd[r] <- 1
                    } else {
                        break
                    }
                }
            }

            cwithbestsol[cctoadd == 1] <- 1
            
            ## Store the values of c associated to the best solution
            Bestsols[[z, 3]] <- cc[which(cwithbestsol == 1)]

            ##  The interval of values of c for which the solution is uniformly
            ##  better has been found, however, now we have to make sure that we
            ##  do not consider anymore the solutions for the same k which have a
            ##  value of R index adjacent to those which have already been found
            ##  greater than a certain threshold
            ##  Remark. If minBICk =1

            if (minBICk == 1) {
                Bestsols[[z, 4]] <- cc
                Xcmod[, minBICk] <- Inf
            } else {
                intc <- seqcc[which(cwithbestsol == 1)]
                outc <- cc[which(is.na(cwithbestsol))]
                cctoadd <- rep(0, length(cc))
                
                if(length(outc) > 0) {

                    ##  candcbelow - indexes of the set of candidate values for c which are smaller than
                    ##  those associated with best solutions found so far. 
                    ##  For example, if candc =(1 2) it means that we check whether for the same value of k
                    ##  solution with cc(2) and cc(1) have an R index greater than a certain threshold.
                    ##  If it is the case this means that these solutions do not have to be considered anymore
                
                    candcbelow <- seqcc[seqcc < min(intc)]
                    if (length(candcbelow) > 0) {
                        for (r in rev(candcbelow)) {
                            if (ARI[r + 1, minBICk] > thresholdRI) {
                                Xcmod[r, minBICk] <- Inf
                                cctoadd[r] <- 1
                            } else {
                                break
                            }
                        }
                    }
            
                    candcabove <- seqcc[seqcc > max(intc)]
                    if (length(candcabove) > 0) {
                        for (r in candcabove) {
                            if (ARI[r, minBICk] > thresholdRI) {
                                Xcmod[r, minBICk] <- Inf
                                cctoadd[r] <- 1
                            } else {
                                break
                            }
                        }
                    }
            
                    Bestsols[[z, 4]] <- cc[which(cctoadd == 1)]
                }
            }

            ##  get out of the loop because in this case there is a new candidate
            ##  solution (with the same k because the R index was smaller than a
            ##  certain threshold, or with a different k)
                    
            ##  Before getting out of the loop, check if the solution
            ##  which has just been found, has a Rand index greater than a
            ##  certain threshold with those which have already been found
            if(z > 1) {
                k_idx <- which(kk == Bestsols[[z, 1]])
                c_idx <- which(cc == Bestsols[[z, 2]])
                idxcurrentz <- IDX[[k_idx, c_idx]]
                
                for(j in 1:(z - 1)) {
                    prev_k_idx <- which(kk == Bestsols[[j, 1]])
                    prev_c_idx <- which(cc == Bestsols[[j, 2]])
                    idxpreviousz <- IDX[[prev_k_idx, prev_c_idx]]
                    
                    if(randIndex(idxpreviousz, idxcurrentz, 0)$AR > thresholdRI && Bestsols[[j, 5]] == 'true') {
                        Bestsols[[z, 5]] <- 'spurious'
                        break
                    } else {
                        Bestsols[[z, 5]] <- 'true'
                    }
                }
            }
        } else {    ## if Xcmod is full of inf than get out of the loop
            NumberOfExistingSolutions <- z - 1
            Bestsols <- Bestsols[1:NumberOfExistingSolutions, , drop = FALSE]
            endofloop <- TRUE
            break
        }
    }

    if(endofloop && trace) {
        cat("There are at most", NumberOfExistingSolutions, "different solutions\n")
        
        ##  break
    }

    ##  Find matrix of ARI for the z solutions which have been found
    ARIbest <- matrix(0, NumberOfExistingSolutions, NumberOfExistingSolutions)
    for (i in 1:NumberOfExistingSolutions) {
        for (j in 1:NumberOfExistingSolutions) {
            k_idx_i <- which(kk == Bestsols[[i, 1]])
            c_idx_i <- which(cc == Bestsols[[i, 2]])
            idxi <- IDX[[k_idx_i, c_idx_i]]
            
            k_idx_j <- which(kk == Bestsols[[j, 1]])
            c_idx_j <- which(cc == Bestsols[[j, 2]])
            idxj <- IDX[[k_idx_j, c_idx_j]]

            idxi[idxi < 0] <- 0
            idxj[idxj < 0] <- 0
            ARIbest[i, j] <- randIndex(idxi, idxj, 0)$AR
        }
    }

    return(list(Bestsols = Bestsols, ARIbest = ARIbest))
} 

#' @export
print.tclustICsol <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    cat("\nInformation criteria for TCLUST:", x$whichIC, "\n", paste0("Trimming = ", x$alpha))
    cat("\nNumber of mixture components (clusters):", x$kk)
    cat("\nvalues of the restriction factor:", x$cc, "\n")
    if(!is.null(x$MIXMIXbs)) {
        cat("\n\nBest solutions for MIXMIX (BIC):\n")
        ##  print.default(format(x$MIXMIXbs, digits = digits), print.gap = 2, quote = FALSE)
        print(x$MIXMIXbs)
    }
    if(!is.null(x$CLACLAbs)) {
        cat("\n\nBest solutions for CLACLA (CLA):\n")
        ##  print.default(format(x$CLACLAbs, digits = digits), print.gap = 2, quote = FALSE)
        print(x$CLACLAbs)
    }
    if(!is.null(x$MIXCLAbs)) {
        cat("\n\nBest solutions for MIXCLA (IC):\n")
        ##  print.default(format(x$MIXCLAbs, digits = digits), print.gap = 2, quote = FALSE)
        print(x$MIXCLAbs)
    }

     invisible(x)
}

#' @export
summary.tclustICsol <- function (object, ...) {
    ans <- list(tclustICobj=object)
    class(ans) <- "summary.tclustICsol"
    ans
}

#' @export
print.summary.tclustICsol <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",
	paste(deparse(x$tclustICobj$call), sep = "\n", collapse = "\n"), "\n", sep = "")

    cat("\nInformation criteria for TCLUST:", x$tclustICobj$whichIC, "\n", paste0("Trimming = ", x$tclustICobj$alpha))
    cat("\nNumber of mixture components (clusters):", x$tclustICobj$kk)
    cat("\nValues of the restriction factor:", x$tclustICobj$cc, "\n")

    invisible(x)
}

