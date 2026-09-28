#' @title Matrix bending to positive-definite
#'
#' @description Bending a symmetric non-positive-definite matrix to positive-definite, using weighted or unweighted methods.
#'
#' @param inmat : The \code{matrix} to be bent.
#'
#' @param wtmat : The weight \code{matrix} for weighted bending. If no input is provided, the unweighted method (default) is used.
#'
#' @param reciprocal : If \code{TRUE}, reciprocal of the weighting factors are used. If no input is provided, default = \code{FALSE}.
#'
#' @param max.iter : Maximum number of iterations. If no input is provided, default = 10000.
#'
#' @param small.positive : Eigenvalues smaller than this value are replaced with this value. If no input is provided, default = 0.0001.
#'
#' @param method : \code{"hj"} (Jorjani et al., 2003) or \code{"lrs"} (Schaeffer, 2014), default = \code{"hj"}
#'
#' @return bent : The bent \code{matrix}.
#'
#' @return init.ev : Eigenvalues of the initial (\code{inmat}) matrix.
#'
#' @return final.ev : Eigenvalues of the \code{bent} matrix.
#'
#' @return min.dev : \code{min(bent - inmat)}.
#'
#' @return max.dev : \code{max(bent - inmat)}.
#'
#' @return loc.min.dev : Location (indices) of \code{min.dev} element.
#'
#' @return loc.max.dev : Location (indices) of \code{max.dev} element.
#'
#' @return ave.dev : Average deviation (\code{bent - inmat}) of the upper triangle elements (excluding diagonal elements for correlation matrices).
#'
#' @return AAD : Average absolute deviation of the upper triangle elements (excluding diagonal elements for correlation matrices) of \code{bent} and \code{inmat}.
#'
#' @return Cor : Correlation between the upper triangle elements (excluding diagonal elements for correlation matrices) of \code{bent} and \code{inmat}.
#'
#' @return RMSD : Root of mean squared deviation of the upper triangle elements (excluding diagonal elements for correlation matrices) of \code{bent} and \code{inmat}.
#'
#' @return w_gt_0 : Number of weight elements greater than 0, in the upper triangle of \code{wtmat} (for weighted bending).
#'
#' @return wAAD : Weighted \code{AAD} (for weighted bending).
#'
#' @return wCor : Weighted \code{Cor} (for weighted bending).
#'
#' @return wRMSD : Weighted \code{RMSD} (for weighted bending).
#'
#' @importFrom stats cor cov.wt
#'
#' @examples
#' # Test data
#' V = matrix(nrow=5, ncol=5, c( # matrix to be bent
#' 100,  95,  80,  40,  40,
#'  95, 100,  95,  80,  40,
#'  80,  95, 100,  95,  80,
#'  40,  80,  95, 100,  95,
#'  40,  40,  80,  95, 100))
#' W = matrix(nrow=5, ncol=5, c( # matrix of weights
#' 1000,  500,   20,   50,  200,
#'  500, 1000,  500,    5,   50,
#'   20,  500, 1000,   20,   20,
#'   50,    5,   20, 1000,  200,
#'  200,   50,   20,  200, 1000))
#'
#' # Example 1: Unweighted bending
#' bend(V)
#' ## The default method (Jojani et al. 2003) is used.
#'
#' # Example 2: Weighted bending using reciprocal of the weighting factors
#' bend(inmat=V, wtmat=W, reciprocal=TRUE)
#'
#' # Example 3: Bending with fixed elements
#' ## Assume we want to keep V[1:2, 1:2] constant.
#' W2 = W; W2[1:2, 1:2] = 0
#' bend(inmat=V, wtmat=W2, reciprocal=TRUE)
#'
#' # Example 4: Bending a correlation matrix
#' V2 = cov2cor(V)
#' bend(V2, W, reciprocal=TRUE)
#'
#' # Example 5: Bending using the method of Schaeffer (2014)
#' bend(inmat=V, method="lrs")
#'
#' # Example 6: Bending a correlation matrix using the method of Schaeffer (2014)
#' bend(V2, method="lrs")
#'
#' # Example 7: Bending the same correlation matrix using a weighted development of Schaeffer (2014)
#' bend(V2, W, reciprocal=TRUE, method="lrs")
#'
#' # Example 8: Bending a covariance matrix using a weighted development of Schaeffer (2014)
#' bend(V, W, reciprocal=TRUE, method="lrs")
#'
#' @export
bend = function(inmat, wtmat, reciprocal=FALSE, max.iter=10000, small.positive=0.0001, method="hj") {
  if(!is.matrix(inmat)) stop("inmat object is not a matrix.")
  if(!identical(rownames(inmat), colnames(inmat))) stop("rownames and colnames of the input matrix are not identical.")
  if(!isSymmetric(inmat)) stop("The matrix is asymmetric.")
  bent = inmat
  eigenval = eigen(bent, symmetric=TRUE)$values
  N = nrow(inmat)
  # Report what you've got.
  if(missing(wtmat)) {
    message("Unweighted bending")
    wtmat = matrix(1, nrow=N, ncol=N)
  } else {
    message("Weighted bending")
    if(!is.logical(reciprocal)) stop("reciprocal should be a logical variable.")
    message("reciprocal = ", reciprocal)
    if(!is.matrix(wtmat)) stop("wtmat object is not a matrix.")
    if(!identical(rownames(wtmat), colnames(wtmat))) stop("rownames and colnames of the weight matrix are not identical.")
    if(!isSymmetric(wtmat)) stop("The weight matrix is asymmetric.")
    if(ncol(inmat)!=ncol(wtmat)) stop("The matrices are incompatible in size.")
    if(length(wtmat[wtmat <0]) > 0) stop("Found negative elements in the weight matrix.")
  }
  message("max.iter = ", max.iter)
  if(method=="hj") message("small.positive = ", small.positive)
  message("method = ", method)
  if(min(eigenval) < 0) { # If non-PD
    # Get reciprocal of the weighting factors if needed
    if(reciprocal) {
      for(i in 1:nrow(wtmat))
      {
        for(j in 1:ncol(wtmat))
        {
          if(wtmat[i,j]!=0) wtmat[i,j] = 1/wtmat[i,j]
        }
      }
    }
    # Check if inmat is a correlation matrix
    correl = FALSE
    if(identical(unname(diag(inmat)), rep(1, N))) {
      message("Found a correlation matrix.")
      correl = TRUE
      diag(wtmat) = 0
    }
    init.ev = eigenval
    wtmat = wtmat/max(wtmat)
    if(!method %in% c("hj","lrs")) stop("Use method \"hj\" or \"lrs\".")
    eigenvec = eigen(bent, symmetric=TRUE)$vectors
    m = length(eigenval[eigenval < 0])
    i = 0
    while(i < max.iter & m > 0)
    {
      # Updating eigenvalues
      if(method=="hj") {
        if(small.positive <= 0 | small.positive > 0.1) stop("0 < small.positive <= 0.1 is not met.")
        eigenval[eigenval < small.positive] <- small.positive
      } else { # method="lrs"
        v = as.numeric(eigenval < 0)
        S = sum(v*eigenval)*2
        W = (S*S*100)+1
        small.positive = eigenval[N - m]
        k = N - m + 1
        for(j in k:N) eigenval[j] = small.positive * (S-eigenval[j])*(S-eigenval[j])/W
        # Make sure that UDtU is PD.
        UDtU = eigenvec %*% diag(eigenval) %*% t(eigenvec)
        eigenval = eigen(UDtU, symmetric=TRUE)$values
        eigenval[eigenval < 0] <- small.positive/10
      }
      UDtU = eigenvec %*% diag(eigenval) %*% t(eigenvec)
      bent = bent - (bent - UDtU) * wtmat
      eigenval = eigen(bent, symmetric=TRUE)$values
      eigenvec = eigen(bent, symmetric=TRUE)$vectors
      m = length(eigenval[eigenval < 0])
      i = i + 1
      if(round(i %% 1000)==0) message(paste('Iteration =', i))
    }
    if(length(eigenval[eigenval < 0]) > 0) {
      message(paste("WARNING: convergence was not met after", i, "iterations."))
    } else {
      message(paste("Convergence met after", i, "iterations."))
    }
  } else {
    stop("No action was required. The matrix is positive-definite.")
  }
  dev = bent - inmat
  loc.min.dev = which(dev==min(dev), arr.ind=TRUE)[1,]
  loc.max.dev = which(dev==max(dev), arr.ind=TRUE)[1,]
  if(correl) {
    dev = dev[upper.tri(dev, diag=FALSE)]
    w = wtmat[upper.tri(wtmat, diag=FALSE)]
    cor.dat = data.frame(V1=inmat[upper.tri(inmat, diag=FALSE)], V2=bent[upper.tri(bent, diag=FALSE)])
  } else {
    dev = dev[upper.tri(dev, diag=TRUE)]
    w = wtmat[upper.tri(wtmat, diag=TRUE)]
    cor.dat = data.frame(V1=inmat[upper.tri(inmat, diag=TRUE)], V2=bent[upper.tri(bent, diag=TRUE)])
  }
  outlist = list(bent = bent,
                 init.ev = init.ev,
                 final.ev = eigenval,
                 min.dev = min(dev),
                 max.dev = max(dev),
                 loc.min.dev = loc.min.dev,
                 loc.max.dev = loc.max.dev,
                 ave.dev = mean(dev),
                 AAD = mean(abs(dev)),
                 Cor = cor(cor.dat)[1,2],
                 RMSD = sqrt(mean(dev^2)))
  if(any(!w %in% 0:1)) {
    dev = dev[w > 0]
    cor.dat = cor.dat[w > 0,]
    w = w[w > 0]
    outlist[["w_gt_0"]] = length(w)
    outlist[["wAAD"]] = sum(abs(dev/w))/sum(1/w)
    outlist[["wCor"]] = cov.wt(cor.dat, wt=1/w, cor=TRUE)$cor[1,2]
    outlist[["wRMSD"]] = sqrt(sum((dev/w)^2)/sum(w^-2))
  }
  return(outlist)
}
