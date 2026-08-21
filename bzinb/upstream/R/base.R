##########################################################################################
## 0. Common functions
##########################################################################################


#' @title Pairwise underlying correlation based on bivariate zero-inflated negative binomial (BZINB) model
#' 
#' @description For each pair of rows in the data, underlying corelation (\eqn{\rho}) is calculated based on
#'    bivariate zero-inflated negative binomial (BZINB) model.
#'    
#' @param data a matrix with nonnegative integers. rows represent the feature (or gene), and
#'    columns represent the sample. If not integers, rounded to the nearest integers.
#' @param nonzero.prop logical. If \code{TRUE}, proportion of nonzero for each of the pair is returned.
#' @param fullParam logical. If \code{TRUE}, estimates of all parameters are returned.
#' @param showFlag logical. If \code{TRUE}, for each pair, the estimates are printed out.
#' @param nsample positive integer. If provided, \code{nsample} random pairs will only be considered for 
#' correlation. A non-integer will be rounded to the nearest integer.
#' @param ... Other arguments passed on to \code{bzinb} function.
#' 
#' @return a table of pairwise underlying correlation (\eqn{\rho}) and related statistics.
#'  \itemize{
#'    \item \code{1} row number of the first vector of the pair
#'    \item \code{2} row number of the second vector of the pair
#'    \item \code{pair} row numbers of the pair
#'    \item \code{rho} underlying correlation estimate
#'    \item \code{se.rho} standard error of the underlying correlation estimate
#'    \item \code{nonzero.1, nonzero.2} non-zero proportion of the first and the second vector. 
#'          Returned if \code{nonzero.prop} is \code{TRUE}.
#'    \item \code{nonzero.min} pairwise minimum of non-zero proportions
#'          Returned if \code{nonzero.prop} is \code{TRUE}.
#'    \item \code{a0, a1, ..., p4} parameter estimates
#'    \item \code{se.a0, se.a1, ..., se.p4} standard error of the parameter estimates
#'    \item \code{logLik} log-likelihood of the maximum likelihood estimates
#'  }
#'  
#' @examples
#' # generating four random vectors
#' set.seed(7)
#' data1 <- rbzinb(n = 20, a0 = 0.5, a1 = 1, a2 = 1, 
#'                 b1 = 1, b2 = 1, p1 = 0.5, p2 = 0.2, 
#'                 p3 = 0.2, p4 = 0.1)
#' set.seed(14)
#' data2 <- rbzinb(n = 20, a0 = 0.5, a1 = 1, a2 = 1, 
#'                 b1 = 2, b2 = 2, p1 = 0.5, p2 = 0.2, 
#'                 p3 = 0.2, p4 = 0.1)
#' data3 <- t(cbind(data1, data2))
#' 
#' # calculating all pairwise underlying correlations
#' \dontrun{pairwise.bzinb(data3, showFlag = TRUE)}
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references 
#'  Cho, H., Liu, C., Preisser, J., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#' @import Rcpp
#' @export
#' 
pairwise.bzinb <- function(data, nonzero.prop = TRUE, fullParam = FALSE, 
                         showFlag = FALSE, nsample = NULL, ...) {
  # data: a dataframe of which pairs are to be analyzed. rows = genes, cols = sample
  # nonzero.prop: include nonzero proportion in the result
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  } else if (!is.matrix(data)) {
    stop ("data is neither matrix nor data.frame.")
  }
  if(!all(is.finite(data))) {stop("All elements in data must be finite and non-NA.")}
  if(any(data < 0)) {stop("All elements in data must be non-negative")}
  if (!is.null(nsample)) { 
    if (!is.numeric(nsample)) stop("nsample must be numeric, if it is not NULL.")
    nsample = round(nsample, 0)
    if (nsample < 1) stop("nsample must be greater than or equal to 1.")
  }
  
  dim.p <- dim(data)[1]
  comb <- expand.grid(1:dim.p, 1:dim.p)
  rw <- which(comb[, 1] > comb[, 2])
  comb <- data.frame(comb[rw, c(2,1)])
  comb$pair <- apply(comb, 1, function(x) paste(x, collapse = "-"))
  names(comb) <- c(1, 2, "pair")
  
  if (!is.null(nsample)) {
    if (nsample >= dim(comb)[1]) {
      warning("nsample is greater than the number of all possible pairs. No sampling is done.")
    } else {
      # random sampling nsample out of total number of pairs
      samp = sample(seq(dim(comb)[1]), nsample)
      comb = comb[samp, ]
    }
  } 
  
  MLE <- apply(comb[,1:2], 1, function(s) {
    # if (s[1] <= 6 | s[2] <=23) {return(data.frame(matrix(NA,1,4)))} # debug #7,24 has problem
    x <- data[s[1], ]
    y <- data[s[2], ]
    result <- bzinb(xvec = x, yvec = y, ...)
    
    if (showFlag) message("pair ", s[1], "-", s[2], ": ", "(rho, se.rho) = (", 
                          formatC(result$rho[1,1], digits = 5, format = "f", flag = "0"), ", ", 
                          formatC(result$rho[1,2], digits = 5, format = "f", flag = "0"), ")")
    result2 <- c(rho = result$rho["rho", "Estimate"], 
                 se.rho = result$rho["rho", "Std.err"],
                 if (fullParam) result$coefficients[, "Estimate"],  
                 if (fullParam) result$coefficients[, "Std.err"], 
                 if (fullParam) c(result$lik, result$iter))
    if (fullParam) {
      names(result2) <- c("rho", "se.rho", abp.names, 
                          paste0("se.", abp.names), "logLik", "num.iter")
    }
    return(result2)
  })
  comb <- cbind(comb, t(MLE))
  
  if (nonzero.prop == TRUE) {
    n <- dim(data)[2]
    comb$nonzero.1 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,1],] != 0) / n})
    comb$nonzero.2 <- sapply(1:dim(comb)[1], function(i) {sum(data[comb[i,2],] != 0) / n})
    comb$nonzero.min <- pmin(comb$nonzero.1, comb$nonzero.2)
  }
  
  return(comb)
}

# binary profile function: for BZIP.B
bin.profile <- function(xvec, yvec) {
  xvec[xvec != 0] = 1
  yvec[yvec != 0] = 1
  vec <- cbind(xvec,yvec)
  a <- rep(0,4)
  a[1] <- sum(apply(vec,1,prod)) # 1,1
  a[2] <- sum(xvec) - a[1]       # 1,0
  a[3] <- sum(yvec) - a[1]       # 0,1
  a[4] <- length(xvec) - sum(a)  # 0,0
  return(a)
}

.check.input <- function(xvec, yvec) {
  if(!all(is.finite(xvec)) | !all(is.finite(yvec))) {stop("xvec, yvec must be finite and non-NA.")}
  if(any(xvec < 0) | any(yvec < 0)) {stop("xvec, yvec must be non-negative")}
  if(!length(xvec) == length(yvec)){stop("The length of xvec and yvec must be the same.")}
}

.check.ab <- function(ab){
  # if (length(ab) != 5) {
  #   stop("Length(ab) must be 5 (a0, a1, a2, b1, b2).")
  # }
  if(any(!is.finite(ab))) {
    stop('a0, a1, a2, b1, b2 must be finite and non-NA.')
  }
  if(any(ab <= 0)) {
    stop('a0, a1, a2, b1, b2 must be greater than 0.')
  }
}
.check.m <- function(m){
  # if (length(m) != 3) {
  #   stop("Length(m) must be 3 (m0, m1, m2).")
  # }
  if(any(!is.finite(m))) {
    stop('m0, m1, m2 must be finite and non-NA.')
  }
  if(any(m <= 0)) {
    stop('m0, m1, m2 must be greater than 0.')
  }
}
.check.p <- function(p){
  # if (length(p) != 4) {
  #   stop("Length(p) must be 4 (p1, p2, p3, p4).")
  # }
  if(any(!is.finite(p))) {
    stop('p1, p2, p3, p4 must be finite and non-NA.')
  }
  if (any(p < 0) | any(p > 1)){
    stop('p1, p2, p3, p4 must be in [0, 1] inclusively.')
  }
  if(!isTRUE(all.equal(1, sum(p), tolerance = .Machine$double.eps^.5))) {
    stop(paste('sum of p1 to p4 must be 1.'))
  }
}
