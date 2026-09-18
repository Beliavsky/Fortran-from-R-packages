##' Calculates Rand type Indices to compare two partitions
##' @param c1 labels of the first partition or contingency table. A numeric vector
##' or factor containining the class labels of the first partition or a 
##' 2-dimensional numeric matrix which contains the cross-tabulation of 
##' cluster assignments.
##' @param c2 labels of the second partition. A numeric vector or a factor 
##' containining the class labels of the second partition. The length of 
##' the vector \code{c2} must be equal to the length of the vector \code{c1}. 
##' The second parameter is required only if c1 is not a 2-dimensional numeric matrix.
##' @param noisecluster label or number associated to the 'noise class' or 'noise level'.
##' Number or character label which denotes the points which do not belong to any
##' cluster. These points are not takern into account for the computation of the 
##' Rand type indexes. The default is to consider all points.
##'
##' @return A list with Rand type indexes:
##' \itemize{
##'     \item AR Adjusted Rand index. A number between -1 and 1. The adjusted 
##'         Rand index is the corrected-for-chance version of the Rand index. 
##'     \item RI Rand index (unadjusted). A number between 0 and 1. Rand index 
##'         computes the fraction of pairs of objects for which both 
##'         classification methods agree. RI ranges from 0 (no pair classified 
##'         in the same way under both clusterings) to 1 (identical clusterings).
##'     \item MI Mirkin's index. A number between 0 and 1. Mirkin's index computes 
##'         the percentage of pairs of objects for which both classification 
##'         methods disagree. \code{MI=1-RI}.
##'     \item HI Hubert index. A number between -1 and 1. HI index is equal to the 
##'         fraction of pairs of objects for which both classification methods 
##'         agree minus the fraction of pairs of objects for which both classification 
##'         methods disagree. \code{HI= RI-MI}. 
##' }
##' @examples
##' ##  1. randindex with the contingency table as input.
##' T <- matrix(c(1, 1, 0, 1, 2, 1, 0, 0, 4), nrow=3)
##' (ARI <- randIndex(T))
##' 
##' ##  2. randindex with the two vectors as input.
##' c <- matrix(c(1, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3), ncol=2, byrow=TRUE)
##' ## c1 = numeric vector containing the labels of the first partition
##' c1 <- c[,1]
##' ## c2 = numeric vector containing the labels of the second partition
##' c2 <- c[,2]
##' 
##' (ARI <- randIndex(c1,c2))
##' 
##' ##  3. Compare ARI for iris data (true classification against tclust classification)
##' library(tclust)
##' c1 <- iris$Species  # first partition c1 is the true partition
##' out <- tclust(iris[, 1:4], k=3, alpha=0, restr.fact=100)
##' c2 <- out$cluster   # second partition c2 is the output of tclust clustering procedure
##' 
##' randIndex(c1,c2)
##' 
##' ##  4. Compare ARI for iris data (exclude unassigned units from tclust).
##' 
##' c1 <- iris$Species      # first partition c1 is the true partition
##' out <- tclust(iris[,1:4], k=3, alpha=0.1, restr.fact=100)
##' c2 <- out$cluster       #  second partition c2 is the output of tclust clustering procedure
##' 
##' ## Units inside c2 which contain number 0 are referred to trimmed observations
##' noisecluster <- 0
##' randIndex(c1, c2, noisecluster=0)

randIndex <- function(c1, c2=NULL, noisecluster=NULL) {
  if(is.null(c2)) {
    if(ncol(as.matrix(c1)) < 2) {
      stop("RandIndex requires a contingency table with at least two columns")
    }
    C <- c1
  } else {
    if(min(dim(as.matrix(c1))) > 1 || min(dim(as.matrix(c2))) > 1) {
      stop("RandIndex requires two vector arguments")
    }
    
    if(!is.null(noisecluster)) {
      boo1 <- c1 == noisecluster
      boo2 <- c2 == noisecluster
      boo <- !(boo1 | boo2)
      c1 <- c1[boo]
      c2 <- c2[boo]
    }
    
    C <- table(c1, c2)  # Form contingency matrix
  }
  
  n <- sum(C)
  nis <- sum(rowSums(C)^2)  # Sum of squares of sums of rows
  njs <- sum(colSums(C)^2)  # Sum of squares of sums of columns
  
  totcomp <- choose(n, 2)   # Total number of pairs of entities
  t2 <- sum(C^2)            # Sum over rows and columns of nij^2
  t3 <- 0.5 * (nis + njs)
  
  # Expected value of the index (for adjustment)
  nc <- (n * (n^2 + 1) - (n + 1) * nis - (n + 1) * njs + 2 * (nis * njs) / n) / (2 * (n - 1))
  
  A <- totcomp + t2 - t3    # Total number of agreements
  D <- -t2 + t3             # Total number of disagreements
  
  if(totcomp == nc) {
    AR <- 0  # Avoid division by zero
  } else {
    AR <- (A - nc) / (totcomp - nc)  # Adjusted Rand Index
  }
  
  RI <- A / totcomp             # Rand Index
  MI <- D / totcomp             # Mirkin Index
  HI <- (A - D) / totcomp       # Hubert Index
  
  return(list(AR = AR, RI = RI, MI = MI, HI = HI))
}


