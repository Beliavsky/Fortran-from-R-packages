#' Build a Rank Correlation Using the Iman and Conover Method
#'
#' Builds a rank correlation structure between columns of a matrix or between \samp{mcnode}
#' objects using the Iman and Conover (1982) method.
#'
#' The function accepts for \samp{data}:
#' \itemize{
#'   \item some \samp{"V" mcnode} objects separated by a comma;
#'   \item some \samp{"U" mcnode} objects separated by a comma;
#'   \item some \samp{"VU" mcnode} objects, building the structure column by column;
#'   \item one \samp{"V" mcnode} as first element and some \samp{"VU" mcnode} objects.
#' }
#'
#' \samp{target} should be a scalar (two columns only) or a real symmetric positive-definite
#' square matrix. Only the upper triangular part is used (see \code{\link{chol}}).
#'
#' The final correlation structure should be checked, as it is not always possible to achieve
#' the target exactly. The order of values within each \samp{mcnode} (except the first) will
#' be changed by this function. The \samp{outrank} option may help to rebuild prior links.
#'
#' @param ... a matrix (each of its \samp{n} columns but the first one will be reordered) or
#'   \samp{n} \samp{mcnode} objects (each element but the first will be reordered). Arguments
#'   should be named.
#' @param target a scalar (only if \samp{n=2}) or a \samp{(n x n)} matrix of correlations.
#' @param outrank should the order be returned?
#' @param result should the correlation eventually obtained be printed?
#' @param seed the random seed for building the correlation. If \samp{NULL} the seed is unchanged.
#' @return If \samp{outrank=FALSE}: the matrix or a list of rearranged \samp{mcnode}s.
#'   If \samp{outrank=TRUE}: the order to be used to rearrange to build the desired correlation.
#' @references Iman, R. L., & Conover, W. J. (1982). A distribution-free approach to inducing
#'   rank correlation among input variables. \emph{Communication in Statistics - Simulation and
#'   Computation}, 11(3), 311-334.
#' @keywords multivariate
#' @examples
#' x1 <- rnorm(1000)
#' x2 <- rnorm(1000)
#' x3 <- rnorm(1000)
#' mat <- cbind(x1, x2, x3)
#' (corr <- matrix(c(1,0.5,0.2, 0.5,1,0.2, 0.2,0.2,1), ncol=3))
#' cor(mat, method="spearman")
#' matc <- cornode(mat, target=corr, result=TRUE)
#' all(matc[,1] == mat[,1])
#' @export
cornode <- function(...,target, outrank=FALSE, result=FALSE, seed=NULL)
{
  iman <- function(mat){  #Fonction de base de calcul
    R <- sapply(1:p,function(i) sample(phi,nvar))
    Rr <- apply(R,2,rank,ties.method="random")
    Ret <- Rr %*% P
    Rret <- apply(Ret,2,rank,ties.method="random")
    Xr <- apply(mat,2,order)
    rang <- sapply(1:p, function(i) (Xr[,i])[Rret[,i]])
    rang <- rang[order(rang[,1]),]
    if(outrank) return(rang)
    return(sapply(1:p, function(i) mat[rang[,i],i]))
  }

  data <- list(...)
  if(!is.null(seed)) set.seed(seed)
  if(!is.matrix(target)) target <- matrix(c(1,target,target,1),ncol=2)

  if(length(data)==1) { # A matrix
    data <- data[[1]]
    if(!is.matrix(data)) stop("data should be a matrix or a list of mcnodes")
    p <- ncol(data)
    nvar <- nrow(data)

    if(p < 2) stop("the matrix should have at least 2 columns")
    if(!is.matrix(target) || any(dim(target)!=c(p,p))) stop("target should be a matrix of correct dimension")

    phi <- qnorm((1:nvar)/(nvar+1))
    P <- chol(target)
    resu <- iman(data)

    colnames(resu) <- colnames(data)
    if(result)
    {cat("Spearman Rank Correlation Post Function\n")
      print(cor(resu,method="spearman"))}
    return(resu)
  }

  # mcnodes
  list.names <- function(...) {
    l <- as.list(substitute(list(...)))[-1]
    nm <- names(l)
    fixup <- if (is.null(nm)) seq(along = l) else nm == ""
    dep <- sapply(l[fixup], function(x) if (is.symbol(x)) as.character(x) else "")
    if (is.null(nm))
      return(dep)
    else {
      nm[fixup] <- dep
      return(nm)
    }
  }

  noms <- list.names(...)
  p <- length(data)
  if(p < 2) stop("the list should have at least 2 mcnodes")
  mcn <- sapply(data,inherits,"mcnode")
  if(!all(mcn))
    stop("the list should be a list of mcnode objects")

  if(!is.matrix(target) || any(dim(target)!=c(p,p)))
    stop("target should be a matrix of dimension l*l where l is the number of mcnodes in the list")

  tmcn <- sapply(data,attr,which="type")
  if(!(all(tmcn=="U") | all(tmcn %in% c("V","VU"))))
    stop("incorrect combination of mcnode: either 'U' or a combinaison of 'V' and 'U'")

  if(tmcn[1]=="U") {
    tmcd <- sapply(data,dim)
    for(i in 1:p) {
      data[[i]][] <- aperm(data[[i]],perm=c(2,1,3))
      dim(data[[i]]) <- tmcd[c(2,1,3),i] }                            # permute without change of structure
  }
  else
    if(any(tmcn=="VU") && any(tmcn=="V")){
      if(result){
        warning("impossible to provide the correlation result")
        result <- FALSE
      }
      if(sum(tmcn=="V")!=1 && tmcn[1]!="V")
        stop("Valid if only one 'V' node in the first position is combined with 'VU' nodes")
    }

  tmcd <- sapply(data,dim)
  nvar <- max(tmcd[1,])
  nunc <- max(tmcd[2,])
  nvariates <- max(tmcd[3,])

  if( !all(tmcd[1,] %in% c(1,nvar)) || !all(tmcd[2,] %in% c(1,nunc)))
    stop("incorrect dimension of mcnodes")
  if( any(tmcd[3,] != nvariates)    )
    stop("incorrect number of variates in mcnodes")

  res   <- data

  datai <- matrix(NA,ncol = p,nrow = nvar)
  phi <- qnorm((1:nvar)/(nvar+1))
  P <- chol(target)

  if(result) {corobs <- array(NA,dim=c(p*p,nunc,nvariates))}

  for(k in 1:nvariates){
    for(i in 1:nunc){
      ir <- ifelse(tmcd[2,1]==1, 1, i)
      datai[,1] <- data[[1]][,ir,k]
      for(j in 2:p) datai[,j] <- data[[j]][,i,k]
      sortim <- iman(datai)
      for(j in 2:p) res[[j]][,i,k] <-  sortim[,j]
      if(result) corobs[,i,k] <- cor(sortim,method="spearman",use="pairwise")
    }
  }

  if(result){
    if(nunc > 1)
    { corobs <- apply(corobs,c(1,3),function(x) c(mean=mean(x),quantile(x,c(0.5,0,1))))
    cat("summary of output Rank Correlation obtained accross the uncertainty dimension for each variates\n")}
    else
    {corobs <- aperm(corobs,perm=c(2,1,3))
    cat("output Rank Correlation per variates\n")}
    for(i in 1:nvariates){
      cat("variates:",i,"\n")
      print(corobs[,,i])}
  }

  if(outrank) res[[1]][] <- 1:nvar

  if(tmcn[1]=="U") {
    for(i in 1:p) {
      res[[i]][] <- aperm(res[[i]],perm=c(2,1,3))
      dim(res[[i]]) <- tmcd[c(2,1,3),i] }                            # repermute without change of structure
  }

  names(res) <- noms

  return(res)}
