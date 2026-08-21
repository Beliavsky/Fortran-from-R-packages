#' Computes Correlation Between Inputs and Output in the Variability Dimension (Tornado)
#'
#' Provides statistics for a tornado chart. Evaluates correlations between output and inputs
#' of a \samp{mc} object in the variability dimension.
#'
#' The tornado function computes Spearman's rho statistic to estimate a rank-based measure
#' of association between one set of random variables of a \samp{mc} object (the output)
#' and the others (the inputs).
#'
#' \samp{tornado} may be applied on a \samp{mccut} object if a \samp{tornado} function was
#' used in the third block of the \code{\link{evalmccut}} call.
#'
#' Type rules for the output node:
#' \itemize{
#'   \item \samp{"V" mcnode}: correlations are only provided for other \samp{"V" mcnode}s.
#'   \item \samp{"VU" mcnode}: correlations are provided for other \samp{"VU"} and \samp{"V"} nodes.
#' }
#'
#' @param mc a \code{\link{mc}} or \code{\link{mccut}} object.
#' @param x a \samp{tornado} object.
#' @param output the rank or name of the output to be considered. By default: the last element.
#' @param use an optional character string for computing covariances in the presence of
#'   missing values: \samp{"all.obs"}, \samp{"complete.obs"} or \samp{"pairwise.complete.obs"}.
#' @param method a character string for the correlation coefficient: \samp{"spearman"}
#'   (default), \samp{"kendall"} or \samp{"pearson"}.
#' @param lim a vector of quantiles used to compute the credible interval.
#' @param ... further arguments to be passed to the final print function.
#' @return an invisible object of class \samp{tornado}, containing: \samp{value} (correlation
#'   coefficients), \samp{output} (name of output), \samp{method}, \samp{use}.
#' @seealso \code{\link{cor}}, \code{\link{plot.tornado}} to draw the results.
#' @keywords univar
#' @examples
#' data(total)
#' tornado(total, 2, "complete.obs", "spearman", c(0.025, 0.975))
#' (y <- tornado(total, 10, "complete.obs", "spearman", c(0.025, 0.975)))
#' plot(y)
#' @export
tornado <- function(mc,output = length(mc),use = "all.obs",	method=c("spearman","kendall","pearson"),lim=c(0.025,0.975))
{
  if(inherits(mc,"mccut")){
    summ <- function(x) {
      y <-c(mean=mean(x,na.rm=TRUE),quantile(x,probs=c(0.5,lim),na.rm=TRUE))
      y[1:2] <- y[2:1]
      return(y)}
    mc <- mc[[which(sapply(mc,inherits,what="tornado.mccut"))[1]]]
    if(length(mc)==0) stop("tornado was not evaluated in evalmccut : impossible to evaluate")
    mc$value <- lapply(mc$value,function(x) drop(apply(x,c(1,3),summ)))
    class(mc) <- "tornado"
    return(mc)
  }

	method <- match.arg(method)
 	na.method <- pmatch(use, lesmet <- c("all.obs", "complete.obs", "pairwise.complete.obs"))

  if(!inherits(mc,"mc")) stop("tornado is not implemented for this object")
  if(is.numeric(output)) output <- names(mc)[output]
	if(!(output %in% names(mc)))	stop("Output is not a valid value")

  typen <- sapply(mc,attr,which="type")
  outm <-  lapply(mc,attr,which="outm")
  typeout <- typen[output]
  nameout <- names(mc[[output]])
	if(typeout!="V" && typeout!="VU")	stop("Output is not a 'V' or a 'VU' node")

  #Select nodes according to the type of the output
  quel <- outm!="none" & (typen == "V" | (typeout=="VU" & typen=="VU"))
  if(sum(quel) < 2) stop("No valid pairs of mcnode")      # Verify latter if one is the output

  # Select the nodes
  mc <- mc[quel]
  outm <- outm[quel]
  typen <- typen[quel]
  nom <- names(mc)
  dimm <- sapply(mc,dim)
  nvariates <- dimm[3,]
  nom <- names(mc)
  nomi <- nom[nom != output]



  # Which row are complete
  quelk <- !apply(sapply(mc,function(x) apply(x,1,function(x) any(is.na(x)))),1,any)

  if(!all(quelk)){
    if(na.method==1) stop("NA values. Change 'use' option")
    if(na.method==2){
      mc <- lapply(mc,"[",quelk,,,drop=FALSE)    # ! perte structure
      use <-  "all.obs"}
  }

  lesnom <- function(outm,nvariates,nom){
    if(outm[1] == "each"){
      if(nvariates==1) nomsortie <- nom
      else nomsortie <- paste(nom,1:nvariates,sep=".")}
    else nomsortie <- paste(nom,": ",outm," of variates",sep="")
    return(nom=nomsortie)}

  nomsortie <- mapply(lesnom,outm,nvariates,nom,SIMPLIFY=FALSE)
  rnbinput <- length(nomsortie)

  gerout <- function(node,outm,nvariates){
    if(outm[1]=="each"){
      if(nvariates==1) return(node)
      return(lapply(1:nvariates,function(x) node[,,x,drop=FALSE]))
      }
    res <- mapply(function(outm){
                  func <- get(outm,mode="function")
                  node <- apply(node,c(1,2),func)
                  dim(node) <- c(dim(node),1)
                  node},outm,SIMPLIFY=FALSE)
    if(length(outm)==1) return(res[[1]])
    return(res)
  }

  mc <- mapply(gerout,mc,outm,nvariates,SIMPLIFY=FALSE)

  # Deal with outm for the output
  nomout <- nomsortie[[output]]
  out <- mc[[output]]
  nco <-  dimm[2,output]

  nomin <- unlist(nomsortie[nomi])
  mc <- mc[nomi]


  # Build correlations
  calcorr <- function(inp,out){
    nunci <- dim(inp)[2]
   sapply(1:nco,function(y) cor(out[,y,],inp[,ifelse(nunci==1,1,y),],method=method,use=use))
  }

  if(!is.list(out)) out <- list(out)
  res <- rapply(out,function(x) matrix(rapply(mc,calcorr,how="unlist",out=x),nrow=nco,dimnames=list(NULL,nomin)),how="replace")

  if(typeout=="VU") {
    res <- lapply(res,function(x) {
                  tmp <- rbind(mean=colMeans(x,na.rm=TRUE),
                                  apply(x,2,quantile,probs=c(0.5,lim),na.rm=TRUE))
                  tmp[1:2,] <- tmp[2:1,]
                  rownames(tmp)[1:2] <- c("median","mean")
                  return(tmp)})
    }
  else   res <- lapply(res,"rownames<-",value="corr")

  names(res) <- nomout


  tc <- list(value = res, output = output, method = method, use = use)
	class(tc) <- "tornado"
  return(tc)
	}

#' @rdname tornado
#' @export
print.tornado <- function(x, ...)
{	tmethod <- c("Spearman's rho statistic","Kendall's tau statistic","Pearson correlation")
	tmethod <- tmethod[x$method==c("spearman","kendall","pearson")]
  cat(tmethod,"\n")
	cat("Output: ",x$output,"\n")
	print(x$value,...)
 }
