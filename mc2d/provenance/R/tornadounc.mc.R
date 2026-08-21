#' Computes Correlation Between Inputs and Output in the Uncertainty Dimension (Tornado)
#'
#' Provides statistics for a tornado chart. Evaluates correlations between output and inputs
#' of a \samp{mc} object in the uncertainty dimension.
#'
#' The function computes Spearman's rho statistic between values or statistics calculated in
#' the variability dimension of inputs and outputs. The statistics are the mean, the standard
#' deviation and the quantiles specified by \samp{quant}.
#'
#' \samp{tornadounc.mccut} may be applied on a \code{\link{mccut}} object if a
#' \samp{summary.mc} function was used in the third block of the \code{\link{evalmccut}} call.
#'
#' @param mc a \samp{mc} object.
#' @param x a \samp{tornadounc} object.
#' @param output the rank or name of the output. Should be a \samp{"VU"} or \samp{"U"} node.
#'   By default: the last element.
#' @param quant the vector of quantiles used in the variability dimension.
#' @param use an optional character string for computing covariances in presence of missing
#'   values: \samp{"all.obs"}, \samp{"complete.obs"} or \samp{"pairwise.complete.obs"}.
#' @param method a character string for the correlation coefficient: \samp{"spearman"}
#'   (default), \samp{"kendall"} or \samp{"pearson"}.
#' @param ... further arguments to be passed to the final print function.
#' @return an invisible object of class \samp{tornadounc}.
#' @seealso \code{\link{cor}}, \code{\link{tornado}} for variability dimension,
#'   \code{\link{plot.tornadounc}} to draw the results.
#' @keywords univar
#' @name tornadounc
#' @examples
#' data(total)
#' tornadounc(total, 3)
#' (y <- tornadounc(total, 10, use="complete.obs"))
#' plot(y, 1, 1)
#' @export
tornadounc <- function(mc,...)
 UseMethod("tornadounc")

#' @rdname tornadounc
#' @export
tornadounc.default <- function(mc,...)
 tornadounc.mc(mc,...)

#' @rdname tornadounc
#' @export
tornadounc.mc <- function(mc,output = length(mc), quant=c(0.5,0.75,0.975),use = "all.obs",	method=c("spearman","kendall","pearson"),...)
{
	method <- match.arg(method)
 	na.method <- pmatch(use, lesmet <- c("all.obs", "complete.obs", "pairwise.complete.obs"))

  if(length(output) > 1) stop("Only one output permitted")
  if(is.numeric(output)) output <- names(mc)[output]
	if(!(output %in% names(mc)))	stop("Output is not a valid value")

  typen <- sapply(mc,attr,which="type")
  outm <-  lapply(mc,attr,which="outm")
  typeout <- typen[output]
	if(typeout!="U" && typeout!="VU")	stop("Output is not a 'U' or a 'VU' node")
  if(outm[output]=="none") stop("Output has a 'none' outm attribute")

  #Select nodes according to the type of the output
  quel <- outm!="none" & (typen == "U" | (typeout=="VU" & typen=="VU"))
  if(sum(quel) < 2) stop("No valid pairs of mcnode")

  # Select the nodes
  mc <- mc[quel]
  outm <- outm[quel]
  typen <- typen[quel]
  dimm <- sapply(mc,dim)
  nvariates <- dimm[3,]
  nom <- names(mc)
  nomi <- nom[nom != output]

   # Which col are complete
  quelk <- !apply(sapply(mc,function(x) apply(x,2,function(x) any(is.na(x)))),1,any)
  nco <- sum(quelk)              # remaining dimension of uncertainty

  if(!all(quelk)){
    if(na.method==1) stop("NA values. Change 'use' option")
    if(na.method==2){
      mc <- lapply(mc,"[",,quelk,,drop=FALSE)    # ! perte structure
      use <-  "all.obs"}
  }

  # Real name of input / output

  yaprob <- length(quant) > 0
  if(yaprob) {nom1 <- paste(quant*100,"%",sep="")
  	nom1[quant==0] <- "Min"
  	nom1[quant==1] <- "Max"} else nom1 <- NULL
  nom1 <- c("mean","sd",nom1)

  lesnom <- function(nom,outm,nvariates,typen){
    if(outm[1] == "each"){
      if(nvariates==1) nomsortie <- nom
      else nomsortie <- paste(nom,1:nvariates,sep=".")}
    else nomsortie <- paste(nom,": ",outm," of variates",sep="")

    if(typen == "VU")
      nomsortie <- lapply(nomsortie, function(x) paste(nom1,x))

    return(nom=nomsortie)}

  nomsortie <- mapply(lesnom,nom,outm,nvariates,typen,SIMPLIFY=FALSE)
  nomlistfin <- lesnom(output,outm[[output]],nvariates[[output]],"U")           # Name final list

  # Deal with multivariate nodes

  gerout <- function(node,name,outm,nvariates){
    if(outm[1]=="each")
      res <- lapply(1:nvariates,function(x) node[,,x,drop=FALSE])
    else
        res <- mapply(function(outm){
                  func <- get(outm,mode="function")
                  node <- apply(node,c(1,2),func)
                  dim(node) <- c(dim(node),1)
                  node},outm,SIMPLIFY=FALSE)

  return(res)
  }

  mc <- mapply(gerout,mc,nom,outm,nvariates,SIMPLIFY=FALSE)

  # Build Statistics on inputs and output

  stat <- function(node,outm,typen,nomsortie){
    if(typen=="U")
      tmp <- lapply(node,function(x) x[1,,])
    else {
      tmp <- lapply(node,function(x)
                              apply(x,c(2,3),function(x) c(mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE), # Compliance with the deprecation of sd(<matrix>)
                                              quantile(x,na.rm=TRUE,prob=quant))))
                                              }
      node <- lapply(tmp,matrix,nrow=nco,byrow=TRUE)

  return(node)
  }

  mc <- mapply(stat,mc,outm,typen,nomsortie,SIMPLIFY=FALSE)

  nomout <- nomsortie[[output]]
  out <- mc[[output]]

  nomin <- unlist(nomsortie[nomi])
  mc <- matrix(unlist(mc[nomi]),nrow=nco)

  lescorr <- mapply(function(x) as.matrix(cor(x,mc,method=method,use=use)),out,SIMPLIFY=FALSE)
  lescorr <- lapply(lescorr,"colnames<-",value=nomin)
  lescorr <- mapply("rownames<-",lescorr,value=nomout,SIMPLIFY=FALSE,USE.NAMES=TRUE)
  names(lescorr) <- nomlistfin

  tc <- list(value = lescorr, output = output, method = method, use = use)
	class(tc) <- "tornadounc"
  return(tc)
	}

#' @rdname tornadounc
#' @export
print.tornadounc <- function(x, ...)
{
tmethod <- c("Spearman's rho statistic","Kendall's tau statistic","Pearson correlation")
	tmethod <- tmethod[x$method==c("spearman","kendall","pearson")]
  cat("Tornado on uncertainty\n")
  cat(tmethod,"\n")
	cat("Output: ",x$output,"\n")
	print(x$value,...)
 }
