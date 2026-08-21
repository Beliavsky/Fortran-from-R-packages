#' @rdname mccut
#' @param lim a vector of quantiles for the uncertainty dimension. For \code{print.mccut}
#'   defaults to \code{c(0.025, 0.975)}; for \code{plot.mccut} defaults to
#'   \code{c(0.025, 0.25, 0.75, 0.975)}.
#' @param digits number of digits in the print.
#' @param ... additional arguments passed to the print or plot function.
#' @export
print.mccut <- function(x, lim=c(0.025,0.975), digits=3,...)
{
  summ <- function(x) c(mean=mean(x,na.rm=TRUE),quantile(x,probs=c(0.5,lim),na.rm=TRUE),Nas=sum(is.na(x)))[c(2,1,3:(length(lim)+3))]
  nbl <- length(x)

  summarize <- function(x){
   if(is.list(x)){
      for(i in 1:length(x)){
        x[[i]] <- summarize(x[[i]])
      }
   }
   else if(inherits(x,"mcnode")) return(x)
   else if(is.numeric(x)){
      dimm <- 1:length(dim(as.array(x)))
      dimm <- dimm[-2]
      x <- apply(x,dimm,summ)

   x <- drop(x)
   if(is.vector(x)) {
    x <- as.matrix(x)
    colnames(x) <- "NoVar"}
   }
   return(x)}

  x <- summarize(x)

  class(x) <- "listof"
  print(x,digits=digits,...)
  return(invisible(x))}
