#' @method plot multnet
#' @rdname plot.glmnet
#' @export
plot.multnet=function(x, xvar=c("lambda","norm","dev"),
                      label=FALSE, sign.lambda=-1, type.coef=c("coef","2norm"),...){
  xvar=match.arg(xvar)
  type.coef=match.arg(type.coef)
    beta=x$beta
    if(xvar=="norm"){
      cnorm1=function(beta){
        which=nonzeroCoef(beta)
        beta=as.matrix(beta[which,])
        apply(abs(beta),2,sum)
      }
      norm=apply(sapply(x$beta,cnorm1),1,sum)
    } else norm = NULL
    dfmat=x$dfmat
  if(type.coef=="coef"){
    ncl=nrow(dfmat)
    clnames=rownames(dfmat)
    for( i in seq(ncl)){
      plotCoef(beta[[i]],norm,x$lambda,dfmat[i,],x$dev.ratio,label=label,xvar=xvar,sign.lambda=sign.lambda,ylab=paste("Coefficients: Response",clnames[i]),...)
    }
  }
  else {
    dfseq=round(apply(dfmat,2,mean),1)
    plotCoef(coefnorm(beta,2),norm,x$lambda,dfseq,x$dev.ratio,label=label,xvar=xvar,sign.lambda=sign.lambda, ylab="Coefficient 2Norms",...)
  }

  }
