IRR=function(cf0,cf,times,plot=FALSE){

  all=list(cf0,cf,times,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(list(cf0,plot),length) != 1)==T) stop("cf0 and plot must be of length 1.")
  #Numeric
  if(!is.vector(cf) | !is.numeric(cf)) stop("cf must be a numeric vector.")
  if(!is.vector(times) | !is.numeric(times)) stop("times must be a numeric vector.")
  if(!is.numeric(cf0)) stop("cf0 must be numeric.")
  #NA
  if(any(is.na(cf)) | any(is.na(times)) | is.na(cf0) | is.na(plot)) stop("Cannot input NA for any variables.")
  #Infinite
  if(any(cf==Inf) | any(times==Inf) | cf0==Inf) stop("Cannot input infinite for any variables.")
  #Logical
  stopifnot(is.logical(plot))

  if(length(cf)==0 | length(times)==0 ) stop("Not enough cash flow information.")
  if(length(cf) != length(times)) stop("Number of cash flows not equal to number of period values.")

  cf0=round(abs(cf0),2)
  if(any(times<=0)) stop("Cannot have negative period.")
  if(any(times%%1!=0)) stop("Periods must be integers.")
  d=unique(times[which(duplicated(times)==T)])
  if(length(d)>0){ cfd=rep(0,length(d))
  for (r in 1:length(d)) cfd[r]=sum(cf[which(times==d[r])])}
  p=rep(0,max(times)+1)
  p[times+1]=cf
  p[1]=-cf0
  if(length(d)>0) {for(r in 1:length(d)) p[d[r]+1]=cfd[r]}
  irr=1/Re(polyroot(p))-1
  if(any(irr<=1 & irr>=0)) irr=irr[which(irr<=1 & irr>=0)]
  if(any(irr>0)) irr=irr[which(irr>0)]
  if(plot==T){
    plot.new()
    plot(0,0,type="n",axes=F,ann=F,xlim = c(0,max(times)),ylim=c(0, 10))
    axis(1,at=seq(0,max(times)),labels=seq(0,max(times)),line=-5)
    text(seq(1,max(times)),rep(5.5,max(times)),labels=p[-1],cex=.85)
    text(max(times)/2,10,labels="Time Diagram",cex=1.2)
    text(max(times)/2,8.5,bquote("IRR "== .(round(irr,4))),cex=.85)
    text(0,5.5,labels=-cf0,cex=.75)  }
  return(irr)
}
