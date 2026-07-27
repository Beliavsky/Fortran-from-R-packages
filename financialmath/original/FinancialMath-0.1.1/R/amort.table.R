amort.table=function(Loan=NA,n=NA,pmt=NA,i,ic=1,pf=1,plot=FALSE) {

  all=list(Loan,n,pmt,i,ic,pf,plot)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(Loan,n,pmt,i,ic,pf)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("Loan, n, pmt, i, ic, and pf must be numeric.")
  #NA
  nalist=list(i,ic,pf,plot)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input i, ic, pf, and plot as NA.")
  #Logical
  stopifnot(is.logical(plot))

  num=c(Loan,n,pmt,i,ic,pf)
  NA.Neg=array(num)
  NA.Neg.Str=c("Loan","n","pmt","i","ic","pf")
  app=apply(NA.Neg,1,is.na)
  #Positive
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  test=c(n,pmt,Loan)
  if(length(test[which(is.na(test))])!=1) stop("Must have one NA (unknown) variable.")
  if(!is.na(n) & n !=round(n)) stop("n must be a positive integer.")

  eff.i=(1+i/ic)^(ic)-1
  int=(1+eff.i)^(1/pf)-1
  table=NA;table2=NA;balloon=NA;drop=NA;last.r=NA;last.t=NA;drop.t=NA;balloon.t=NA
  nom1=NA;nom2=NA;nk=0
  if(ic !=1) nom1=((1+eff.i)^(1/ic)-1)*ic
  if(pf != ic & pf !=1) nom2=((1+eff.i)^(1/pf)-1)*pf

  if(!is.na(Loan) & !is.na(pmt) & pmt<=Loan*(int)) stop("Too small of pmt.")
  if(!is.na(Loan) & !is.na(pmt) & pmt>=Loan*(1+int)) stop("Single payment will pay off loan.")

  if(pf!=1) {c=5;cn=c("Year","Payment","Interest Paid","Principal Paid","Balance")}
  else {c=4;cn=c("Payment","Interest Paid","Principal Paid","Balance")}

  if(is.na(Loan)) {Loan=pmt*(1-(1+int)^(-n))/int
  l=T
  } else l=F

  if(is.na(pmt) | l==T) {if(l==F) pmt=Loan/((1-(1+int)^-(n))/(int))
  table=matrix(c(rep(0,(c*(n)))),nrow=n,ncol=c)
  colnames(table)=cn
  rownames(table)=c(1:n)
  for (r in 1:n){
    if(pf!=1) table[r,"Year"]=r/pf
    table[r,"Payment"]=pmt
    if(r>=2)  table[r,"Interest Paid"]=table[r-1,"Balance"]*int
    if(r==1) table[r,"Interest Paid"]=Loan*int
    table[r,"Principal Paid"]=table[r,"Payment"]-table[r,"Interest Paid"]
    if(r==1)  table[r,"Balance"]=Loan-table[r,"Principal Paid"]
    if(r>=2)  table[r,"Balance"]=table[r-1,"Balance"]-table[r,"Principal Paid"]
  }
  if(plot==TRUE) {
    plot(1:n,table[,"Interest Paid"]/table[,"Payment"],
         xlab="Period",ylab="Percent",
         main="Percentage of Payment\nToward Interest",pch=20)
    lines(1:n,table[,"Interest Paid"]/table[,"Payment"],lty=2)
  }
  total.i=round(sum(table[,"Interest Paid"]),2)
  total.p=round(sum(table[,"Payment"]),2)
  }

  if(is.na(n)){nk=log(-(Loan/pmt*int-1))/log(1/(1+int))
  n=floor(nk);k=nk-n
  table2=matrix(c(rep(0,(c*(n+1)))),nrow=n+1,ncol=c)
  colnames(table2)=cn
  rownames(table2)=c(1:n,round(nk,2))
  for (r in 1:n){
    table2[r,"Payment"]=pmt
    if(pf!=1) table2[r,"Year"]=r/pf
    if(r>=2)  table2[r,"Interest Paid"]=table2[r-1,"Balance"]*int
    if(r==1) table2[r,"Interest Paid"]=Loan*int
    table2[r,"Principal Paid"]=table2[r,"Payment"]-table2[r,"Interest Paid"]
    if(r==1)  table2[r,"Balance"]=Loan-table2[r,"Principal Paid"]
    if(r>=2)  table2[r,"Balance"]=table2[r-1,"Balance"]-table2[r,"Principal Paid"]
  }

  table2[n+1,"Payment"]=pmt*((1+int)^k-1)/(int)
  table2[n+1,"Interest Paid"]=table2[n,"Balance"]*int*k
  table2[n+1,"Principal Paid"]=table2[n+1,"Payment"]-table2[n+1,"Interest Paid"]
  table2[n+1,"Balance"]=table2[n,"Balance"]-table2[n+1,"Principal Paid"]
  if(pf!=1) table2[n+1,"Year"]=nk/pf
  if(plot==TRUE) {
    plot(c(1:n,nk),table2[,"Interest Paid"]/table2[,"Payment"],
         xlab="Period",ylab="Percent",main="Percentage of Payment\nToward Interest")
    lines(c(1:n,nk),table2[,"Interest Paid"]/table2[,"Payment"],col="blue",lty=2)
  }
  balloon=Loan*(1+int)^n-pmt*((1+int)^n-1)/(int)+pmt
  drop=Loan*(1+int)^(1+n)-pmt*((1+int)^n-1)/(int)*(1+int)
  last.r=table2[n+1,"Payment"]
  last.t=nk
  drop.t=n+1
  balloon.t=n
  total.i=round(sum(table2[,"Interest Paid"]),2)
  total.p=round(sum(table2[,"Payment"]),2)
  if(round(k,2)==0) table2=table2[-(n+1),]
  }

  out1=if(any(is.na(table))) table2 else table
  out1=round(out1,2)
  out2=matrix(c(Loan,total.p,total.i,balloon,drop,last.r,eff.i,nom1,nom2),nrow=9)
  rownames(out2)=c("Loan","Total Paid","Total Interest",
                  "Balloon PMT","Drop PMT","Last Regular PMT",
                   "Eff Rate",paste("i^(",round(ic,2),")",sep=""),
                   paste("i^(",round(pf,2),")",sep=""))
  na=apply(out2, 1, function(x) all(is.na(x)))
  out2=as.matrix(out2[!na,])
  colnames(out2)=c("Details")
  out=list(Schedule=out1,Other=out2)
  return(out)
}

