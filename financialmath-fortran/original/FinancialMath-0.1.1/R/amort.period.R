amort.period=function(Loan=NA,n=NA,pmt=NA,i,ic=1,pf=1,t=1){

  all=list(Loan,n,pmt,i,ic,pf,t)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(any(lapply(all,length) != 1)==T) stop("All inputs must be of length 1.")
  #Numeric
  num2=list(Loan,n,pmt,i,ic,pf,t)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("Loan, n, pmt, i, ic, pf, and t must be numeric.")
  #NA
  nalist=list(i,ic,pf,t)
  if(any(lapply(nalist,is.na)==T)) stop("Cannot input i, ic, pf, and t as NA.")
  #Logical

  NA.Neg=array(c(Loan,n,pmt,i,ic,pf,t))
  NA.Neg.Str=c("Loan","n","pmt","i","ic","pf","t")
  #Positive
  app=apply(NA.Neg,1,is.na)
  na.s=which(app==F & NA.Neg<=0)
  if(length(na.s)>0) {errs=paste(NA.Neg.Str[na.s],collapse=" & ")
  stop(cat("Error: '",errs, "' must be positive real number(s).\n"))}
  #Infinite
  na.s2=which(app==F & NA.Neg==Inf)
  if(length(na.s2)>0) {errs=paste(NA.Neg.Str[na.s2],collapse=" & ")
  stop(cat("Error: '",errs, "' cannot be infinite.\n"))}

  test=c(Loan,n,pmt)
  if(length(test[which(is.na(test))])!=1) stop("Must have one NA (unknown) variable.")

  if(!is.na(n) & n!=round(n)) stop("n must be a positive integer.")
  #if(t !=round(t) & n.as.years==F) stop("t must be a positive integer")
  #if(t*pf != round(t*pf) & n.as.years==T) stop("t*pf must be an integer")

  nom1=NA;nom2=NA

  eff.i=(1+i/ic)^(ic)-1
  int=(1+eff.i)^(1/pf)-1
  if(ic !=1) nom1=((1+eff.i)^(1/ic)-1)*ic
  if(pf !=1 & ic !=pf) nom2=((1+eff.i)^(1/pf)-1)*pf

  if(!is.na(Loan) & !is.na(pmt) & pmt<=Loan*(int)) stop("Too small of pmt.")
  if(!is.na(Loan) & !is.na(pmt) & pmt>=Loan*(1+int)) stop("Single payment will pay off loan.")

  if(is.na(n)){n=log(-(Loan/pmt*int-1))/log(1/(1+int))}

  if(t>n) stop("t cannot be greater than n.")

  if(is.na(Loan)) Loan=pmt*(1-(1+int)^(-n))/int
  if(is.na(pmt))pmt=Loan/((1-(1+int)^-(n))/(int))

  t.int=pmt*int*(1-(1+int)^(-(n-t+1)))/int
  t.princ=pmt*(1/(1+int))^(n-t+1)
  t.bal=pmt*(1-(1+int)^(-(n-t)))/int

  if(pf==1) {years=n;n=NA} else years=n/pf
  out=c(Loan,pmt,eff.i,nom1,nom2,n,years,t,t.int,t.princ,t.bal)
  m.out=matrix(out,nrow=length(out))
  rownames(m.out)=c("Loan","PMT","Eff Rate",paste("i^(",round(ic,2),")",sep=""),
                    paste("i^(",round(pf,2),")",sep=""),"Periods","Years",
                    paste("At Time ",round(t,2),":",sep=""),
                    "Int Paid","Princ Paid","Balance")
  na=apply(m.out, 1, function(x) all(is.na(x)))
  m.out=as.matrix(m.out[!na,])
  colnames(m.out)=c("Amortization")
  return(round(m.out,6))
}
