yield.dollar = function(cf, times, start, end, endtime){

  all=list(cf,times,start,end,endtime)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  one=list(start,end,endtime)
  if(any(lapply(one,length)!= 1)) stop("start, end, and endtime must be of length 1.")
  if(length(cf) != length(times)){
    stop('length of vectors "cf" and "times" must be equal')}
  #NA
  if(any(is.na(cf) | any(is.na(times)) | any(is.na(c(start,end,endtime)))))
    stop("Cannot input any variables as NA.")
  #Numeric
  if(!is.vector(cf) | !is.vector(times) | !is.numeric(cf) | !is.numeric(times))
    stop('variable "cf" and "times" must be numeric vector')
  num2=list(start,end,endtime)
  na.num2=num2[which(lapply(num2,is.na)==F)]
  if(any(lapply(na.num2,is.numeric)==F)) stop("All variables must be numeric.")
  #Infinite
  if(any(cf==Inf) | any(times==Inf) | any(c(start,end,endtime)==Inf))
    stop("Cannot input any variables as infinite.")
  #Positive
  if(any(times<0) | endtime <=0)
    stop('"times" vector and "endtime" must be positive')
  if(any(times>endtime))
    stop('time of comparison (endtime) must be larger than any number in vector "times"')
  if(length(cf) != length(times))
    stop("length of cashflow vector and times vector must be equal")

  I = end - start - sum(cf)
  runsum=0
  for(i in 1:length(cf)){
    runsum = runsum + cf[i]*(endtime-times[i])}
  yield = I/(start*endtime+runsum)
  return(yield)
}
