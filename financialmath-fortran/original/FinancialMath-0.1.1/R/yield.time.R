yield.time=function(cf,bal){

  #NULL
  if(any(lapply(list(cf,bal),is.null)==T)) stop("Cannot input any variables as NULL.")
  #Numeric
  if(!is.numeric(cf) | !is.vector(cf)) stop("cf must be a numeric vector")
  if(!is.numeric(bal) | !is.vector(bal)) stop("cf must be a numeric vector")
  #NA
  if(any(is.na(c(cf,bal)))) stop("Cannot input any variables as NA.")
  #Infinite
  if(any(c(cf,bal)==Inf)) stop("Cannot input any variables as infinite.")
  #Length
  if(length(cf)==length(bal)) cf=cf[-length(cf)]
  if(length(cf) != (length(bal)-1)) stop("length of cf must be one less than length of bal")

  runprod=1
  for(i in 1:length(cf)){
    runprod=runprod*bal[1+i]/(bal[i]+cf[i])
  }
  yield=runprod-1
  return(yield)
}
