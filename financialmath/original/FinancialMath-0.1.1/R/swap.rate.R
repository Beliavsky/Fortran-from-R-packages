swap.rate = function(rates, type="spot_rate"){

  all=list(rates,type)
  #Null
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #Length
  if(length(type)!=1) stop("type must be of length 1")
  #Numeric
  if(!is.vector(rates) | !is.numeric(rates))
    stop('variable "rates" must be numeric vector')
  #NA
  if(any(is.na(rates)) | is.na(type)) stop("Cannot input NA for any variables.")
  #Infinite
  if(any(rates==Inf)) stop("Rates cannot be infinite.")
  #Positive
  if(any(rates<0)) stop("All rates must be positive.")

  if(type != "spot_rate" & type != "zcb_price"){
    stop('type must be "spot_rate" or "zcb_price"')}

  if(type == "spot_rate"){
    for(i in 1:length(rates)){
      if(i == 1){runsum = 0}
      runsum = runsum + (1+rates[i])^(-i)
      if(i == length(rates)){out=(1-(1+rates[i])^(-i))/runsum}
    }}
  if(type == "zcb_price"){
    for(i in 1:length(rates)){
      if(i == 1){runsum = 0}
      runsum = runsum + rates[i]
      if(i == length(rates)){out=(1-rates[i])/runsum}
    }}
  return(out)
}
