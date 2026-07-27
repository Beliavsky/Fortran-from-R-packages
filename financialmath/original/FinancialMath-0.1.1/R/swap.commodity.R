swap.commodity = function(prices, rates, type="spot_rate"){
  #Error Handling
  all=list(prices,rates,type)
  #NULL
  if(any(lapply(all,is.null)==T)) stop("Cannot input any variables as NULL.")
  #NA
  if(any(is.na(prices)) | any(is.na(rates))) stop("Cannot input any variables as NA.")
  #Numeric
  if(!is.vector(prices) | !is.vector(rates) |
     !is.numeric(prices) | !is.numeric(rates)){
    stop('variable "prices" and "rates" must be numeric vectors')}
  #Infinite
  if(any(prices==Inf) | any(rates==Inf)) stop("Variables cannot be infinite.")
  #Length
  if(length(type)!=1) stop("type must be of length 1")
  if(length(prices) != length(rates)){
    stop('length of price vector must be same length as rate vector')}

  if(type != "spot_rate" & type != "zcb_price"){
    stop('type must be one of "spot_rate" or "zcb_price"')}

  #Finding Values
  runsum = 0; runsum2 = 0
  if(type == "spot_rate"){
    for(i in 1:length(prices)){
      runsum = runsum + prices[i]/((1+rates[i])^i)
      runsum2 = runsum2 + 1/((1+rates[i])^i)}
    out = runsum/runsum2}
  if(type == "zcb_price"){
    for(i in 1:length(prices)){
      runsum = runsum + prices[i]*rates[i]
      runsum2 = runsum2 + rates[i]}
    out = runsum/runsum2}
  return(out)
}
