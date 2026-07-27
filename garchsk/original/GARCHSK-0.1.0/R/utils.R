#' This function calculates skewness of given data.
#' @param data vector or T by 1 matrix
#'
#' @return skewness of given data

skewness<-function(data){
  data_m<-data-base::mean(data)
  data_sd<-stats::sd(data)
  data_std<-data_m/data_sd
  return( base::mean(data_std^3) )
}

#' This function calculates kurtosis of given data.
#' @param data vector or T by 1 matrix
#'
#' @return kurtosis of given data

kurtosis<-function(data){
  data_m<-data-base::mean(data)
  data_sd<-stats::sd(data)
  data_std<-data_m/data_sd
  return( base::mean(data_std^4) )
}
