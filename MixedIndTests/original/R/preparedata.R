#'@title Computes unique values, cdf and pdf
#'
#'@description This function computes the unique values, cdf and pdf for a series of data.
#'
#'@param x        Vector
#
#'
#'
#'@return \item{values}{Unique (sorted) values}
#'@return \item{m}{Number of unique values}
#'@return \item{Fn}{Empirical cdf of the unique values}
#'@return \item{fn}{Empirical pdf of the unique values}
#'
#'@references B.R. Nasri (2022). Tests of serial dependence for arbitrary distributions
#'@references C. Genest, J.G. Neslehova, B.N. Remillard and O. Murphy (2019). Testing for independence in arbitrary distributions.
#'

#'#'@examples
#'x = c(0,0,0,2,3,1,3,1,2,0)
#'out <- prepare_data(x)

#'@keywords internal
#'
#'@export
#'
#'

preparedata =function(x){

  values = sort(unique(x))
  n = length(x)
  m = length(values)


  out0 = .C("prepare_data",
            as.double(x),
            as.integer(n),
            as.double(values),
            as.integer(m),
            Fn = double(m),
            fn = double(m)

  )
  out0$values=values
  out0$m = m
  return(out0)
}






