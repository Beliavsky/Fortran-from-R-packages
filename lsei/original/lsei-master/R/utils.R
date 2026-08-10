##'Index-finding in a Sorted Vector
##'
##'
##'For each of given values, \code{indx} finds the index of the value in a
##'vector sorted in ascending order that the given value is barely greater than
##'or equal to.
##'
##'For each x[i], the function returns integer j such that \deqn{v_j \le x_i <
##'v_{j+1}}{v[j] <= x[i] < v[j+1]} where \eqn{v_0 = - \infty \mathrm{ and }
##'v_{n+1} = \infty}{v[0] = -Inf and v[n+1] = Inf}.
##'
##'@param x vector of numeric values, the indices of which are to be found.
##'@param v vector of numeric values sorted in ascending order.
##'@return Returns a vector of integers, that are indices of x-values in vector
##'v.
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@keywords array algebra
##'@examples
##'
##'indx(0:6,c(1:5,5))
##'indx(sort(rnorm(5)),-2:2)
##'
##'@export indx

indx = function(x, v) {
  m = length(x)
  n = length(v)
  x = pmax(v[1] - 1e300, pmin(v[n] + 1e300, x))  # x may contain -Inf or Inf
  storage.mode(x) = storage.mode(v) = "double"
  ind = integer(m)
  .C("indx", x, m, v, n, ind=ind, PACKAGE="lsei")["ind"]$ind
}

## # Much slower than indx().
## 
## indx2 = function(x, v) {
##   ox = order(x)
##   vx = c(v,x)
##   o = order(vx)
##   indo = c(rep(1, length(v)), rep(0, length(x)))[o]
##   indx = double(length(x))
##   indx[ox] = cumsum(indo)[indo == 0] 
##   indx
## }



##'Row or Column Maximum Values of a Matrix
##'
##'
##'Finds either row or column maximum values of a matrix.
##'
##'Matrix \code{x} may contain \code{Inf} or \code{-Inf}, but not \code{NA} or
##'\code{NaN}.
##'
##'@param x numeric matrix.
##'@param dim \code{=1}, for row maximum values; \code{=2}, for column maximum
##'values.
##'@return Returns a numeric vector with row or column maximum values.
##'
##'The function is very much the same as using \code{apply(x, 1, max)} or
##'\code{apply(x, 2, max)}, but faster.
##'@author Yong Wang <yongwang@@auckland.ac.nz>
##'@keywords array algebra
##'@examples
##'
##'x = cbind(c(1:4,Inf), 5:1)
##'matMaxs(x)
##'matMaxs(x, 2)
##'
##'@export matMaxs

matMaxs = function(x, dim=1) {
  if(length(x) == 0) return(NULL)
  x.mode = storage.mode(x)
  n = nrow(x)
  m = ncol(x)
  v = if(dim == 1) double(n) else double(m)
  storage.mode(dim) = storage.mode(n) = storage.mode(m) = "integer"
  x[x == Inf] = 1e308
  x[x == -Inf] = -1e308
  storage.mode(x) = "double"
  v = .C("matMaxs", x, n, m, v=v, dim, PACKAGE="lsei")["v"]$v
  v[v > 9e307] = Inf
  v[v < -9e307] = -Inf
  storage.mode(v) = x.mode
  v
}
