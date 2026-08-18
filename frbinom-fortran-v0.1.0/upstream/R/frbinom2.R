#' Fractional binomial distribution II
#' @description
#' Generating random variables and computing density, cumulative distribution, and quantiles of the fractional binomial distribution II with the parameters `size`, `h`, `c`, `la`.
#'
#' @param h A number specifying the strength of the dependency among trials.
#' @param c A number specifying the dispersion of distributions.
#' @param size A number specifying the total number of trials.
#' @param la A number related to the probability of success in a trial; the default value is `c`/2.
#' @param n A number of random variables to be simulated.
#' @param p A numeric vector specifying probabilities at which quantiles of the fractional binomial distribution II are computed.
#' @param x A numeric vector specifying values of the fractional binomial random variable II at which the pmf or cdf is computed.
#' @param start logical; if TRUE, the starting point is changed after the first success in the generalized Bernoulli process II.  The default is FALSE.
#'
#' @returns A numeric vector of random variables (rfrbinom2) or pmf (dfrbinom2) or cdf (pfrbinom2) or quantile (qfrbinom2) of the fractional binomial distribution II.
#' @details
#' In the fractional binomial distribution II,  the number of successes is counted in the generalized Bernoulli process II (GBP II).
#' In GBP II, each trial has the constant probability of success `la`*`size`^\{2`h`-2\}, and the strength of the dependency among the trials is determined by the parameters, `h` and `c`.
#' The parameters `h`=\eqn{H}, `c`=\eqn{c}, `la`=\eqn{\lambda} should satisfy
#' \eqn{H \in (0.5,1)}, \eqn{0< c< 2^{2H-2}}, and  \eqn{0<\lambda<c}.
#' With the number of trials (`size`) =\eqn{n}, the mean of the fractional binomial random variable II is \eqn{E(X)=n\lambda^{2H-1}},
#' and the k-th moment is asymptotically proportional to \eqn{n^{(2H-1)k}} for \eqn{k\geq 2}.
#'
#' @references
#' Lee, J. (2023). Generalized Bernoulli process and fractional Poisson process. arXiv:2209.01516.
#'
#' @examples
#' # 10 random variables of a fractional binomial distribution II.
#' rfrbinom2(n=10, size=50, h=.8, c=.2)
#' # The probability that the fractional binomial random variable II equals 22.
#' dfrbinom2(x=22, size=50, h=.8, c=.2)
#' # The probability that the fractional binomial random variable II is less than or equal to 22.
#' pfrbinom2(x=22, size=50, h=.8, c=.2)
#'  # The 80th percentile of the fractional binomial distribution.
#' qfrbinom2(p=.8, size=50, h=.8, c=.2)
#' @importFrom stats dbinom runif
#' @export dfrbinom2
#' @export pfrbinom2
#' @export qfrbinom2
#' @export rfrbinom2
#'
#' @describeIn frbinom2 The pmf of fractional binomial distribution II.
#' @export
dfrbinom2<-function(x, size, h, c, la=c/2, start = FALSE) {

  if( size==1){ if(start==TRUE){dbinom(x, size, c)}else{
    pn<-la*size^(2*h-2);
    dbinom(x,size,pn) } } else{

      if (1 < h || h < 0.5) {
        stop("Invalid value for h parameter.") }

      if ( 2^(2*h-2) <= c
           || c <= 0) {
        stop("Inavlid value for c parameter.") }

      if ( c<=la || la <= 0) {
        stop("Invalid value for la parameter.")  }
      max2<-size-1

      theo.p<-theo.p_0<-NULL
      PPm<-matrix(0, ncol=size, nrow=size); PP1<-c()

      yes_in<-(x %in% seq(0, size,1))
      final<-rep(0,length(x))
      if(start==FALSE){
        p.fun2<-function(X){
          H<-X[1]
          c<-X[2]
          P<-rep(0,size)
          P[1]<-c
          d<-0
          for(i in 2:size){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(c*(i-j)^(2*H-2))*P[j]}
          P[i]<-c*i^(2*H-2)-d }; P;   }


        p.fun20<-function(X){
          H<-X[1]
          c<-X[2]
          la<-X[3]
          p<-la*(size^(2*H-2))
          P<-rep(0,size)
          P[1]<-p
          d<-0
          for(i in 2:size){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(c*(i-j)^(2*H-2))*P[j]}
          P[i]<-p-d }; P;   }

        theo.p<-p.fun2(c(h, c))
        theo.p_00<-p.fun20(c(h, c, la))
        theo.p_0<-1-cumsum(theo.p)
        PPm[1,]<-theo.p_00
        for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  }
        PP1<-PPm[,size]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
        final[yes_in]<-c(1-sum(theo.p_00),PP1)[(x + 1)[yes_in]]
        final
      }
      else{  p.fun2<-function(X){
        H<-X[1]
        c<-X[2]
        P<-rep(0,size)
        P[1]<-c
        d<-0
        for(i in 2:size){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(c*(i-j)^(2*H-2))*P[j]}
        P[i]<-c*i^(2*H-2)-d }; P;   }
      theo.p<-p.fun2(c(h, c))
      theo.p_0<-1-cumsum(theo.p)
      PPm[1,]<-theo.p
      for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  }
      PP1<-PPm[,size]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)

      final[yes_in]<-c(1-sum(theo.p),PP1)[(x + 1)[yes_in]]
      final }
    }}



#' @describeIn frbinom2 The cdf of fractional binomial distribution II.
#' @export


pfrbinom2<-function(x, size, h, c, la=c/2, start = FALSE) {

  if (1 < h || h < 0.5) {
    stop("Invalid value for h parameter.") }

  if ( 2^(2*h-2) <= c
       || c <= 0) {
    stop("Inavlid value for c parameter.") }

  if ( c<=la || la <= 0) {
    stop("Invalid value for la parameter.")  }
  final<-rep(0, length(x))
  final[which(x>size)]<-1
  belong<-which(x>=0 & x<=size)
  x2<-as.integer(x+1)[belong]
  if(start==FALSE){ dist <- dfrbinom2(0:size, size, h, c, la)

  final[belong]<-cumsum(dist)[x2]
  final}
  else{ dist <- dfrbinom2(0:size, size, h, c, start=TRUE)

  final[belong]<-cumsum(dist)[x2]
  final}

}


#' @describeIn frbinom2 The quantiles  of fractional binomial distribution II.


#' @export
qfrbinom2<-function(p, size, h, c, la=c/2 , start = FALSE) {
  if (1 < h || h < 0.5) {
    stop("Invalid value for h parameter.") }

  if ( 2^(2*h-2) <= c
       || c <= 0) {
    stop("Inavlid value for c parameter.") }

  if ( c<=la || la <= 0) {
    stop("Invalid value for la parameter.")  }

  if(start==FALSE){
    cum.dist <- pfrbinom2(0:size, size, h, c, la)
  }else{  cum.dist <- pfrbinom2(0:size, size,  h, c, start=TRUE) }
  q <- c()

  for (cur in 1:length(p)) {
    q[cur] <- min(which(cum.dist>=p[cur]))-1
  }

  q

}

#' @describeIn frbinom2 The generation of fractional binomial random variables II.

#' @export

rfrbinom2 <- function(n, size, h, c, la=c/2, start = FALSE) {

  if (1 < h || h < 0.5) {
    stop("Invalid value for h parameter.") }

  if ( 2^(2*h-2) <= c
       || c <= 0) {
    stop("Inavlid value for c parameter.") }

  if ( c<=la || la <= 0) {
    stop("Invalid value for la parameter.")  }

  r<-c()

  r<-runif(n,0,1)

  r1<-c()
  if(start==FALSE){
    cum.dist <- pfrbinom2(0:size, size,  h, c, la)}
  else{cum.dist <- pfrbinom2(0:size, size,  h, c, start=TRUE) }

  for(i in 1:n) {
    r1[i]<-min(which(r[i] <= cum.dist))-1
  }

  r1
}


