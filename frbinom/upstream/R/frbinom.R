#' Fractional binomial distributions
#' @description
#' Generating random variables and computing density, cumulative distribution, and quantiles of the fractional binomial distribution with the parameters `size`, `prob`,  `h`, `c`.
#'
#' @param prob A number specifying the probability of success in each trial.
#' @param h A number specifying the strength of the dependency among trials; it determines  the skewness of the distribution.
#' @param c A number specifying the overdispersion of distributions.
#' @param size A number specifying the total number of trials.
#' @param n A number of random variables to be simulated.
#' @param p A numeric vector specifying probabilities at which quantiles of the fractional binomial distribution are computed.
#' @param x A numeric vector specifying values of the fractional binomial random variable at which the pmf or cdf is computed.
#' @param start logical; if TRUE, the starting point is changed after the first success in the generalized Bernoulli process. The default is FALSE.
#'
#' @returns A numeric vector of random variables (rfrbinom) or pmf (dfrbinom) or cdf (pfrbinom) or quantile (qfrbinom) of the fractional binomial distribution.
#' @details
#' The regular binomial random variable counts the number of successes in i.i.d. Bernoulli trials.
#' In the fractional binomial distribution,  the number of successes is counted among  dependent Bernoulli trials that are from the generalized Bernoulli process (GBP).
#' In GBP, each trial has the constant probability of success `prob`, and the strength of the dependency among the trials is determined by the parameters, `h` and `c`.  The parameters `c` and `h` are related to the overdispersion and skewness of the fractional binomial distribution.
#' The parameters `prob`=\eqn{p}, `h`=\eqn{H}, `c`=\eqn{c} should satisfy
#' \eqn{p,H \in (0,1)} and
#' \deqn{0\leq c<\min\{1-p, \frac{1}{2} ( -2p+2^{2H-2} +\sqrt{4p-p2^{2H}+2^{4H-4}})\}.}
#' With the number of trials (`size`) =\eqn{n}, the mean of the fractional binomial random variable is \eqn{E(X)=np},
#' and the variance is asymptotically proportional to \eqn{n^{2H}}, if \eqn{H\in (0.5, 1)};
#' \eqn{n\ln n}, if \eqn{H=0.5}; and \eqn{n}, if \eqn{H\in (0, .5)}.
#' If `c`=0, it becomes the regular binomial distribution.
#'
#' @references
#' Lee, J. (2021). Generalized Bernoulli process with long-range dependence and fractional binomial distribution. Dependence Modeling, 9(1), 1-12.
#'
#' @examples
#' # 10 random variables of a fractional binomial distribution.
#' rfrbinom(n=10, size=50, prob=.6, h=.8, c=.2)
#' # The probability that the fractional binomial random variable equals 22.
#' dfrbinom(x=22, size=50, prob=.6, h=.8, c=.2)
#' # The probability that the fractional binomial random variable is less than or equal to 22.
#' pfrbinom(x=22, size=50, prob=.6, h=.8, c=.2)
#'  # The 80th percentile of the fractional binomial distribution.
#' qfrbinom(p=.8, size=50, prob=.6, h=.8, c=.2)
#'
#' @export dfrbinom
#' @export pfrbinom
#' @export qfrbinom
#' @export rfrbinom
#'
#' @describeIn frbinom The pmf of fractional binomial distribution.
#' @export
dfrbinom<-function(x, size, prob, h, c, start = FALSE) {
  if( size==1){dbinom(x, size, prob)} else{
    if (1 < prob || prob < 0) {
      stop("Invalid value for probability parameter.")
    }

    if (1 < h || h < 0) {
      stop("Invalid value for h parameter.")
    }

    if (min(.5*(-2*prob+2^(2*h-2)+(4*prob-prob*2^(2*h)+2^(4*h-4))^(1/2)),1-prob) < c
        || c < 0) {
      stop("Inavlid value for c parameter.")
    }
    max2<-size-1

    p.fun<-function(X){
      p<-X[1]
      H<-X[2]
      c<-X[3]
      P<-rep(0,size)
      P[1]<-p+c
      d<-0
      for(i in 2:size){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
      P[i]<-p+c*i^(2*H-2)-d }; P;   }

    p.fun_0<-function(X){
      p<-X[1]
      H<-X[2]
      c<-X[3]
      P<-rep(0,size)
      P[1]<-p
      d<-0
      for(i in 2:size){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
      P[i]<-p-d }; P;}

    theo.p<-theo.p_0<-NULL
    PPm<-matrix(0, ncol=size, nrow=size); PP1<-c()

    yes_in<-(x %in% seq(0, size,1))
    final<-rep(0,length(x))
    if(start==FALSE){ theo.p_00<-p.fun_0(c(prob, h, c))
    theo.p<-p.fun(c(prob, h, c))
    theo.p_0<-1-cumsum(theo.p); PPm[1,]<-theo.p_00
    for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  }
    PP1<-PPm[,size]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
    final[yes_in]<-c(1-sum(theo.p_00),PP1)[(x + 1)[yes_in]]
    final
    }
    else{
      theo.p<-p.fun(c(prob, h, c))
      theo.p_0<-1-cumsum(theo.p) ;PPm[1,]<-theo.p
      for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  }
      PP1<-PPm[,size]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)
      final[yes_in]<-c(1-sum(theo.p),PP1)[(x + 1)[yes_in]]
      final }
  }}



#' @describeIn frbinom The cdf of fractional binomial distribution.
#' @export


pfrbinom<-function(x, size, prob, h, c, start = FALSE) {
  if (1 < prob || prob < 0) {
    stop("Invalid value for probability parameter.")
  }

  if (1 < h || h < 0) {
    stop("Invalid value for h parameter.")
  }

  if (min(.5*(-2*prob+2^(2*h-2)+(4*prob-prob*2^(2*h)+2^(4*h-4))^(1/2)),1-prob) < c
      || c < 0) {
    stop("Inavlid value for c parameter.")
  }
  final<-rep(0, length(x))
  final[which(x>size)]<-1
  belong<-which(x>=0 & x<=size)
  x2<-as.integer(x+1)[belong]
  if(start==FALSE){ dist <- dfrbinom(0:size, size, prob, h, c)

  final[belong]<-cumsum(dist)[x2]
  final}
  else{ dist <- dfrbinom(0:size, size, prob, h, c, start=TRUE)

  final[belong]<-cumsum(dist)[x2]
  final}

}


#' @describeIn frbinom The quantiles  of fractional binomial distribution.


#' @export
qfrbinom<-function(p, size, prob, h, c, start = FALSE) {
  if (1 < prob || prob < 0) {
    stop("Invalid value for probability parameter.")
  }

  if (1 < h || h < 0) {
    stop("Invalid value for h parameter.")
  }

  if (min(.5*(-2*prob+2^(2*h-2)+(4*prob-prob*2^(2*h)+2^(4*h-4))^(1/2)),1-prob) < c
      || c < 0) {
    stop("Inavlid value for c parameter.")
  }

  if(start==FALSE){
    cum.dist <- pfrbinom(0:size, size, prob, h, c)
  }else{  cum.dist <- pfrbinom(0:size, size, prob, h, c, start=TRUE) }
  q <- c()

  for (cur in 1:length(p)) {
    q[cur] <- min(which(cum.dist>=p[cur]))-1
  }

  q

}

#' @describeIn frbinom The generation of fractional binomial random variables.

#' @export

rfrbinom <- function(n, size, prob, h, c, start = FALSE) {
  if (1 < prob || prob < 0) {
    stop("Invalid value for probability parameter.")
  }

  if (1 < h || h < 0) {
    stop("Invalid value for h parameter.")
  }

  if (min(.5*(-2*prob+2^(2*h-2)+(4*prob-prob*2^(2*h)+2^(4*h-4))^(1/2)),1-prob) < c
      || c < 0) {
    stop("Inavlid value for c parameter.")
  }

  r<-c()

  r<-runif(n,0,1)

  r1<-c()
  if(start==FALSE){
    cum.dist <- pfrbinom(0:size, size, prob, h, c)}
  else{cum.dist <- pfrbinom(0:size, size, prob, h, c, start=TRUE) }

  for(i in 1:n) {
    r1[i]<-min(which(r[i] <= cum.dist))-1
  }

  r1
}


