#' @title Calculates Credit Default Swap rates
#'
#' @description
#' Calculates CDS rates starting form default intensities.
#'
#' @param t premium timetable.
#' @param int deterministic default intensities vector.
#' @param r spot interest rates.
#' @param R constant premium payments, value that the buyer pays in each \eqn{t_i}.
#' @param RR recovery rate on the underline bond, default value is 40\%.
#' @param simplified logic argument. If FALSE calculates the CDS rates using the
#' semplified version of calculations, if TRUE use the complete version.
#'
#' @details
#' \itemize{
#' \item  Premium timetable is \eqn{t_i;  i=1,...,T}. The vector starts from
#' \eqn{t_1\le 1}, i.e. the first premium is payed at a year fraction in the possibility that
#'  the bond is not yet defaulted. Since premium are a postponed payment (unlike usual insurance
#'   contracts).
#' \item Intensities timetable have domains \eqn{\gamma_i;  i=t_1,...,T}.
#' \item  spot interest rates of bond have domain \eqn{r_i; i=t_1,...,T}. The function transforms
#'  spot rates in forward rates. If we specify that we want to calculate CDS rates with the
#'  simplified alghoritm, in each period, the amount of the constant premium payment
#'  is expressed by: \deqn{\pi^{pb}=\sum_{i=1}^Tp(0,i)S(0,i)\alpha_i}
#'  and the amount of protection, assuming a recovery rate \eqn{\delta}, is:
#'  \deqn{\pi^{ps}=(1-\delta)\sum_{i=1}^Tp(0,i)\hat{Q}(\tau=i)\alpha_i}
#'  If we want to calculate same quantities with the complete version, that evaluate premium in the continous,
#'  the value of the premium leg is calculated as:
#'  \deqn{\pi^{pb}(0,1)=-\int_{T_a}^{T_b}P(0,t)\cdot(t-T_{\beta(t)-1}) d_t Q
#'  (\tau\geq t)+\sum_{i=a+1}^bP(0,T_i)\cdot\alpha_i * Q(\tau\geq T_i)} and the protection leg
#'  as:
#'  \deqn{\pi_{a,b}^{ps}(1):=-\int_{t=T_a}^{T_b}P(0,t)d*Q(\tau\geq t)}
#'  In both versions the forward rates and intensities are supposed as costant stepwise functions
#'   with discontinuity in \eqn{t_i}
#'}
#'
#' @return
#' \code{cds} returns an object of class \code{data.frame} with columns, for esch date
#' \eqn{t_i} the value of survival probability, the premium and protection leg, CDS rate
#' and CDS price.
#'
#' @references
#' David Lando  (2004) Credit risk modeling.
#'
#' Damiano Brigo, Massimo Morini, Andrea Pallavicini (2013)
#' Counterparty Credit Risk, Collateral and Funding.
#' With Pricing Cases for All Asset Classes
#'
#'
#' @examples
#' cds(t = seq(0.5, 10, by = 0.5), int = seq(.01, 0.05, len = 20),
#' r = seq(0,0.02, len=20), R = 0.005, RR = 0.4, simplified = FALSE)
#'
#' @export
#'
cds <- function(t, int, r, R = 0.005, RR = 0.4, simplified = FALSE){

	if( !is.numeric(t) ) stop( "Time t must be a numeric vector" )
	if( !is.numeric(int) ) stop( "Intensity int must be a numeric vector" )
	if( !is.numeric(r) ) stop( "Interest spot rates r must be a numeric vector" )
	if( RR == 1 ) stop( "Recovery rate cannot be 100%" )
	if( length(t) != length(int)) stop( "t and int must have same length" )
	if( R == 0 ) stop ( "Do you pretend any protection paying 0?" )

	n <- length(t)
	dt <- diff( c(0, t) )
	fs <- exp(-r*t)
	intcum <- cumsum(int*dt)
	Q <- exp(-intcum)

	frw = 0;

	if( simplified == F ){

		frw <- diff(c(0,r*t))/dt

		w <- int /(int+frw)^2
		d <- int/(int+frw)
		k <- dt + ( int*( 1+dt*(int+frw) )/(int+frw)^2)

		B <- fs*Q*k;
 		A	 <- rep(0, n)
		A[1]	 <- w[1]*1*1
		A[2:n] <- w[2:n]*fs[1:n-1]*Q[1:n-1]
      	C<- -d*( diff( c(1, fs*Q) ) )

	} else {
		A <- rep(0,n)
		B <- fs*Q*dt
		C <- fs*( -diff( c(1, Q) ) )
	}

	premium <- cumsum(B-A)
	protection <- (1-RR)*cumsum(C)
	rate <- protection/premium
	premium <- premium*R
	price <- -premium + protection

	out<-data.frame(t, Q, premium, protection, rate, price)
	names(out) <- c("t","Survival","PremiumLeg","ProtectionLeg","Rate","Price")
      
	return(out)

}
