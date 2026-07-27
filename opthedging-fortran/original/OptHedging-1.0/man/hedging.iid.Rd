\name{hedging.iid}
\alias{hedging.iid}
\alias{HedgingIID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Value and optimal hedging strategy for a call or a put option using simulations. 
}
\description{
Computes the value of C and the optimal hedging strategy for a call or a put option on a grid at discrete time intervals, using optimal hedging and simulations. The continuous time model is
assumed to be a Levy process, so the periodic returns are i.i.d. Only the returns at the first period need to be simulated. For values of the asset not on the grid, interpolation is needed.
The optimal number of shares phi to be bought at period i-1, when the discounted price is s and the discounted value of the hedging portfolio is P, is given by
phi = (interpol1d(s,a[i,],minS,maxS)-P*rho)/s and the change in the discounted portfolio is phi * s.
At time 0, P = interpol1d(s,C[1,],minS,maxS). 
}
\usage{
hedging.iid(R,T,K,r,put,n,m,minS,maxS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{R}{Simulated iid excess periodic returns for the first period.}
\item{T}{Maturity of the option (in years).}
\item{K}{Strike price.}
\item{r}{Annual (continuous) interest rate.}
\item{put}{1 (default) for a put and 0 for a call.}
\item{n}{Number of hedging periods.}
\item{m}{Number of points of the grid.}
\item{minS}{Minimum value of the grid.}
\item{maxS}{Maximum value of the grid.}
}


\value{
\item{S}{Points on the grid at which the option is evaluated.}
\item{C}{C(i,j) represents the value of the option at period i-1 for point S(j) on the grid.}
\item{a}{a(i,j) represents a value needed to compute the optimal hedging strategy at period i-1 for point S(j) on the grid.}
\item{phi1}{phi1(j) is the initial number of shares of the asset to be bought if its price is S(j).}
\item{rho}{Constant needed for the computation of the hedging strategy.}
}
\references{
Chapter 3 of 'Statistical Methods for Financial Engineering, B. Remillard, CRC Press, (2013).
}
\author{
Bruno Remillard
}


\examples{
# Computes the price of a one-year maturity put option when hedged 5 times 
# at regular time intervals. 
# The model is assumed to be Black-Scholes with parameters mu and sigma, 
# so the excess periodic returns are Gaussian.

n = 5;         # number of hedging periods
m = 5001;      # number of points of the grid
minS = 80.0;   # minimum value of the grid
maxS = 120.0;  # maximum value of the grid
S0 = 100.0;    # initial value
K  = 100.0;    # strike price
T  = 1.0;      # maturity of the option
r  = 0.05;     # annual (continuous) rate 
put = 1;       # Put = 0 implies call!
  

#Simulation of excess periodic returns
sigma = 0.06;  # annual volatility of the returns
mu = 0.09;     # annual mean of the returns
Tp = T/n;
rp = r*Tp;
sigmap = sigma*sqrt(Tp);
Kp = K*exp(-r*T);
mup = mu*Tp-0.5*sigmap*sigmap;

#Gaussian excess returns
N = 10000;     # number of simulated returns
R = mup -rp +sigmap*rnorm(N);


# Computation
out0 = hedging.iid(R,T,K,r,put,n,m,minS,maxS)

C = out0$C;
a = out0$a;
rho = out0$rho;
S = out0$S;
phi1 = out0$phi1;

# Initial value of the option computed from interpolating C
C0 = interpol1d(S0,C[1,],minS,maxS);

# Initial value of the option computed from interpolating C
phi = (interpol1d(S0,a[1,],minS,maxS)-C0*rho)/S0;

par(mfrow=c(2,1))
plot(S,C[1,],type='s',main=bquote('Put values ' * C[0] * ' at time 0  for n' ==.(n)  ))


plot(S,phi1,type='s',main=expression('Number of shares ' *phi[1] * ' at start'))

par(new=TRUE)
C0
phi}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hedging}
\keyword{Monte Carlo}
\keyword{Call option}
\keyword{Put option}% __ONLY ONE__ keyword per line
