\name{nlsr.control}
\alias{nlsr.control}
\title{nlsr.control}
\description{
Set and provide defaults of controls for package \code{nlsr}
}
\usage{
nlsr.control(control)
}
\arguments{
\item{control}{A list of controls. 
       If missing, the defaults are provided. See below.
       If a named control is provided,  e.g., via a call
       ctrllist<- nlsr.control(japprox="jand"), then
       that value is substituted for the default of the control
       in the FULL list of controls that is returned.

NOTE: at 2022-6-17 there is NO CHECK FOR VALIDITY

The set of possible controls to set is as follows, and is returned.}
}
\value{
\item{femax}{INTEGER limit on the number of evaluations of residual function
                          Default  10000.}
                          
\item{japprox}{CHARACTER name of the Jacobian approximation to use
           Default NULL, since we try to use analytic gradient}
                 
\item{jemax}{INTEGER limit on the number of evaluations of the Jacobian
           Default 5000}
                          
\item{lamda}{REAL initial value of the Marquardt parameter
           Default 0.0001 
           Note: mis-spelling as in JNMWS, kept as historical serendipity.}
           
\item{lamdec}{ REAL multiplier used to REDUCE \code{lambda} (0 < \code{lamdec} < \code{laminc})
           Default 4, so \code{lamda <- lamda * (lamdec/laminc)} }
                           
\item{laminc}{ REAL multiplier to INCREASE \code{lambda} (1 < \code{laminc}
           Default 10}

\item{nbtlim}{ if stepredn > 0, then maximum number of backtrack loops
           (in addition to default evaluation); Default 6}

\item{ndstep}{ REAL  stepsize for numerical Jacobian approximation
           Default 1e-7}
                          
\item{offset}{ REAL A value used to test for numerical equality, i.e. \code{a} and
                  \code{b} are taken equal if \code{(a + offset) == (b + offset)}
                          Default 100.}
                          
\item{phi}{REAL Factor used to add unit Marquardt stabilization matrix in Nash
              modification of LM method.  Default 1}
                          
\item{prtlvl}{ INTEGER The higher the value, the more intermediate output is provided.
           Default 0}
                          
\item{psi}{REAL Factor used to add scaled Marquardt stabilization matrix in Nash
              modification of LM method.  Default 0} 
                          
\item{rofftest}{ LOGICAL If TRUE, perform (safeguarded) relative offset convergence test
                          Default TRUE}
                          
\item{scaleOffset}{ a positive constant to be added to the denominator sum-of-squares in
              the relative offset convergence criteria. Default 0}

\item{smallsstest}{ LOGICAL. If TRUE tests sum of squares and terminates if very small.
                          Default TRUE}
                          
\item{stepredn}{REAL Factor used to reduce the stepsize in a Gauss-Newton algorithm (Hartley's
       		method). 0 means NO backtrack. Default 0}
                          
\item{watch}{LOGICAL to provide a pause at the end of each iteration for user to monitor
           progress.   Default FALSE}
}
\author{
J C Nash 2014-7-16 revised 2022-11-22   nashjc _at_ uottawa.ca
}
