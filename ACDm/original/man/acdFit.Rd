\name{acdFit}
\alias{acdFit}
\title{
ACD (Autoregressive Conditional Duration) Model Fitting
}

\description{
This function estimates various ACD models with various assumed error term distributions,
using Maximum Likelihood Estimation.

The currently available models (conditional mean specifications) are: Standard ACD,
Log-ACD (two alternative specifications), AMACD, ABACD, BACD, BCACD, AACD, EXACD,
SNIACD, LSNIACD, TACD, and TAMACD.

The currently available distributions are: Exponential (also used for QML), Weibull,
Burr, generalized Gamma, generalized F, Q-Weibull, and Mixture inverse Gaussian.
}


\usage{
acdFit(durations = NULL, model = "ACD", dist = "exponential", 
    order = NULL, startPara = NULL,  dailyRestart = FALSE, optimFnc = "optim",
    method = "Nelder-Mead", output = TRUE, bootstrapErrors = FALSE, 
    forceErrExpec = TRUE, fixedParamPos = NULL, bp = NULL, 
    exogenousVariables = NULL, optimFncArgs = list(),
    control = list())
}

\arguments{

  \item{durations}{
either (1) a data frame including, at least, a column named 'durations' or 'adjDur' (for adjusted durations), or (2) a vector of durations
}
  \item{model}{
the conditional mean model specification. Must be one of \code{"ACD"}, \code{"LACD1"}, \code{"LACD2"}, \code{"AMACD"}, \code{"BACD"}, 
\code{"ABACD"}, \code{"SNIACD"} or \code{"LSNIACD"}. See 'Details' for detailed model specification. 
}
  \item{dist}{
the assumed error term distribution. Must be one of \code{"exponential"}, \code{"weibull"}, \code{"burr"}, \code{"gengamma"}, \code{"genf"}, \code{"qweibull"}, \code{"mixqwe"}, \code{"mixqww"}, or \code{"mixinvgauss"}. See 'Details' for detailed model specification. 
}
  \item{order}{
a vector detailing the order of the particular ACD model. For example an ACD(p, q) specification should have \code{order = c(p, q)}.
}
  \item{startPara}{
a vector with parameter values to start the maximization algorithm from. Must be in the correct order according to the model specification (see Details).
}
  \item{dailyRestart}{
if \code{TRUE} the conditional duration will start fresh every new trading day. Can only be used if the \code{durations} arguments included the clock time of the durations, or if the \code{time} argument was provided.
}
  \item{optimFnc}{
Specifies which optimization function to use for the estimation. \code{"optim"}, \code{"nlminb"}, \code{"solnp"}, and \code{"optimx"} are available.
}
  \item{method}{
Argument passed to the optimization function if \code{optimFnc = "optim"} or \code{optimFnc = "optimx"} were chosen. Specifies the optimization algorithm. See the help files for \link{optim}, \link{nlminb} or \link[Rsolnp]{solnp}.
}
  \item{output}{
if \code{FALSE} the estimation results won't be printed.
}
  \item{bootstrapErrors}{
if \code{TRUE} the standard errors will be computed by using bootstrap simulations. Currently only works with the standard ACD model.
}
  \item{forceErrExpec}{
if \code{TRUE} the expectation of the error terms' distribution will be forced to be 1, otherwise the distribution parameter specifying the mean will be set to 1 to ensure identification.
}
  \item{fixedParamPos}{
a logical vector of \code{TRUE} and \code{FALSE}. Can only be used if the argument \code{startPara} was provided, and should be of the same length. Each element represents the respective start parameter and if \code{TRUE}, this parameter will be held fixed when estimating the other parameters.
}
  \item{bp}{
a numeric vector of breakpoints, used if the \code{model} is any of \code{"SNIACD", "LSNIACD", "TACD", "TAMACD"}. The default is a single breakpoint at 1. 
}
  \item{exogenousVariables}{
specifies the columns in the \code{durations} data.frame that should be used as exogenous variables when fitting the model. Must be a vector, either with the column positions or the names of the columns. It is highly recommended to standardize the exogenous variables before running the estimation. The models use the current exogenous variable and not their lagged values. 
}
  \item{optimFncArgs}{
a list of values passed to the optimaization function. 
}
  \item{control}{
a list of control values:
\describe{
   \item{maxit}{
	maximum number of iterations performed by the numerical maximization algorithm.
   }

   \item{trace}{
	An integer. If set to a non-zero value, the parameter values at each call of the log-likelihood function during optimization will be recorded, and the resulting search path of the maximum-likelihood estimate will be plotted. This argument is also passed on to the underlying optimizer; see the help files for \link{optim}, \link{nlminb} or \link[Rsolnp]{solnp} for details. }

   \item{B}{
	number of bootstrap samples
   }
}
}
}
\details{

The \code{startPara} argument is  a vector of the parameter values to start from. The length of the vector naturally depends on the \code{model} and \code{distribution}. The first elements represent the model parameters, and the last elements the distribution parameters. For example for an ACD(1,1) with Weibull errors the first 3 elements are \eqn{\omega, \alpha_1, \beta_1} for the model, and the last is \eqn{\gamma} for the Weibull distribution.

The family of ACD models are
\deqn{x_i = \mu_i \epsilon_i,}
where different specifications of the conditional mean duration \eqn{\mu_i} and the error term \eqn{\epsilon_i} give rise to different models as shown below. 

When exogenous variables are used, they are added in the form of 
\deqn{\sum_{j=1}^{k} \xi_j z_j} 
to the right hand side of the equations, where \eqn{z_j} are the exogenous variables.

\strong{Conditional mean duration \eqn{\mu_i} specifications according to the \code{model} argument:}
\cr

\bold{ACD(p, q) specification:} (Engle and Russell,  1998)
\deqn{\mu_i = \omega + \sum_{j=1}^{p} \alpha_j x_{i-j} + \sum_{j=1}^{q} \beta_j \mu_{i-j}}
The element order of the \code{startPara} vector is \eqn{(\omega, \alpha_j...,\beta_j...)}.
\cr

\bold{LACD1(p, q):} (Bauwens and Giot,  2000)
\deqn{\ln\mu_i = \omega + \sum_{j=1}^{p} \alpha_j \ln \epsilon_{i-j} + \sum_{j=1}^{q} \beta_j \ln \mu_{i-j}}
The element order of the \code{startPara} vector is \eqn{(\omega, \alpha_j...,\beta_j...)}.
\cr

\bold{LACD2(p, q):} (Lunde,  1999)
\deqn{\ln\mu_i = \omega + \sum_{j=1}^{p} \alpha_j \epsilon_{i-j} + \sum_{j=1}^{q} \beta_j \ln \mu_{i-j}}
The element order of the \code{startPara} vector is \eqn{(\omega, \alpha_j...,\beta_j...)}.
\cr

\bold{AMACD(p, r, q) (Additive and Multiplicative ACD):} (Hautsch , 2012)
\deqn{\mu_i = \omega + \sum_{j=1}^{p} \alpha_j x_{i-j} + \sum_{j=1}^{r} \nu_j \epsilon_{i-j} + \sum_{j=1}^{q} \beta_j \mu_{i-j}}
The element order of the \code{startPara} vector is \eqn{(\omega, \alpha_j...,\nu_j...,\beta_j...)}.
\cr

\bold{ABACD(p, q) (Augmented Box-Cox ACD):} (Hautsch,  2012)
\deqn{\mu_i^{\delta_1} = \omega + \sum_{j=1}^{p} \alpha_j \left( |\epsilon_{i-j}-\nu|+c|\epsilon_{i-j}-\nu| \right)^{\delta_2} + \sum_{j=1}^{q} \beta_j \mu_{i-j}^{\delta_1}}
The element order of the \code{startPara} vector is \eqn{(\omega, \alpha_j..., \beta_j..., c, \nu, \delta_1, \delta_2)}.
\cr

\bold{BACD(p, q) (Box-Cox ACD):} (Hautsch,  2003)
\deqn{\mu_i^{\delta_1} = \omega + \sum_{j=1}^{p} \alpha_j \epsilon_{i-j}^{\delta_2} + \sum_{j=1}^{q} \beta_j \mu_{i-j}^{\delta_1}}
The element order of the \code{startPara} vector is \eqn{(\omega, \alpha_j..., \beta_j...)}.
\cr

\bold{EXACD(p, q) (Exponential ACD):} (Dufour and Engle, 2000)
\deqn{\ln\mu_i = \omega + \sum_{j=1}^{p} \left( \alpha_j \epsilon_{i-j} + \delta_j |\epsilon_{i-j} - 1| + \sum_{j=1}^{q} \beta_j \mu_{i-j} \right)}
The element order of the \code{startPara} vector is \eqn{(\omega, \alpha_j..., \beta_j..., \delta_j...)}.
\cr

\bold{SNIACD(p,  q, M) (Spline News Impact ACD):} (Hautsch,  2012, with a slight difference)
\deqn{\mu_i = \omega + \sum_{j=1}^{p} (\alpha_{j-1}+c_0) \epsilon_{i-j} + \sum_{j=1}^{p} \sum_{k=1}^{M} (\alpha_{j-1}+c_k)1_{(\epsilon_{i-j} \le  \bar{\epsilon_k})}+\sum_{j=1}^{q} \beta_j \mu_{i-j},}
where \eqn{1_{()}} is an indicator function and \eqn{\alpha_0=0}. \cr
The element order of the \code{startPara} vector is \eqn{(\omega, c_k..., \alpha_j..., \beta_j...)} (The number of \eqn{\alpha}-parameters are \eqn{p-1}).
\cr

\strong{The distribution of the error term \eqn{\epsilon_i} specifications according to the \code{dist} argument:}
\cr \cr
\bold{Exponential distribution, \code{dist = "exponential"}:}

\deqn{f(\epsilon)=\exp(-\epsilon)}

\bold{Weibull distribution, \code{dist = "weibull"}:}

\deqn{f(\epsilon)=\theta \gamma \epsilon^{\gamma-1}e^{-\theta \epsilon^{\gamma}} ,}
where \eqn{\theta=[\Gamma(\gamma^{-1}+1)]^{\gamma}} if \code{forceErrExpec = TRUE}.

Note that this is a diffrent parameterization then the one in \code{stats::dweibull}. While the \code{shape} parameter in \code{stats::dweibull} is the same as \eqn{\gamma}, the \code{scale} parameter is equal to \eqn{\theta^{-\gamma}}.

\bold{Burr distribution, \code{dist = "burr"}:}

\deqn{f(\epsilon)= \frac{\theta \kappa \epsilon^{\kappa-1}}{(1+\sigma^2 \theta \epsilon^{\kappa})^{\frac{1}{\sigma^2}+1}},}
where,
\deqn{\theta= \sigma^{2 \left(1+\frac{1}{\kappa}\right)} \frac{\Gamma \left(\frac{1}{\sigma^2}+1\right)}{\Gamma \left(\frac{1}{\kappa}+1\right) \Gamma \left(\frac{1}{\sigma^2}-\frac{1}{\kappa}\right)},}
if \code{forceErrExpec = TRUE}. \cr
The element order of the \code{startPara} vector is \eqn{(model parameters, \kappa, \sigma^2)}.
\cr

\bold{Generalized Gamma distribution, \code{dist = "gengamma"}:}

\deqn{f(\epsilon)=\frac{\gamma \epsilon^{\kappa \gamma - 1}}{\lambda^{\kappa \gamma}\Gamma (\kappa)}\exp \left\{{-\left(\frac{\epsilon}{\lambda}\right)^{\gamma}}\right\}}
where \eqn{\lambda=\frac{\Gamma(\kappa)}{\Gamma(\kappa+\frac{1}{\gamma})}} if \code{forceErrExpec = TRUE}.
The element order of the \code{startPara} vector is \eqn{(model parameters, \kappa, \gamma)}.
\cr

\bold{Generalized F distribution, \code{dist = "genf"}:}
\deqn{f(\epsilon)= \frac{\gamma \epsilon^{\kappa \gamma -1}[\eta+(\epsilon/\lambda)^{\gamma}]^{-\eta-\kappa}\eta^{\eta}}{\lambda^{\kappa \gamma}B(\kappa,\eta)},}
where \eqn{B(\kappa,\eta)=\frac{\Gamma(\kappa)\Gamma(\eta)}{\Gamma(\kappa+\eta)}}, and if \code{forceErrExpec = TRUE},
\deqn{\lambda=\frac{\Gamma(\kappa)\Gamma(\eta)}{\eta^{1/\gamma}\Gamma(\kappa+1/\gamma)\Gamma(\eta-1/\gamma)}.}
\cr
The element order of the \code{startPara} vector is \eqn{(model parameters, \kappa,  \eta, \gamma)}.
\cr

\bold{q-Weibull distribution, \code{dist = "qweibull"}:}
\deqn{f(\epsilon) = (2-q)\frac{a}{b^a} \epsilon^{a-1} \left[1-(1-q)\left(\frac{\epsilon}{b}\right)^a\right]^{\frac{1}{1-q}}}
where if \code{forceErrExpec = TRUE},
\deqn{b = \frac{(q-1)^{\frac{1+a}{a}}}{2-q}\frac{a\Gamma(\frac{1}{q-1})}{\Gamma(\frac{1}{a}) \Gamma(\frac{1}{q-1}-\frac{1}{a}-1)}.}
\cr
The element order of the \code{startPara} vector is \eqn{(model parameters, a,  q)}.
\cr
}
\value{

a list of class \code{"acdFit"} with the following slots: 

\item{durations}{
        the durations object used to fit the model.
    }
\item{muHats}{
        a vector of the estimated conditional mean durations
    }
\item{residuals}{
        the residuals from the fitted model, calculated as durations/mu
    }
\item{model}{
        the model for the conditional mean durations
    }
\item{order}{
        the order of the model
    }
\item{distribution}{
        the assumed error term distribution
    }
\item{distCode}{
        the internal code used to represent the distribution
    }
\item{mPara}{
        a vector of the estimated conditional mean duration parameters
    }
\item{dPara}{
        a vector of the estimated error distribution parameters
    }
\item{Npar}{
        total number of parameters
    }
\item{goodnessOfFit}{
        a data.frame with the log likelihood, AIC, BIC, and MSE calculated as the mean squared deviation of the durations and the estimated conditional durations.
    }
\item{parameterInference}{
        a data.frame with the estimated coefficients and their standard errors and p-values
    }
\item{forcedDistPara}{
        the value of the unfree distribution parameter. If \code{forceErrExpec = TRUE} were used, this parameter is a function of the other distribution parameters, to force the mean of the distribution to be one. Otherwise the parameter was fixed at 1 to ensure identification.
    }
\item{comments}{
        
    }
\item{hessian}{
        the numerical hessian of the log likelihood evaluated at the estimate
    }
\item{N}{
        number of observations
    }
\item{evals}{
        number of log-likelihood evaluations needed for the maximization algorithm
    }
\item{convergence}{
        if the maximization algorithm converged, this value is zero. (see the help file \link{optim}, \link{nlminb} or \link[Rsolnp]{solnp})
    }
\item{estimationTime}{
        time required for estimation
    }
\item{description}{
        who fitted the model and when 
    }
\item{robustCorr}{
        only available for QML estimation (choosing the exponential distribution) for the standard ACD(p, q) model. The robust correlation matrix of the parameter estimates. 
    }

}
\references{
Bauwens, L., and P. Giot (2000) 
\emph{The logarithmic ACD model: an application to the bid-ask quote process of three NYSE stocks.} Annales d'Economie et de Statistique, 60, 117-149.


Dufour, A., & Engle, R. F. (2000)
\emph{The ACD Model: Predictability of the Time Between Consecutive Trades} (Discussion Papers in Finance: 2000-05). ISMA Centre, University of Reading.
    

Engle R.F, Russell J.R. (1998)
\emph{Autoregressive Conditional Duration: A New Model for Irregularly Spaced Transaction Data}, Econometrica, 66(5): 1127-1162.
    
Grammig, J., and Maurer, K.-O. (2000) 
\emph{Non-monotonic hazard functions and the autoregressive conditional duration model.} Econometrics Journal 3: 16-38.

Hautsch, N. (2003) 
\emph{Assessing the Risk of Liquidity Suppliers on the Basis of Excess Demand Intensities.} Journal of Financial Econometrics (2003) 1 (2): 189-215

Hautsch, N. (2012) 
\emph{Econometrics of Financial High-Frequency Data.} Berlin, Heidelberg: Springer.

Lunde, A. (1999): 
\emph{A generalized gamma autoregressive conditional duration model,} Working paper, Aalborg University.
}
\author{
 Markus Belfrage
}

\examples{

data(adjDurData)
fitModel <- acdFit(durations = adjDurData, model = "ACD", 
            dist = "exponential", order = c(1,1), dailyRestart = TRUE)



}

