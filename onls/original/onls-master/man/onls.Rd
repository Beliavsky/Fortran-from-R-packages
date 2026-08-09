\name{onls}
\alias{onls}
\encoding{latin1}

\title{Orthogonal Nonlinear Least-Squares Regression}

\description{
Determines the orthogonal nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model. This is accomplished by minimizing the residual sum-of-squares of the orthogonal distances using Levenberg-Marquardt minimization in an outer loop and one-dimensional optimization for each \eqn{(x_i, y_i)} in an inner loop. See 'Details' for a description of the algorithm.
}

\usage{
onls(formula, data = parent.frame(), start, jac = NULL, 
     control = nls.lm.control(), lower = NULL, upper = NULL, 
     trace = FALSE, subset, weights, na.action, window = 12, 
     extend = c(0.2, 0.2), fixed = NULL, verbose = TRUE, ...) 
}

\arguments{
\item{formula}{a nonlinear model \code{\link{formula}} including variables and parameters. Will be coerced to a formula if necessary.}
\item{data}{an optional data frame in which to evaluate the variables in \code{formula} and \code{weights}. Can also be a list or an environment, but not a matrix.}
\item{start}{a named list or named numeric vector of starting estimates.}
\item{jac}{A function to return the Jacobian.}
\item{control}{an optional list of control settings. See \code{\link[minpack.lm]{nls.lm.control}} for the names of the settable control values and their effect.}
\item{lower}{A numeric vector of lower bounds on each parameter. If not given, the default lower bound for each parameter is set to \code{-Inf}.}
\item{upper}{A numeric vector of upper bounds on each parameter. If not given, the default lower bound for each parameter is set to \code{Inf}.}
\item{trace}{logical value indicating if a trace of the iteration progress should be printed. Default is \code{FALSE}. If \code{TRUE}, the orthogonal residual (weighted) sum-of-squares and the parameter values are printed at the conclusion of each iteration.}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{weights}{an optional numeric vector of (fixed) weights. When present, the objective function is orthogonal weighted least squares.}
\item{na.action}{a function which indicates what should happen when the data contain NAs. The default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The 'factory-fresh' default is \code{\link{na.omit}}.}
\item{window}{a numeric value for a window around the predictor values in which the optimization takes place. Used when \eqn{n > 25}.} 
\item{extend}{a two-element vector defining an extended range of predictor values, so that \eqn{(x_{0i}, y_{0i})} corresponding to \eqn{(x_i, y_i)} values at the beginning and end of the data can reside outside the original predictor range. See 'Details'.}
\item{fixed}{a logical vector defining the parameters to be fixed during optimization. See 'Examples'.}
\item{verbose}{logical. If \code{TRUE}, steps of the algorithm are printed to the console.}
\item{...}{Additional optional arguments. None are used at present.}
}

\details{
In a standard nonlinear regression setup, we estimate parameters \eqn{\boldsymbol{\beta}} in a nonlinear model \eqn{y_i = f(x_i, \beta) + \varepsilon_i}, with \eqn{\varepsilon_i \sim \mathcal{N}(0, \sigma^2)}, by minimizing the residual sum-of-squares of the vertical distances \deqn{\min_{\beta} \sum_{i=1}^{n}(y_i - f(x_i, \beta))^2}
In contrast, orthogonal nonlinear regression aims to estimate parameters \eqn{\boldsymbol{\beta}} in a nonlinear model \eqn{y_i = f(x_i + \delta_i, \beta) + \varepsilon_i}, with \eqn{\varepsilon_i, \delta_i \sim \mathcal{N}(0, \sigma^2)}, by minimizing the sum-of-squares of the orthogonal distances \deqn{\min_{\beta, \delta} \sum_{i=1}^{n}([y_i - f(x_i + \delta_i, \beta)]^2 + \delta_i^2)}
We do this by using the orthogonal distance \eqn{D_i} from the point \eqn{(x_i, y_i)} to some point \eqn{({\color{red}x_{0i}}, {\color{red}y_{0i}})} on the curve \eqn{f(x_i, \hat{\beta})} that minimizes the Euclidean distance \deqn{\min\| D_i \| \equiv \min\sqrt{(x_i - {\color{red}x_{0i}})^2 + (y_i - {\color{red}y_{0i}})^2}}
The minimization of the Euclidean distance is conducted by using an inner loop on each \eqn{(x_i, y_i)} with the \code{\link{optimize}} function that finds the corresponding \eqn{({\color{red}x_{0i}}, {\color{red}y_{0i}})} in some window \eqn{[a, b]}:
\deqn{\mathop{\rm argmin}\limits_{x_{0i} \in [a, b]} \sqrt{(x_i - {\color{red}x_{0i}})^2 + (y_i - {f({\color{red}x_{0i}}, \hat{\beta}))^2}}}

In detail, \code{onls} conducts the following steps:\cr
1) A normal (non-orthogonal) nonlinear model is fit by \code{\link[minpack.lm]{nls.lm}} to the data. This approach has been implemented because parameters of the orthogonal model are usually within a small window of the normal model. The obtained parameters are passed to 2).\cr
2) Outer loop: LM-minimization of the orthogonal distance sum-of-squares \eqn{\sum_{i=1}^N \| D_i \|^2} using \code{\link[minpack.lm]{nls.lm}}, optimization of \eqn{\boldsymbol{\beta}}.\cr
3) Inner loop: For each \eqn{(x_i, y_i)}, find \eqn{({\color{red}x_{0i}}, f({\color{red}x_{0i}}, \hat{\beta}))},  \eqn{x_{0i} \in [a, b]}, that minimizes \eqn{\| D_i \|}, using \code{\link{optimize}}. Return vector of orthogonal distances \eqn{\| \vec{D} \|}.\cr

The outer loop (\code{[minpack.lm]{nls.lm}}) scales with the number of parameters \eqn{p} in the model, probably with \eqn{\mathcal{O}(p)} for evaluating the 1-dimensional Jacobian and \eqn{\mathcal{O}(p^2)} for the two-dimensional Hessian. The inner loop has \eqn{\mathcal{O}(n)} for finding \eqn{\min \| D_i \|}. Simulations with different number of \eqn{n} showed that the processor times scale linearly.

Two arguments of this function are VERY important.\cr
1) \code{extend}: By default, it is set to \code{c(0.2, 0.2)}, which means that \eqn{(x_{0i}, y_{0i})} in the inner loop are also optimized in an extended predictor value \eqn{x} range of \eqn{[\min(x) - 0.2 \cdot \mathrm{range}(x), \max(x) + 0.2 \cdot \mathrm{range}(x)]}. This is important for the values at the beginning and end of the data, because the resulting model can display significantly different curvature if \eqn{(x_{0i}, y_{0i})} are forced to be within the predictor range, often resulting in a loss of orthogonality at the end points. See \code{\link{check_o}} for an example. If the model should be forced to be within the given predictor range, supply \code{extend = c(0, 0)}.\cr
2) \code{window}: defines the window \eqn{[x_{i - w}, x_{i + w}]} in which \code{\link{optimize}} searches for \eqn{(x_{0i}, y_{0i})} to minimize \eqn{\| D_i \|}. The default of \code{window = 12} works quite well with sample sizes \eqn{n > 25}, but may be tweaked.\cr

\bold{IMPORTANT:} If not all points are orthogonal to the fitted curve, \code{\link{print.onls}} gives a "FAILED: Only X out of Y fitted points are orthogonal" message. In this case, it is suggested to conduct a more detailed analysis using \code{\link{check_o}}. If any points are non-orthogonal, tweaking with either \code{extend} or \code{window} should in many cases resolve the problem.\cr 

The resulting orthogonal model houses information in respect to the (classical) vertical residuals as well as the (minimized Euclidean) orthogonal residuals.\cr
The following functions use the \bold{vertical} residuals:\cr
\code{\link{deviance}}\cr
\code{\link{fitted}}\cr
\code{\link{residuals}}\cr
\code{\link{logLik}}\cr

The following functions use the \bold{orthogonal} residuals:\cr
\code{\link{deviance_o}}\cr
\code{\link{residuals_o}}\cr
\code{\link{logLik_o}}\cr
}

\value{
An orthogonal fit of class \code{onls} with the following list items:\cr
\item{data}{the name of the \code{data} argument.} 
\item{call}{the matched call.}
\item{convInfo}{a list with convergence information.}
\item{model}{if \code{model = TRUE}, the model frame.}
\item{parNLS}{the converged parameters of the normal (not orthogonal!) nonlinear model, as passed to the ONLS routine.}
\item{parONLS}{the converged parameters of the not orthogonal nonlinear model, as obtained from the ONLS routine.}
\item{x0}{the \eqn{{\color{red}x_{0i}}} values from \eqn{({\color{red}x_{0i}}, {\color{red}y_{0i}})} that minimize \eqn{\|D_i\|}.}
\item{y0}{the \eqn{{\color{red}y_{0i}}} values from \eqn{({\color{red}x_{0i}}, {\color{red}y_{0i}})} that minimize \eqn{\|D_i\|}.}
\item{fittedONLS}{the fitted values \eqn{\hat{y_i}} corresponding to \eqn{x_i} and the orthogonal nonlinear model.}
\item{fittedNLS}{the fitted values \eqn{\hat{y_i}} corresponding to \eqn{x_i} and the non-orthogonal nonlinear model, as passed to the ONLS routine.}
\item{residONLS}{the (vertical) residual values \eqn{y_i - \hat{y_i}} corresponding to \eqn{x_i} and the orthogonal nonlinear model.}
\item{residNLS}{the (vertical) residual values \eqn{y_i - \hat{y_i}} corresponding to \eqn{x_i} and the non-orthogonal nonlinear model, as passed to the ONLS routine.}
\item{resid_o}{the orthogonal residual values \eqn{\| D_i \|}.}
\item{pred}{the orginal predictor values of \code{data}.}
\item{resp}{the orginal response values of \code{data}.}
\item{grad}{the numeric gradient of the model.}
\item{QR}{the QR decomposition of the model.}
\item{weights}{the weights used for fitting.}
\item{control}{the control list used, see the \code{control} argument.}
\item{ortho}{a vector of logicals for orthogonality of all points, as obtained from \code{\link{check_o}}.}
}

\author{
Andrej-Nikolai Spiess
}   

\references{
Nonlinear Perpendicular Least-Squares Regression in Pharmacodynamics.\cr
Ko HC, Jusko WJ and Ebling WF.\cr
Biopharm Drug Disp (1997), \bold{18}: 711-716.

Orthogonal Distance Regression.\cr
Boggs PT and Rogers JE.\cr
NISTIR (1990), \bold{89-4197}: 1-15.\cr
\url{https://docs.scipy.org/doc/external/odr_ams.pdf}.

User's Reference Guide for ODRPACK Version 2.01\cr
Software for Weighted Orthogonal Distance Regression.\cr
Boggs PT, Byrd RH, Rogers JE and Schnabel RB.
NISTIR (1992), \bold{4834}: 1-113.\cr
\url{https://docs.scipy.org/doc/external/odrpack_guide.pdf}.

ALGORITHM 676 ODRPACK: Software for Weighted Orthogonal Distance Regression.\cr
Boggs PT, Donaldson JR, Byrd RH and Schnabel RB.\cr
ACM Trans Math Soft (1989), \bold{15}, 348-364.\cr
\doi{10.1145/76909.76913}.
}

\examples{
## 1A. The DNase data from 'nls',
## use all generic functions.
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
print(mod1)
plot(mod1)
summary(mod1)
predict(mod1, newdata = data.frame(conc = 6))
logLik(mod1)
deviance(mod1)
formula(mod1)
weights(mod1)
df.residual(mod1)
fitted(mod1)
residuals(mod1)
vcov(mod1)
coef(mod1)

DNase2 <- DNase1
DNase2$conc <- DNase2$conc * 2
mod2 <- update(mod1, data = DNase2)
print(mod2)

\donttest{
## 1B. Same model as above, but using the restricted
## predictor range which results in non-orthogonality
## of some points.
onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
     data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
     extend = c(0, 0))

## 2. Example from odrpack_guide.pdf, 2.C.i, pages 39ff.
x <- c(0, 0, 5, 7, 7.5, 10, 16, 26, 30, 34, 34.5, 100)
y <- c(1265, 1263.6, 1258, 1254, 1253, 1249.8, 1237, 1218, 1220.6, 
1213.8, 1215.5, 1212)
DAT <- data.frame(x, y)
mod3 <- onls(y ~ b1 + b2 * (exp(b3 * x) -1)^2, data = DAT, 
             start = list(b1 = 1500, b2 = -50, b3 = -0.1))
deviance_o(mod3) # 21.445 as in page 47
coef(mod3) # 1264.65481/-54.01838/-0.08785 as in page 48

## 3. Example from Algorithm 676: ODRPACK, page 355 + 356.
x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100, 105)
y <- c(4.14, 8.52, 16.31, 32.18, 64.62, 98.76, 151.13, 224.74, 341.35, 
       423.36, 522.78, 674.32, 782.04, 920.01)
DAT <- data.frame(x, y)
mod4 <- onls(y ~ b1 * 10^(b2 * x/(b3 + x)), data = DAT, 
             start = list(b1 = 1, b2 = 5, b3 = 100))
coef(mod4) # 4.4879/7.1882/221.8383 as in page 363
deviance_o(mod4) # 15.263 as in page 363

## 4. Example with bounds from simple_example.f90
## in https://www.netlib.org/toms/869.zip.
x <- c(0.982, 1.998, 4.978, 6.01)
y <- c(2.7, 7.4, 148.0, 403.0)
DAT <- data.frame(x, y)
mod5 <- onls(y ~ b1 * exp(b2 * x), data = DAT, 
            start = list(b1 = 2, b2 = 0.5), 
            lower = c(0, 0), upper = c(10, 0.9))
coef(mod5) # 1.4397(1.6334)/0.9(0.9) ## Different to reference!
deviance_o(mod5) # 0.1919 (0.2674) => but lower RSS than original ODRPACK!

## 5. Example with a fixed parameter
## => Asym = 3.
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod6 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1), 
             fixed = c(TRUE, FALSE, FALSE))
print(mod6)

## 6. Example to show that one can even conduct
## linear orthogonal regression (Deming regression).
## Comparison to XLstat
## https://help.xlstat.com/6650-run-deming-regression-compare-methods-excel
x <- c(9.8, 9.7, 10.7, 10.9, 12.4, 12.5, 12.8, 12.8, 12.9, 13.3, 
       13.4, 13.5, 13.7, 14.9, 15.2, 15.5)
y <- c(10.1, 11.4, 10.8, 11.3, 11.8, 12.1, 12.3, 13.6, 14.2, 14.4,
       14.6, 15.3, 15.5, 15.8, 16.2, 16.5)
DAT <- data.frame(x, y)
mod7 <- onls(y ~ a + b * x, data = DAT, 
            start = list(a = 2, b = 3))
print(mod7) ## -1.909/1.208 as on webpage
plot(mod7)
}
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

