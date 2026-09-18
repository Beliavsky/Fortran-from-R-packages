\name{pbs}
\alias{pbs}
\alias{pbs.pbs}
\title{
Periodic B-Spline Basis for Polynomial Splines
}
\description{
Generate the periodic B-spline basis matrix for a polynomial spline.
}
\usage{
pbs(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
Boundary.knots = range(x), periodic = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the predictor variable.  Missing values are allowed.}

  \item{df}{
          degrees of freedom; one can specify 'df' rather than 'knots';
          'pbs()' then chooses 'df - 1' knots at suitable quantiles
          of 'x' (which will ignore missing values).
  }

  \item{knots}{
          the _internal_ breakpoints that define the spline.
	  The number of internal knots
          must be greater than or equal to degree
          polynomial regression. See also
          'Boundary.knots'.
  }

  \item{degree}{
          degree of the piecewise polynomial-default is 3 for cubic
          splines.
  }

  \item{intercept}{
          if 'TRUE', an intercept is included in the basis; default is
          'FALSE'.
  }

  \item{Boundary.knots}{
          boundary points at which to set the period of the perodic 
	  B-spline basis(default the range of the data). If both 'knots' 
	  and 'Boundary.knots' are supplied, the basis parameters do not
          depend on 'x'. Data CAN NOT be extended beyond 'Boundary.knots'.
	  Typical Bourday knots are start and end values of period.
  }

  \item{periodic}{
          if "TRUE", periodic basis is generated. Default is TRUE. 
	  This is disabled. DO NOT USE. 
  }
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
     A matrix of dimension 'length(x) * (df)', where either 'df' was
     supplied or if 'knots' were supplied, 'df = length(knots) +
     intercept'.  Attributes are returned that  correspond to the
     arguments to 'pbs', and explicitly give the 'knots',
     'Boundary.knots' etc for use by 'predict.pbs()'.

     'pbs()' is based on the function 'spline.des()' in package splines. 
     It generates a basis matrix for representing the family of piecewise 
     polynomials with the specified interior knots and degree, evaluated 
     at the values of 'x'.  A primary use is in modeling formulas to 
     directly specify a piecewise polynomial term in a model.
}

\references{
     Package splines in Base R. pbs is actually modified from bs function 
}

%%\author{
%%}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    'bs', 'ns', 'predict.pbs'
}

\examples{
     require(stats); require(graphics); require(splines)
     x = seq(1,628)/100
     z = rep(seq(1, 314)/100, 2)
     
     
     pbs(x, df = 5, Boundary.knots = c(0, 2*pi))
     pbs(x, knots=c(pi/2, pi, pi*3/2), Boundary.knots = c(0, 2*pi))
     
     #### example of one periodic functions
     y= sin(x) + cos(2*x) + 
        rnorm(628, 0, 0.1) ## x has a period of 2*pi

     ## df method, need to use large enough df to get a better fit.  
     ## May use max loglik to choose optimal df
     summary( fm1 <- lm(y ~ pbs(x, df = 10, Boundary.knots = c(0, 2*pi))) )
     plot(x, y, xlab = "x", ylab = "sin(x)", pch="x", cex=.5)
     
     lines(x, predict(fm1, data.frame(x=x, z=z)), col='blue')
     lines(x, sin(x) + cos(2*x), col='red') 
     
     ## knots methods, usually selected at turning points
     summary( fm2 <- lm(y ~ pbs(x, knots=c(pi/2, pi, pi*3/2), 
              Boundary.knots = c(0, 2*pi)))
     )
     plot(x, y, xlab = "x", ylab = "sin(x)", pch="x", cex=.5)
     
     lines(x, predict(fm2, data.frame(x=x, z=z)), col='blue')
     lines(x, sin(x) + cos(2*x), col='red') 

     #### example of two periodic functions
     x0 = seq(1,628, by=4)/100
     z0 = seq(1, 314, by=3)/100
     x = rep(x0, each=length(z0))
     z = rep(z0, length(x0))
     y = sin(x) + cos(2*z) + 
        rnorm(length(x), 0, 0.1) ## x has a period of 2*pi and z of pi

     summary( fm3 <- lm(y ~ pbs(x, df = 5, Boundary.knots = c(0, 2*pi))+ 
                            pbs(z, df = 5, Boundary.knots = c(0, pi))) 
     )

     plot(sin(x) + cos(2*3), predict(fm3, data.frame(x=x, z=3)))
     summary(sin(x) + cos(2*3)- predict(fm3, data.frame(x=x, z=3)))
     ## End(Not run)
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%%\keyword{}
%%\keyword{}% __ONLY ONE__ keyword per line
