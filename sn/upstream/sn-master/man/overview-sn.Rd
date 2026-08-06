%  file sn/man/overview-sn.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2017-2022 Adelchi Azzalini
%---------------------
% HTML version: R CMD Rdconv overview.Rd --type=html -o overview.html
%
\name{overview-sn}
\docType{package}
\encoding{UTF-8}
\alias{overview-sn}
% \alias(overview}

\title{Package \pkg{sn}: overview of the structure and main commands}

\description{ 
The package provides facilities to build and manipulate probability 
distributions of the skew-normal (\acronym{SN})  and some related families, 
notably the skew-\eqn{t} (\acronym{ST}) and the \sQuote{unified skew-normal} 
(\acronym{SUN}) families. 
For the \acronym{SN}, \acronym{ST} and skew-Cauchy (\acronym{SC}) families,
statistical methods are made available for data fitting and model diagnostics,  
in the univariate and  the multivariate case.
}

\section{Two main sides}{ 
The package comprises two main sides: one side provides  facilities for the 
pertaining probability distributions; the other one deals with related
statistical methods.
Underlying formulation, parameterizations of distributions and terminology 
are in agreement with the monograph of Azzalini and Capitanio (2014).
%
%% The present document refers to version 2.0.0 of the package (March 2021).
}
 
\section{Probability side}{%

There are two layers of support for the probability distributions of interest.
At the basic level, there exist functions which follow the classical \R  
scheme for distributions.
In addition, there exists facilities to build an object which incapsulates 
a probability distribution and then certain operations can be be performed 
on such an object;
these probability objects operate according to the S4 protocol. 
The two schemes are described next.

\describe{
\item{Classical \R scheme}{%
The following functions work similary to \code{{d,p,q,r}norm} and other \R
functions for probability distributions:

\itemize{

\item skew-normal (\acronym{SN}): functions \code{{d,p,q,r}sn} for the 
  univariate case, functions  \code{{d,p,r}msn} for the multivariate case, 
  where in both cases the \sQuote{Extended skew-normal} (\acronym{ESN}) 
  variant form is included;
  
\item skew-\eqn{t} (\acronym{ST}): functions \code{{d,p,q,r}st} for the
  univariate case,  functions  \code{{d,p,r}mst} for the multivariate case;

\item skew-Cauchy (\acronym{SC}): functions \code{{d,p,q,r}sc} for the 
  univariate case, functions \code{{d,p,r}msc} for the multivariate case.
}

In addition to the usual specification of their parameters as a sequence of
individual components, a parameter set can be specified as a single \code{dp}
entity, namely a vector in the univariate case, a list in the multivariate
case; \code{dp} stands for \sQuote{Direct Parameters} (\acronym{DP}).

Another set of parameters is in use, denoted Centred Parameters (\acronym{CP}),
which are more convenient for interpretability, since they correspond to
familiar quantifies, such as the mean and the standard deviation. 
Conversion from the \code{dp} parameter set to the corresponding \acronym{CP}
parameters can be accomplished using the function \code{dp2cp},
while function \code{cp2dp} performs the inverse transformation.

The \acronym{SUN} family is mainly targeted to the multivariate context, 
and this is reflected in the organization of the pertaining functions, 
although univariate \acronym{SUN} distributions are supported.
Density, distribution function and random numbers are handled by
\code{{d,p,r}sun}. Mean value, variance matrix and Mardia's measures
of multivariate skewness and kurtosis are computed by
\code{sun{Mean,Vcov,Mardia}}.

In addition, one can introduce a user-specified density function using
\code{dSymmModulated} and \code{dmSymmModulated}, in the univariate and the
multivariate case, respectively. These densities are of the
\sQuote{symmetry-modulated} type, also called \sQuote{skew-symmetric}, where
one can specify the base density and the modulation factor with high degree of
flexibility.  Random numbers can be sampled using the corresponding functions
\code{rSymmModulated} and \code{rmSymmModulated}. In the bivariate case,
a dedicated plotting function exists.
}

\item{Probability distribution objects: \acronym{SEC} families}{%
Function \code{makeSECdistr} can be used to build a \sQuote{\acronym{SEC}
distribution} object representing a member of a specified parametric family
(among the types \acronym{SN, ESN, ST, SC}) with a given \code{dp} parameter
set. This object can be used for various operations such as plotting  or
extraction of moments and other summary quantities.  
Another way of constructing a \acronym{SEC} distribution object is via
\code{extractSECdistr} which extracts suitable components of an object 
produced by function \code{selm} to be described below.

Additional operations on these objects are possible in the multivariate case,
namely \code{marginalSECdistr} and \code{affineTransSECdistr}
for marginalization and affine trasformations. 
For the multivariate \acronym{SN} family only (but including \acronym{ESN}), 
\code{conditionalSECdistr} performs a conditioning on the values taken 
on by some components of the multivariate variable.}

\item{Probability distribution objects: the \acronym{SUN} family}{%
Function \code{makeSUNdistr} can be used to build a \acronym{SUN}
distribution object representing a member of the \acronym{SUN} 
parametric family.
This object can be used for various operations such as plotting  or
extraction of moments and other summary quantities.

Moreover there are several trasformation operations which can be performed
on a \acronym{SUN} distribution object, or two such objects in some cases:
computing a (multivariate) marginal distribution,  a conditional distribution
(on given values of some components or on one-sided intervals), an affine
trasformation, a convolution (that is, the distribution of the sum of two
independent variables), and joining two distributions under assumption of
independence.
}}
} % end of Section `Probability side'


\section{Statistics side}{%
The main function for data fitting is represented by \code{selm}, which allows
to specify a linear regression model for the location parameter, similarly 
to function \code{lm}, but assuming a \dfn{skew-elliptical} distribution
of the random component;
this explains the name \dfn{selm=(se+lm)}. Allowed types of distributions
are  \acronym{SN} (but not \acronym{ESN}), \acronym{ST} and \acronym{SC}.
The fitted distribution is univariate or multivariate, depending on the nature
of the response variable of the posited regression model. The model fitting
method is either maximum likelihood or maximum penalized likelihood; 
the latter option effectively allows the introduction of a prior distribution 
on the slant parameter of the error distribution, hence leading to a 
\sQuote{maximum a posteriori} estimate.

Once the fitting process has been accomplished, an object of class either
\dfn{selm} (for univariate response) or \dfn{mselm} (for multivariate 
response) is produced.
A number of \sQuote{methods} are available for these objects: \code{show}, 
\code{plot}, \code{summary}, \code{coef}, \code{residuals}, \code{logLik} 
and others.
For univariate \dfn{selm}-class objects, univariate and bivariate profile
log-likelihood functions can be obtained; a \code{predict} method also exists.
These methods are built following the S4 protocol; however, the user must not
be concerned with the choice of the adopted protocol (unless this is wished).

The actual fitting process invoked via \code{selm} is actually performed by a
set of lower-level procedures. These are accessible for direct call, 
if so wished, typically for improved efficiency, at the expense of a little 
additional programming effort. Similarly, functions to compute the Fisher 
information matrix are available, in the expected and the observed form 
(with some restrictions depending on the selected distribution).

The \code{extractSECdistr} function extracts the fitted \acronym{SEC} 
distribution from  \dfn{selm}-class and \dfn{mselm}-class objects, hence
providing a bridge with the probability side of the package.

The facilities for statistical work do not support the \acronym{SUN} family.
}
 
\section{Additional material}{
Additional material is available in the section 
\sQuote{User guides, package vignettes and other documentation}
accessible from the front page of the documentation. 
See especially the document \code{pkg_sn-intro.pdf}
}

\section{Author}{Adelchi Azzalini.
% Dipart. Scienze Statistiche, Università di Padova, Italia.
Please send comments, error reports \emph{et cetera} to the author, 
whose web page is \url{http://azzalini.stat.unipd.it/}.
}
 
\references{
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.
}

\keyword{multivariate}
\keyword{distribution}
\keyword{univar}
\keyword{regression}
