------------------------------------------------------------------------------
version 1.1.10
------------------------------------------------------------------------------

  * fix clang15 warnings
  * fix latex escape note
  * fix version info from new.md note

------------------------------------------------------------------------------
version 1.1.9
------------------------------------------------------------------------------

  * reset license from GPL-2 to the original GPL (>=2)

------------------------------------------------------------------------------
version 1.1.8
------------------------------------------------------------------------------

  * fixed typo in C code causing i386 error in rmutil::int()

------------------------------------------------------------------------------
version 1.1.7
------------------------------------------------------------------------------

  * added Author information of TOMS614/INTHP
  * fixed void / SEXP registration issue (LTO / valgrind additional issues)
  * fixed 0-sized array warning 
  * Fixed call_R Issue https://github.com/swihart/rmutil/issues/11
  * Updated links to https: where appropriate

------------------------------------------------------------------------------
version 1.1.6
------------------------------------------------------------------------------

  * Fixed call_R Issue https://github.com/swihart/rmutil/issues/11
  * Updated links to https: where appropriate


------------------------------------------------------------------------------
version 1.1.5
------------------------------------------------------------------------------

  * Fixed Issue https://github.com/swihart/rmutil/issues/10 -- updated BetaBinomial 
  help documentation to relate `m` and `s` to the standard/traditional `alpha` and `beta`
  parameterization.  Thanks to  @sjenkins20 on github for the suggestion!
  * Fixed Issue https://github.com/swihart/rmutil/issues/9 -- Condition has length > 1 and only the first element will be used warning in `int`. Thank you @ellessenne on github for chasing this one down.

------------------------------------------------------------------------------
version 1.1.4
------------------------------------------------------------------------------


  * In response to K. Hornik email regarding the stringsAsFactors=FALSE default in upcoming R versions, I did the following 3 edits:
  * add `stringsAsFactors=TRUE` to the v <- data.frame() call in the example of `tvctomat.Rd`
  * add `stringsAsFactors=TRUE` to the data.frame() call in the function `tvctomat()` (around line 1160 in `objectrm.r`): (`oldtvcov$tvcov <- data.frame(oldtvcov$tvcov,tvcv, stringsAsFactors = TRUE)`
  * turn `tvcv <- as.data.frame(as.vector(t(as.matrix(tvcov))))` into `tvcv <- data.frame(as.character(as.vector(t(as.matrix(tvcov)))),stringsAsFactors=TRUE)`  in the function `tvctomat()` (around line 937 of `objectrm.r`)
  * Issue https://github.com/swihart/rmutil/issues/8. 

------------------------------------------------------------------------------
version 1.1.3
------------------------------------------------------------------------------

  * Implemented better fixes for https://github.com/swihart/rmutil/issues/5 
    (Thanks @hennerw!) where values below the support will give 0 and values 
    above support give 1 for pbetabinom and dbetabinom as the last fix in 
    1.1.2 could give erroneous results.  
    Let me know if other functions need this
    functionality -- I'd be happy to update.
  * Removed local generic `weights` as per personal communication from 
    K Hornik

------------------------------------------------------------------------------
version 1.1.2
------------------------------------------------------------------------------

  * Implemented https://github.com/swihart/rmutil/issues/5 (Thanks @hennerw!) 
    where values below the support will give 0 and values above support give 1
    for pbetabinom and dbetabinom.  Let me know if other functions need this
    functionality -- I'd be happy to update.

------------------------------------------------------------------------------
version 1.1.1
------------------------------------------------------------------------------

  * fixed an `_R_CHECK_LENGTH_1_CONDITION_=true` error
  * fixed `_R_S3_METHOD_LOOKUP_BASEENV_AFTER_GLOBALENV_=true` errors.  
  * changed Example in DataMethods.Rd to `glm(y~x+z, data=as.data.frame(reps))`


------------------------------------------------------------------------------
version 1.1.0
------------------------------------------------------------------------------

Major changes

Passed CRAN checks and is back on CRAN.  Please see github page 
https://github.com/swihart/rmutil to see the changes required 
to pass CRAN checks.
I'll try to document the changes henceforth here and at 
https://github.com/swihart/rmutil/issues.


------------------------------------------------------------------------------
version 1.0
------------------------------------------------------------------------------
30.11.10 (growth)
elliptic: twins option added for model using covfn with covariance matrix
diagonal constant
28.11.10
elliptic: added an error check when covfn used
15.2.10 (rmutil)
changed s<0 to s<=0 in qlaplace & rlaplace (thanks to Peter Ehlers)
18.11.09 (repeated)
removed a redundant line in gausscop.r that now produced an error
	(thanks to Patrick Lindsey)
7.4.09
removed extra } in biv.binom.Rd (thanks to Christopher Marcum)
20.10.08
discrete q functions: changed trunc to round (thanks to Frederic Gosselin)
3.7.08 (gnlm)
fit.dist: corrected check for negative values with Laplace, Cauchy,
	 and Student t plus error in counts (f -> ni) for Laplace
	 (thanks to Michael Anyadike-Danes)
24.10.07
fnenvir: changed way "=" is handled in gsub() because of error since R2.5.0
8.10.07 (event, gnlm, growth, repeated)
changed typsiz to typsize in nlm() throughout
11.7.07
romberg.c: added missing R.h (thanks to Olivia Lau)
8.2.07
print out name of response variable in elliptic, bnlr, gnlr, gnlr3,
	 gnlmm, gnlmm3, and fmr (thanks to Patrick Lindsey)
qsimplex: corrected search interval updates
27.9.06
qhjorth, qinvgauss, qginvgauss, qboxcox: changed lower limit of search
	 from 0.001 to .Machine$double.xmin (thanks to Franco Mendolia)
8.12.05 (rmutil, repeated, event)
minor modifications to Fortran for F95 compatibility (thanks to Jan de Leeuw)
22.11.05
finterp, objectrm: added na.action=NULL in calls to model.frame
	 (default is na.omit !!!!!) (thanks to Patrick Lindsey)
30.9.05
elliptic: corrected calculation of number of parameters for builtin
	  logistic function (thanks to Tom Van Dooren)
1.8.05
qbetabinom: changed trunc() to round() (thanks to Elias Krainski)
mprofile, iprofile: added check that times are available (thanks to
	Patrick Lindsey)
6.7.05
ksurvb, ksurvg, kcountb: added break; after default: in C code to
	satisfy some compilers (thanks to Patrick Lindsey)
30.6.05
finterp: correction for change in functioning of match()
gnlm family: added coef.gnlm() and vcov.gnlm() (thanks to Bendix Carstensen)
25.4.05 (stable)
rstable: eliminate production of NAs (thanks to Zhu Wang)
26.1.05
finterp: fixed a bug when >= or <= is used in a formula (thanks to
	 Eliot McIntire)
1.11.04
gnlmm3: generalized nonlinear mixed models for three-parameter distributions
28.9.04
catmiss: removed codes() from example in help (thanks to Kjetil Brinchmann)
21.9.04
finterp: fixed if test to remove occasional warning (thanks to Ken Knoblauch)
17.9.04
gnlmix: removed erroneous printing that distribution is censored for
	binomial (thanks to Ken Knoblauch)
28.7.04
gnlmix, hnlmix: fixed printing of results when nonlinear function
	contains a linear part (thanks to Ken Knoblauch)
2.7.04
tvctomat: fixed warning message on removing tm (thanks to Patrick Lindsey)
1.6.04
glmm: changed print.summary.glmm to work under R1.9 (thanks to Spencer Graves)
5.4.04
fnenvir: fixed obscure error when linear used in a function (thanks to
	Ken Knoblauch)
help: corrected truncation of usage formula for certain functions
	(thanks to Patrick Lindsey)
9.1.04
fitdist: fixed typo that stopped geometric distribution from working
6.1.04
ordglm: changed tapply to capply because of change to former (thanks
	to Andrew Criswell)
9.12.03
corgram: start abscissa at 0 for PACF
fnenvir: fixed grep for checking redundant parameters
bnlr, fmr, gnlmm, gnlr, gnlr3, int, nlr, nordr, read.list, stable.mode:
	fixed if() on vector
14.11.03
readdna, survkit, glmm, gnlmm, objectrm, readrm: removed obsolete
	codes() function (thanks to Ken Knoblauch)
carma: give error when ccov=~1 used (thanks to Patrick Lindsey)
21.8.03
elliptic: corrected print function when the dispersion function depends
	on the location function (thanks to Gabrielle Kelly)
31.7.03
hnlmix: corrected options for dnbinom (thanks to Jagat Sheth)
30.6.03
dftorep: corrected check for ordered times to allow for two levels of nesting
25.5.03
ordglm: added a data argument (thanks to Kosuke Imai)
13.5.03
ordglm: corrected order of printing of standard errors (thanks to Kosuke Imai)
25.4.03
gnlr, gnlr3, gnlmm, fmr, nordr: changed test for environment because
	the value returned by parent.frame() has a class in R1.7
22.4.03
cphidden: a function to locate a changepoint in continuous time using
	a two-state hidden Markov model
9.4.03
biv.binom: corrected degrees of freedom printed (thanks to Goran Brostrom)
12.2.03
restovec: fixed handling of delta when the response is a list
16.1.03
kalsurv: fixed typo in print.kalsurv (thanks to Anthony Gichangi)
4.12.02
int: changed eps
2.12.02
fit.dist: added Laplace distribution
1.12.02
glmm: added error message if (codes of) nesting variable not consecutively
	numbered (thanks to Renaud Lancelot)
27.11.02
fit.dist: changed Weibull parametrisation so that mu refers to y and
	not to y^alpha
22.11.02
fit.dist: added Cauchy and Student t distributions
	use (log) density functions instead of writing formulae
18.11.02
fit.dist: added beta-binomial distribution
16.11.02
fit.dist: corrected error in calculation of log likelihood when censor=T
14.11.02
fit.dist: corrected error in calculation of fitted values for zeta distribution
31.10.02
int2: added default limits (thanks to Patrick Lindsey)
8.9.02 (repeated)
gar: corrected recursive fitted values when binomial (thanks to
	Patrick Lindsey)
4.9.02
gausscop: exponential distribution now works (thanks to Patrick Lindsey)
30.8.02
restovec: modified checks for nesting in lists and allow nesting to be
	supplied separately when a list is given (thanks to Patrick Lindsey)
gausscop: for positive-valued distributions, returned predicted values
	without transforming by log link function (thanks to Patrick Lindsey)
18.7.02
ehr: addition checks of data supplied to this suite of functions
rs3: fixed typo
marg.hom: added checks on data
16.7.02
chidden.r, hidden.r: corrected negative binomial check so that 0
	responses are allowed (thanks to Ben Cooper)
10.7.02
modified man pages for changed arguments to rgamma function
rmutil: created dist.h file
11.6.02
hnlmix: corrected AIC for penalty constraint (was too large by one)
	changed calculation of multiplicative random effects
23.5.02
rmutil: added [pdqr]twosidedpower distribution
	added log option to all density (d) functions
gar, gnlr, gnlmix, gnlmm, hnlmix: added two-sided power distribution
gnlr: user-supplied likelihood function works again (thanks to Martin Liermann)
finterp, fnenvir: added option to allow any response to be a covariate
9.5.02
hnlmix: recursive fitted values available
ordglm: fixed error that PearsRes not defined when individual data are
	supplied (thanks to Kamal Desai)
6.5.02
gnlmix, hnlmix: added inverse gamma mixture distribution
gnlmix: handles censored data
gnlmm: finds nesting variable when repeated environment is specified
5.5.02
finterp: modified so that as.factor(times) works
30.4.02
hnlmix: nonlinear random effects models using a modified Lee and
	Nelder h-likelihood
gnlr: modified check on location parameters for Levy distribution
      added check that double Poisson, multiplicative Poisson, gamma
      count, and logarithmic data are not censored
------------------------------------------------------------------------------
version 0.9
------------------------------------------------------------------------------
28.4.02
gnlmix: corrected typo in negative binomial distribution
23.4.02
carma, chidden, elliptic, hidden, kalseries: give error if censored
	data supplied (thanks to Troels Ring)
22.4.02
elliptic: when two levels of nesting, calculate correctly first
	recursive fitted value in each cluster (was plotted correctly
	using iprofile) plus corresponding simplification of
	plot.iprofile (thanks to Troels Ring)
17.4.02 (all libraries)
gnlmix: corrected typo in inverse Gauss mixing distribution
print model methods: added option not to print correlations
15.3.02
restovec, tcctomat, tvctomat: added optional description slot to response,
	tccov, and tvcov objects
13.3.02
glmm: convert repeated object to dataframe if supplied as data
tcctomat, tvctomat: corrected to detect contrast options when dataframe=F
12.3.02
tvctomat: corrected problem for list of factors when dataframe=F
	(thanks to Patrick Lindsey)
finterp.default: give error if members of dataframe named using $
	(thanks to Christof Bigler)
28.2.02
chidden, hidden: added check for correct number of initial estimates
	when list of functions supplied (thanks to Patrick Lindsey)
22.2.02
corgram: added option for PACF
19.2.02
fmr: modified some discrete distributions to avoid overflow with large counts
17.2.02
elliptic: added as.double for y in call to C code because of change in
	  read.table
12.2.02
finterp: give error if offset used in W&R formula
31.1.02
%^%: power of a matrix
elliptic: corrected problem when common parameters in mean and
	variance functions
20.1.02
plot.repeated: added selection of profiles to be plotted by using ccov
14.1.02
gar: added absolute value arch (names of others changed: additive -> square,
	multiplicative -> exponential)
     volatility method for extracting values of nonconstant dispersion
	parameter
Makefiles: removed . for SHLIB_EXT for R1.4.0
dist.c, kcountb.c, romberg.c, stable.c: changed malloc to R_alloc
10.1.02
(dhpqr)ggamma, fmr, gausscop, gnlmix, gnlmm, gnlr, gnlr3, hgamma:
	changed argument of (dpqr)gamma for compatibility with R1.4.0
modified help to work with R1.4.0
18.12.01
contr.mean: provides correct labels on mean constraints (corrects contr.sum)
4.12.01
chidden, hidden: corrected printing out family parameter with AR when
	there is not one
28.11.01
qstable: corrected bug when tail<1 and skew=1 (thanks to Alec Stephenson)
23.11.01
corgram: handles NAs in the series
19.11.01
cprocess: fixed error in checking for list of events (thanks to Troels Ring)
stablereg: changed alpha to allow parameter to be (0,2) instead of (1,2)
18.11.01
chidden: added time-discretized Poisson process
17.11.01
chidden, hidden: added Student t distribution
		 changed Cauchy shape parameter to scale instead of scale^2
15.11.01
gar: added Student t distribution
     added ARCH models
elliptic: when AR, take log of determinant returned by dpodi
	(thanks to Gabrielle Kelly)
13.11.01
elliptic: when series of independent observations, calculate covariance
	determinant as sum of log variances instead of log of product
	(thanks to Gabrielle Kelly)
8.11.01
cmcre: corrected problems when a covariate is used (thanks to Anthony Gichangi)
6.11.01
print.response: do not print mean if nominal (but not binary) or ordinal
	(thanks to Patrick Lindsey)
25.10.01
hidden.r: corrected check for fixed zeros in transition matrix
	  relaxed check for rows of transition matrix summing to one
chidden.r: relaxed check for rows of transition matrix summing to zero
	(all thanks to Patrick Lindsey)
24.10.01
restovec: weights can be logical
14.10.01
gar: fixed output printing when shape is a function of location parameter
	use dnbinom function
	changed negative binomial shape parameter to be same as in gnlr
10.10.01
carma, chidden, gar, hidden, kalcount, kalseries: check for two levels of
	nesting when serial dependence fitted
9.10.01
kalseries: corrected error when torder used with tvcov
8.10.01
hidden, chidden: added observed AR(1)
gnlr, gnlmm, gnlmix: changed parametrization of the shape parameter for the
	beta distribution (thanks to Goran Arnoldsson)
binnest: duplicate variables in Fortran call
model functions using envir: check that response specified is one in
	envir when only one present
3.10.01 
plevy, qlevy: use pnorm and qnorm instead of integrating
26.9.01
elliptic: added second form of asymmetric multivariate Laplace
	distribution with constant asymmetry parameter
25.9.01
elliptic: added asymmetric multivariate Laplace distribution
24.9.01
carma.r: removed unnecessary check that envir is a repeated object
	(thanks to Troels Ring)
11.9.01
fit.dist: added checks that grouped frequency data are supplied
10.9.01
kalsurv: corrected output errors when environment is supplied
gar: use log option in dbinom, dpois
kalcount: set first recursive prediction in series to marginal prediction
6.9.01
gar: added loglog link for binomial data (corrected cloglog which was,
	in fact, loglog)
20.8.01
gnlmix: set undefined sh3 to NULL for one parameter distributions
1.8.01
chidden, gar, gnlr3, hidden: added skew Laplace distribution
27.7.01
corgram: improved presentation of correlogram
25.7.01
d,h,p,q,rskewlaplace: probability functions for the skew Laplace distribution
24.7.01
autointensity.r: plots autointensity function of a point process
12.7.01
plot.repeated: fixed error of unknown type when plotted time-varying covariate
	(thanks to Patrick Lindsey)
carma: clearer error message when incorrect environment supplied
10.7.01
carma: will handle data objects with (one of) multivariate responses
chidden, hidden: handle Jacobian correctly with (one of) multivariate
	responses
6.7.01
cprocess.r: recognizes data objects for events and not just for times
5.7.01
f2c.h included in library for toms614.c (missing in R1.3.0)
27.6.01
iprofile, mprofile: corrected links to other libraries for html help
plot.cum.pergram: corrected confidence interval
pergram: changed calculation of length for odd-lengthed series
22.6.01
gar.r: check that times are supplied and, if not, create if possible
19.6.01
fmr.r, gnlmm.r, gnlr.r, gnlr3.r: linear can be ~1 if mu not supplied
14.6.01
marg.hom.r: modified to handle factor variables with non-numeric levels
8.6.01
ordglm.r: corrected fitted values when weighted observations
	(thanks to Troels Ring)
31.5.01
elliptic: changed check on initial variance function estimates
16.5.01
print.response, print.tvcov, print.repeated: added option to print
	range of numbers of observations per individual instead of
	vector of numbers (thanks to Markus Jantti)
dmultpois, etc: added additional check on dependence parameter
9.5.01
gar.r: corrected printout for generalized gamma parameter
26.4.01
changed F and T to FALSE and TRUE throughout
read.rep: removed col.names option because of changes in read.table
22.4.01
glmm: corrected typo when dataframe used with no offset
20.4.01
finterp: detects functions when given as arguments of other functions
19.4.01
finterp: formulae can be written on several lines, as several
	instructions (e.g. to assign temporary variables)
11.4.01
dburr, pburr, qburr, hidden, chidden, gnlr3: changed parametrization
	of Burr distribution (thanks to Patrick Lindsey)
28.3.01
chidden, hidden: corrected vector length problem in check for ordered
	intercepts in ordinal models (thanks to Niko Speybroeck)
several p and d functions: changed check to y>=0 (thanks to Patrick Lindsey)
22.3.01
glmm: works again with weights and/or offset (thanks to Luc Duchateau)
gnlmix: changed to log dispersion for mixing distribution
21.3.01
int.c: corrected memory allocation problem
GExtVal: changed man page to agree with functions
	 (both thanks to Patrick Lindsey)
20.3.01
use log option in d and p functions for h functions
14.3.01
chidden, hidden: added further checks on ordering of intercepts for
	ordinal data
13.3.01
gnlmix: changed dispersion parameter for normal mixing distribution
	from standard deviation to variance
delta: returns a vector instead of a matrix if only for one variable
11.3.01
gnlmix: correction to work with unbalanced data
9.3.01
gar, gnlr3: added power variance function Poisson distribution
8.3.01
covariates: added expand option
dpvfpois, ppvfpois, qpvfpois, rpvfpois: functions for the
	overdispersed power variance function Poisson distribution
kalcount: corrected for power variance function
7.3.01
plot.response: corrected indexing problem
1.3.01
kalcount, kalseries, kalsurv: removed constraints on family parameter
27.2.01
chidden, fmr, gar, gausscop, gnlmix, gnlmm, gnlr, gnlr3, hidden,
	kalseries, kalsurv, nlr: relaxed type checks for continuous
	and duration data
26.2.01
kalcount, kalseries, kalsurv: added two-parameter power variance family
	mixture including gamma and inverse Gauss mixtures (family) for
	serial dependence
23.2.01
response: if response is univariate, returns a vector instead of a matrix
covariates: if only one covariate requested, return as a vector
chidden, hidden: improved checks for ordered intercepts with ordinal response
		 improved calculation of ordinal probabilities
22.2.01
plot(iprofile()): works for models from kalsurv (thanks to Jacob Bowers)
19.2.01
chidden.r, hidden.r: corrected error in calculating individual profiles when
		tvmu used (thanks to Jacob Bowers)
	ordinal data can be used with multinomial option (thanks to
		Patrick Lindsey)
	work with ordinal data with a list of different formulae
		(thanks to Niko Speybroeck)
31.1.01
glmm.r: works if response is a one-column matrix instead of a vector
	(thanks to Luc Duchateau)
restovec: corrected manual so that arguments section appears (thanks
	to Patrick Lindsey)
30.1.01
finterp, fnenvir: further correction to handle decimal numbers
	(including scientific notation) correctly
finterp: replaced gsub by all.vars
25.1.01
name of response can be a character string when environment is supplied
	(thanks to Patrick Lindsey)
hidden, chidden: added description of pintercept to man page
delta: works properly when name is supplied
plot functions: use 6 different line types instead of 4
gausscop: corrected mean for gamma margin
	  check that only one initial estimate when no dispersion function
18.1.01
transform.response: works when units is NULL
hidden, chidden: reversed order of categories for proportional odds
	and continuation ratio
	replaced dqrcf with dqrsl
nordr: bug fix to work with data objects
8.1.01
cutil.c, romberg.c, toms614.c: changed include for call_R
	(thanks to Dennis Murphy)
7.1.01
model fitting functions check for correct type of response
dftorep, read.rep: modified to handle new "types" of responses
dftorep: now handles two column binomial response
ehr.r: rewrote to conform to other functions
4.1.01
restovec: option to add responses to an old response object and
	  types of responses returned in list instead of as a class
resptype: new method to return types of response variable(s)
finterp.repeated: check that binomial and censored responses are not
	used as covariates
21.12.00
gar.r: corrected error in printing three-parameter distributions
19.12.00
finterp, fnenvir: methods for dataframes
gnlm functions: environment can be a dataframe
18.12.00
changed check for existence of response when environment supplied
bnlr, fmr, gnlmm, gnlr, gnlr3, nordr: fixed calculation of n for null function
17.12.00
various changes because of new R1.2.0 handling of formulae
finterp: check for + or - before ( (change to R1.2.0)
elliptic: removed check on tvcov, so can accept times and individuals
15.12.00
bnlr, fmr, gnlr, gnlr3, gnlmm, nordr: nonlinear formula need not
	contain covariates
changes to cutil.c and coxre for R1.2.0
14.12.00
restovec, tcctomat, tvctomat: added slot for units of measurement
dftorep, read.rep: id can be a factor variable
carma, elliptic, gar, gausscop, kalcount, kalseries, kalsurv: test if name of
	response exists when envir contains multivariate responses
stable: corrected man pages
3.12.00
gnlmix: generalized nonlinear models with one random parameter having
	arbitrary mixing distribution
30.11.00
fitdist.r: calculate log probabilities and likelihoods to avoid underflow
	(thanks to Goran Brostrom)
int2: vectorized two-dimensional Romberg integration
28.11.00
added d, p, q, and r functions for Consul generalized Poisson distribution
added PACKAGE option to .C and .Fortran calls
bnlr.r: added cloglog and loglog links
changed class of gnlr-type functions from gnlr to gnlm
gar.c: corrected calculation of censored Burr, Pareto, power
	exponential distributions (thanks to Patrick Lindsey)
24.11.00
q and r functions: improved calculations and checks
23.11.00
bnlr, fmr, gnlmm, gnlr, gnlr3, nordr: nonlinear formulae can have a
	linear part
qgweibull, qggamma, qgextval: corrected arguments to functions and
	docs (thanks to Patrick Lindsey)
21.11.00
fmobj: find objects referred to in a formula
elliptic.r, fmr.r, gar.r, gausscop.r, gnlmm.r, gnlr.r, gnlr3.r: models with
	common parameters in several regression functions can be specified
	using formulae, not just functions
elliptic.r, gar.r, gausscop.r, gnlmm.r, gnlr.r: models with shape as a
	function of location can be specified using formulae, not just
	functions
20.11.00
finterp.r: formulae can have common parameters and/or depend on functions
16.11.00
hidden, chidden: added recursive predicted values
added q and r functions for distributions in rmutil
14.11.00
kalseries: corrected error in inverse Gaussian distribution (thanks to
	Patrick Lindsey)
bnlr.r: added stable and mixture links
gnlr, gnlmm, gar: added beta and simplex distributions
rmutil: added psimplex and dsimplex
9.11.00
improved checking for multivariate response and choosing one response
	when several present in a data object
6.11.00
fmr.r, printrm.r: corrected so that works with common parameters
	(thanks to Laura Thompson)
29.10.00
gnlr3.r: corrected typo in normal and inverse Gauss distributions
19.10.00
gausscop: multivariate Gaussian copula with arbitrary marginals
elliptic.r: several typing errors corrected
------------------------------------------------------------------------------
version 0.8
------------------------------------------------------------------------------
17.10.00
carma, elliptic, kalseries: handles NULL delta correctly with
	multivariate response in repeated object
restovec: gives names correctly to multivariate 3-dim array
covariates.repeated: calculates number of observations correctly when
	multivariate response
15.10.00
qstable: corrected bug due to change in uniroot (thanks to Gabrielle Kelly)
dstable, pstable, qstable, rstable: added checks that parameter values
	are correct
13.9.00 (growth, repeated, rmutil)
restovec: check for NAs in ordinal responses (thanks to Patrick Lindsey)
elliptic, kalseries: check that torder is not larger than number of
	time points (thanks to Patrick Lindsey)
elliptic: corrected undefined n (thanks to Patrick Lindsey)
12.9.00
kalseries: constant shape parameter for (log) logistic, Cauchy, and Laplace
	distributions previously had square root transform
kalseries: added inverse Gauss distribution
7.9.00
restovec: corrected error (change in R) when censor is all ones
	(thanks to Troels Ring)
17.8.00
removed provide()
rmutil: removed det()
14.8.00
rmna: corrected typo in man page
17.7.00
nordr.r: corrected minor bugs for weights and data object handling
	(thanks to Patrick Lindsey)
5.7.00
as.data.frame.x: added options from default (thanks to Patrick Lindsey)
rmna: removes NAs in weights (thanks to Patrick Lindsey)
restovec: handle correctly option, times=T (thanks to Patrick Lindsey)
covariates.repeated: handle correctly non-repeated observations
	(thanks to Patrick Lindsey)
21.6.00
plotrm.r: plot.residuals corrected so ccov works (thanks to Patrick Lindsey)
14.6.00
carma.r: correction for ccov as one-column matrix (thanks to Patrick Lindsey)
7.6.00
elliptic.f: fixed crash with more than one covariate in tvcov (thanks
	to Bruno Genicot)
1.6.00
elliptic.r: corrected check to allow varfn="identity" or "square"
30.5.00
bnlr.r: binomial regression with various links
22.5.00
fnenvir.r: can handle functions without parameters (thanks to Troels Ring)
11.5.00
fit.dist: corrected exact fit for negative binomial and added default
	options for main and xlab
6.4.00
runge.kutta, lin.diff.eqn: functions to solve differential equations
5.4.00
gar.r: handles censored data correctly when a data object contains
	more than one response
29.3.00
runge.kutta.r: solution of differential equations
20.3.00
nlr: corrected undefined mu1
17.3.00
print.response: check for NAs in times
15.3.00
glmm: obtain nest vector from dataframe if supplied
14.3.00
nordr, ordglm: clearer error message if the response is not a numeric
	vector with integral values starting at 0 (thanks to Troels Ring)
15.2.00
ordglm: corrected bug when more than three categories
12.2.00 (repeated, event)
kalcount, kalseries, kalsurv: autoregression with frailty dependence
9.2.00
kcountb.c, kserieb.c, ksurvb.c, ksurvg.c: changed -log(1-pfn()) to
	-pfn(,,,0,1) and removed inthaz.c
8.2.00
all libraries: corrected C code for R0.99
kalcount: corrected error in recursive predicted values for gamma intensity
1.2.00
restovec: corrected handling of weights when response is a list
kalsurv.r: corrected plotting of profiles for logged distributions
cutil.c: changed Fortran.h to Rconfig.h and moved to rmutil
cgamma.c: replaced by cutil.c
inthaz.c: changed finite() to R_FINITE()
27.1.00
gar: three-parameter distributions work with constant dispersion parameter
kalcount, kalseries, kalsurv: if mu contains time-varying covariates,
	initial estimates must be in ptvc
24.1.00
finterp, fnenvir: changed name of .fn to avoid conflicts
most model functions: check that supplied function was not produced by
	finterp already
20.1.00
as.data.frame: puts binomial and censored responses as two-column matrices
gnlr, fmr, gnlmm, gar: binary response need only be one column for binomial
17.1.00
finterp, fnenvir: handle decimal numbers correctly
most model functions: print out regression function correctly when
	envir supplied
16.1.00
gar: added Consul generalized Poisson distribution
transform: check for nonpositive and infinite values in Jacobian
carma, elliptic, gar, kalseries: sqrt transformation checks for zero
	response values
14.1.00
most model functions: check for NAs in data objects (because of lvna)
gnlr, gnlr3, fmr, gnlmm: possible to fit a model without optimizing
	any parameters in the regressions if they are functions
restovec, dftorep, read.rep: add additional validity checks for responses
times.default: replaces times.repeated
11.1.00
hidden, chidden: handle multinomial count data
nlr: modified to handle data objects correctly
most model functions: changed way of detecting multivariate responses
finterp: correct error for length of response when multivariate
10.1.00
gettvc: works correctly for first observation when ties=FALSE (thanks
	to Patrick Lindsey)
finterp: can find response variables of formula in repeated objects
for most model functions, one of multiple responses in a repeated
	data object can be selected for the model
9.1.00
restovec: handles multivariate matrices, arrays, and lists
dftorep: transform a dataframe to a repeated data object
read.rep: read a rectangular data set from a file and create a
	  repeated data object directly
7.1.00
logitord.f: reformatted to remove tabs for HP compilers (thanks to
	    Osman Buyukisk)
restovec: responses can have more than one type class
------------------------------------------------------------------------------
version 0.7
------------------------------------------------------------------------------
3.1.2000
residuals.elliptic: corrected error in calculation of raw residuals
31.12.99
objectrm.r: can select certain individuals with methods, covariates, delta,
	    nesting, times, weights
transform: handles NAs in response correctly
28.12.99
restovec: added name of response variable to list returned
objectrm.r: added as.data.frame and as.matrix methods for data objects
wr: works with my data objects
nordr: changed sign of coefficients for continuation ratio and
       adjacent categories models so comparable with proportional odds
27.12.99
finterp: with W&R notation, design matrix no longer returned as attribute
	 when ~1 and .envir supplied, returns a function yielding a
		vector of correct length
26.12.99
fit.dist: corrected exact fit of negative binomial
gnlr, gnlr3, fmr, gnlmm: improved speed
nordr, ordglm: ordinal categories numbered from 0 instead of 1
hidden, chidden: multinomial categories numbered from 0 instead of 1
		 handles ordinal models, thanks to Patrick Lindsey
gettvc: now handles NAs in response variable
25.12.99
improved documentation for methods to access data objects and functions
int: call C instead of Fortran for TOMS614
23.12.99
restovec: added additional checks that correct data are supplied
	  (thanks to Troels Ring)
mprofile.carma: corrected bug when no covariates
carma: corrected bug when delta is a scalar
carma, elliptic: added checks for incorrect formulae
21.12.99
hidden, chidden: improved printout and corrected error in checking
	number of parameters in lists of formulae
20.12.99
lvna: creates a repeated object leaving NAs in
hidden, chidden: interactions between time-constant and time-varying
	covariates allowed
17.12.99
hidden, chidden: improved printout
16.12.99
tvctomat: handles lists of factor variables correctly
restovec: value returned has class of type of response as well as "response"
	  added checks
hidden, chidden: can also use formulae if multinomial
12.12.99
hidden, chidden: can use formulae if not multinomial
7.12.99
cmcre: corrected memory leak
6.12.99
cmcre: continuous-time two-state Markov process with random effects
5.12.99
coxre: corrected several errors
1.12.99
stable: fixed plot arguments in help examples
29.11.99
finterp: fixed bug when multiple ('s or ^ before ( when detecting
	 function names
nobs: use method instead of direct access in all functions
      provide default method
covind: provide default method
25.11.99
collapse: changed name to capply because of conflict in nlme
23.11.99
profile: changed to mprofile because of conflict in R0.90
22.11.99
finterp: properly distinguishes unknown parameters from functions
finterp and fnenvir: when no variables found, changed stop to warning
nobs: corrected for independent observations when length is one
18.11.99
stablereg: corrected bug when some parameters are not optimized
	   check for NAs in the hessian
17.11.99
plot.repeated, plot.response: added special call for ordinal responses
	corrected plot for independent observations (thanks to Patrick Lindsey)
14.11.99
removed unneeded aliases in man pages
added aliases to plot.profile and plot.iprofile
11.11.99
added check for Inf (as well as NAs) in hessian to all functions using nlm
kalseries.r: added error message if times not available for Markov dependence
	     changed rep(1,n) to rep(1,nind) when mu function returns scalar
stable.r: moved call to C code into likelihood function for speed
int.r: limits can be specified as -Inf, Inf
4.11.99
kalcount.r, kalseries.r, kalsurv.r: with time-varying covariates in a
	function or formula, initial estimates can be in preg or ptvc
	and changed length(resp$response$y) to n for speed
31.10.99
gar.r: fixed undefined npt3 for autoregression parameter
finterp.r: fixed bug for : in W&R formulae
kalcount.r, kalseries.r, kalsurv.r: added error message when time-varying
	covariates
22.10.99
covind: changed so that it works with carma, elliptic, gar, hidden,
	kalcount, kalseries, and kalsurv objects
18.10.99
gar.r: corrected printing of parameter values for three-parameter distributions
gar.c: corrected calculation of lambda in three-parameter distributions
17.10.99
gar.r: corrected fitted values (due to changes on 12.10.99)
14.10.99
gar.r: corrected undefined variable, tm (due to changes on 12.10.99)
12.10.99
stable: stablereg with nonlinear regression replaces stableglm
finterp, fnenvir: check for factor variables instead of not being a
	numerical or logical vector
gar: allow autoregression parameter to depend on covariates
dist.c, kcountb.c, kserieb.c, ksurvb.c, ksurvg,c, stable.c: added
	#include "Rconfig.h"
4.10.99
ordglm.r: added deviance and corrected for zeros in table
nordr.r: corrected typo
potthoff.r: corrected erroneous calculation of standard errors (thanks
	to Tim Auton)	
1.10.99
finterp, fnenvir: fixed conflict of names by beginning all names with a dot
	(thanks to Patrick Lindsey)
elliptic.r: changed option and title from elliptic to power exponential
30.9.99
ordglm.r: generalized linear ordinal regression
------------------------------------------------------------------------------
version 0.6
------------------------------------------------------------------------------
21.9.99
pkpd.r: changed mu2.1o1cfp to ensure ked>0
20.9.99
resid.f: correction to work with MS-Windows
7.9.99
binnest.r, survkit.r: changed NULLs for Fortran to work with R0.65
6.9.99
ehr.r, kalsurv.r, fmr.r, gnlr.r, gnlr3.r, nlr.r, nordr.r, elliptic.r,
	gar.r, gnlmm.r, kalcount.r, kalseries.r: changed attributes to
	work with R0.65
finterp, fnenvir: variables can be logical as well as numeric
3.9.99
Makefiles:  moved $(FLIBS) to end of line
14.8.99
print.gnlr: corrected errors in printing fmr, gnlr3, and gnlmm output
fnenvir.tvcov: corrected error for undefined ex1a (-> ex2a)
Pareto, gnlmm, hstudent, kalcount, kalseries, pkpd, read.list, read.surv,
	tvctomat: corrected examples
18.7.99
hidden.r, chidden.r: corrected one error message
		     added printout of degrees of freedom
14.7.99
binnest.f: modified comments to compile with standard Fortran (thanks
	to Martin Maechler)
------------------------------------------------------------------------------
version 0.5
------------------------------------------------------------------------------
29.6.99
plot.response: remove NAs when calculating default ylim
28.6.99
gnlr.r, gnlr3.r, fmr.r, nordr.r, gnlmm.r: check if user gives a
	nonlinear formula in linear argument and correctly handle it
27.6.99
finterp: corrected error message when non-numeric vector supplied
restovec: corrected printing of total times when negative times present
added transform methods for response, tccov, and tvcov objects
24.6.99
gar.r: corrected error in printing shape functions
22.6.99
binnest.f: modified to compile with g77
8.6.99
binnest: binary random effects model with two levels of nesting
7.6.99
restovec: added an additional check for nest variable in lists
6.6.99
logitord.f: corrected bug in calculation of Hessian (and s.e.)
1.6.99
elliptic: added multivariate Student t distribution
11.5.99
finterp.r: functions allowed in W&R formulae
carma.r, elliptic.r, kalseries.r, kalcount.r, kalsurv.r: allow factor variables
finterp: can locate and use indices for individuals and nesting as
	factor covariates
10.5.99
tcctomat.r, tvctomat.r: allow factor variables
finterp, fnenvir: changed to check for factor variables
6.5.99
elliptic.r: allow variance to be a function of the mean function
4.5.99
gar.c: changed normal distribution shape parameter from sd to variance
3.5.99
profile and iprofile: fixed to plot correctly with nesting
1.5.99
tcctomat, tvctomat: allow dataframes
28.4.99
tvctomat: time-varying covariates can be factors
elliptic.r, gnlr.r, gnlr3.r, fmr.r, gnlmm.r, gar.r: location and shape
	functions can have common parameters
26.4.99
restovec: weights allowed for lists
finterp, fnenvir: can find the times when envir is a repeated object
gar.r: allow shape to be a function of the location function
23.4.99
gnlr.r, gnlr3.r, fmr.r, nordr.r, nlr.r, elliptic.r, gnlmm.r, gar.r,
	kalseries.r, kalcount.r, kalsurv.r, ehr.r: do not require
	envir if response has class, repeated
corrected bugs in restovec and plot.response (Lorenz Gygax)
22.4.99
generalized plot.residuals
tvctomat: allow calculation of more than one interaction with
	time-constant covariates at a time
finterp and fnenvir: allow variables in environment to have same name
	as a function
21.4.99
correction of 18.1.99 by Brian Ripley wrong: dist.c put back in gnlm
20.4.99
ksurvb.c: corrected bug when time-varying covariates
elliptic.r: added option to ccov and tvcov to give covariate names
	when response has class, repeated
carma.r: added option to ccov to give covariate names when response
	has class, repeated
19.4.99
changed plot.profile to profile and plot.iprofile to iprofile
18.4.99
elliptic.r: added recursive predicted values when AR(1) and/or random effect
16.4.99
gnlr.r, gnlr3.r, fmr.r, nordr.r, gnlmm.r, ehr.r: changed order of
	parameters when function with linear part
15.9.99
ehr: corrected two errors when lambda function with linear part
nordr.r: n changed to nrows
13.4.99
carma.r: corrected predicted values when response is transformed
gar.r, kalseries.r: changed handling of transformed responses
------------------------------------------------------------------------------
version 0.4
------------------------------------------------------------------------------
12.4.99
added dependency on rmutil to DESCRIPTION
11.4.99
elliptic.f: corrected handling of dose for PKPD model when
	    time-varying covariates are present
6.4.99
elliptic.r, gnlmm.r, gar.r, kalseries.r, kalcount.r, kalsurv.r, ehr.r,
	 nordr.r, nlr.r: modified to use fnenvir
5.4.99
gnlr.r, gnlr3.r, fmr.r: modified to use fnenvir
4.4.99
fnenvir: checks functions for covariates and parameters and
	 modifies them to read from data objects
1.4.99
elliptic.r: modified to use model formulae with unknowns
finterp.r: added data objects as environment
tvctomat, tcctomat: can combine two data objects
31.3.99
gar.r: modified to use model formulae with unknowns
30.3.99
rmna: check if a covariate only has one value after NA removal
fixed examples docs so that they work
29.3.99
kalcount.r, kalsurv.r, ehr.r: modified to use model formulae with unknowns
28.3.99
gnlmm.r, kalseries: modified to use model formulae with unknowns
restovec: added coordinates to response class for spatial data
26.3.99
gnlr.r, gnlr3.r, fmr.r, nordr.r, nlr.r: modified to use model formulae
	with unknowns
24.3.99
changed language check to inherits formula in all functions
added methods for extracting elements from data objects
finterp.r: transforms model formulae with unknowns into functions
22.3.99
restovec: times no longer required for clustered data
	type attribute added
carma.r, elliptic.r, kalseries.r kalcount.r: check to see if times available
15.3.99
rmaov.r: wrote documentation
pkpd.r : added two new models and corrected one other
13.3.99
restovec: allow ties in times
23.2.99
gar.c: corrected Laplace cdf and allowed negative values
11.2.99
ehr: corrected for ties
kalsurv.r: prints out "birth process" when applicable instead of
	renewal process
logitord.r: removed DUP=F from Fortran call
8.2.99
km.r: fixed bug in plot.dist.km when several groups are plotted
	(Gareth Ridall)
7.2.99
improved handling of variable names in tcctomat, tvctomat, and
	functions calling them
rmaov.r: split-plot aov from Ralf Goertz
6.2.99
glmm.r: accepts transformed data if dataframe supplied
5.2.99
km.r: fixed bug for double lines with censored observations (Gareth Ridall)
ehr.r: modified handling of language
4.2.99
km.r: added print.km to remove attributes
restovec: accepts all response data, not just repeated measurements
tvctomat: added calculation of interactions
2.2.99
restovec: added adding several column matrices in lists with censoring
kalsurv.r: added delta option
1.2.99
glmm.r: binary response with binomial handled correctly
30.1.99
plot.iprofile.carma: corrected nind argument
restovec, carma, elliptic, kalcount, kalseries: added how to handle `times'
	  for clustered data to docs
28.1.99
bivbinom.r, marghom.r: minor corrections
rs.r: improved printout
26.1.99
readrm.r: corrected lty specification in plot.response
	  added option to plot points
24.1.99
gnlr.r, gnlr3.r, fmr.r, gnlmm.r: y can have classes, response or repeated
added DUP=F to all .C and .Fortran calls
pbirth.r: binomial data distributions
22.1.99
readrm.r: added ... for graphics options in plot.response and plot.repeated
21.1.99
rmna: added checks that ccov and tvcov have correct class
19.1.99
dist.c: changed static romberg to romberg2 and added static interp
carma.r, chidden.r, elliptic.r, gar.r, hidden.r, kalcount.r, kalseries.r,
	kalsurv.r: allow response to have class, repeated
restovec: allow delta to be a dataframe
18.1.99
corrections by Brian Ripley
gnlm: removed redundant dist.c
enclosed library.dynam in .First.lib
potthoff.r: added matrix dimension checks
util.r: removed orth function
potthoff.r: replaced orth by contr.poly
17.1.99
carma.r, chidden.r, elliptic.r, gar.r, hidden.r, kalcount.r, kalseries.r,
	kalsurv.r: copy response vector more efficiently
restovec: added total time for survival data
coxre.r: reorganized for efficiency, eliminating data.frame
cprocess.r: times can have class, response
16.1.99
gnlr.r, gnlr3.r, fmr.r, gnlmm.r: removed -delta/2 in right censoring
	calculation
dist.r, gnlr3.r, gar.c, hidden.f: changed parametrization of Burr to
	agree with kalsurv.r
elliptic.r: use var(y) to get initial estimate of variance
------------------------------------------------------------------------------
version 0.3
------------------------------------------------------------------------------
14.1.99
kalsurv.r: corrected printing of number of subjects and observations
2.1.99
cprocess.r: allow event counts with unequal times
added mode="double" to is.vector
29.12.98
corrected minor bugs in fmr
28.12.98
corrected abs bug for Laplace in kalman C functions
27.12.98
restovec: corrected binary totals when given as a vector
gar: added Levy, Pareto, generalized inverse Gauss, and
     power exponential distributions
hidden and chidden: added various overdispersed and continuous distributions
22.12.98
hidden and chidden: added filter calculation and plots
21.12.98
moved Student t from gnlr to gnlr3
renamed beta as Pareto in kalcount, kalseries, and kalsurv
corrected various minor errors in fmr and gnlr3
20.12.98
dist.r, gnlr.r, fmr.r: added gamma count and Pareto distributions
18.12.98
chidden: continuous-time hidden Markov chain models
7.12.98
dist.r, gnlr.r, fmr.r: added Levy distribution
removed .so from Makefiles and library.dynam 
6.12.98
util.r: added spectral decomposition to mexp
5.12.98
rmutil: added several p and d functions
gnlr3.r: added censored generalized inverse Gaussian and power
	 exponential distributions
2.12.98
int.r: vectorized Romberg integration
1.12.98
int.r: added option for Romberg integration
30.11.98
updated libraries with Brian Ripley's corrections
25.11.98
hidden:	allow values in the transition matrix to be fixed at 0 or 1
24.11.98
hidden: added independence model
23.11.98
inthaz.c: changed header include
bessel: removed function
gnlr3.r: changed to internal bessel function
14.11.98
hidden: added multinomial distribution
12.11.98
hidden.f: corrected Poisson and binomial calculations
5.11.98
carmasub.f and survkit.f: changes for compatibility with g77
------------------------------------------------------------------------------
version 0.2
------------------------------------------------------------------------------
2.11.98
ehr.r: corrected printing coefficients with linear and other parameters
1.11.98
km.r: corrected NaNs in log
29.10.98
carma.r: corrected printout of mean time
km.r: corrected ylab for cdf
26.10.98
rmna: handles NAs in time-constant covariates properly
carma.r and elliptic.r: accept ccov of class, tccov
cprocess.r: added plots from counts of events
19.10.98
changed to inherits() throughout
rationalized printing of gnlr, gnlr3, fmr, gnlmm and moved to rmutil
added delta option to carma and elliptic
18.10.98
carma.r and elliptic.r: added handling of delta when y has class, response
17.10.98
gar.r: added cloglog link
12.10.98
gnlmm.r: corrected handling of delta when y has class, response
11.10.98
replaced tapply() with collapse() in bivbinom, catmiss, glmm, gnlmm, ehr, coxre
10.10.98
ehr.r check for singular covariance matrix
      print names of variables for coefficients when language
8.10.98
kcountb.c: corrected dplogis call
gnlmm.r: corrected calls to ddp, dmp, ddb, and dmb
coxre.r: removed as.is=T in data.frame
corrected printing shape parameters when language used in gnlr, gnlr3,
	fmr, gnlmm
7.10.98
rs.r: put in check that all covariates are positive
gnlmm.r: set censor to F for binomial data
dist.c: changed ddp, dmp, ddb, and dmb to log and introduced weights
6.10.98
kseries.c: corrected error in serial update
kalseries.r: correcting printing error when there is an interaction
kalsurv: added serial update
inthaz.c: put back ihlogis (disappeared with nmath)
renamed wr.r as util.r
moved det and %**% from repeated and growth to rmutil/R/util.r
5.10.98
corrected check in carma, elliptic, gar, and kalseries for nonpositive
	transformed values
4.10.98
glmm.r: corrected two errors
1.10.98
extended residual plots to all of class recursive
kalcount, kalseries, kalsurv: return mean profiles in z$pred
plot.profile: accepts z$pred as well as a mean function
nbkal.r: corrections
corrected and updated a lot of docs
30.9.98
moved kalsurv to event library
renamed rmtools as rmutil
inthaz.c: corrected error from change to nmath
29.8.98
kalsurv.r: added recursive fitted values
kalseries.r: added recursive fitted values
updated plot.residuals for recursive class
27.9.98
corrected docs for plot.profile, plot.iprofile
added covind.default
plot.iprofile: corrected default settings
24.9.98
gettvc.r: allow an option for ties
bessel.r: only calculate one kind of function
20.9.98
gettvc.r: allow NAs in time-varying covariate
	  corrected for ties between response and covariate
tvctomat: allow tvcov to already have class, "tvcov"
added as.double in all Fortran and C calls
18.9.98
plotrm.r: corrected bug in plot.iprofile due to new covind()
16.9.98
pkpd.r: corrected mu2.0o2c and added mu2.0o2cfp
15.9.98
replaced Bessel, Gauss-Hermite, and integration routines
bessel.r: added docs
14.9.98
moved wr to rmtools and added docs
added covind function to rmtools
12.9.98
kalserie.r: added delta option
tcctomat.Rd: corrected alias
created new library, rmtools
11.9.98
dist.r: added beta binomial
dist.c: simplified calls to overdispersion functions
autocor.r: corrected pergram
kalserie.r: corrected error in printing parameters with torder>0
kserieb.c: corrected error when mu function used
10.9.98
readlist.r: corrected binomial totals for lists in restovec
fmr.r: removed unnecessary code
gar.r: added overdispersed binomial data
dist.r: allow dispersion to be a scalar when mean is a vector
	created documentation for p and d functions
9.9.98
nordr.r: corrected weights for adjacent categories model
	 test for p>1 in proportional odds
gar.r: added checks on times and mu arguments
       added binomial data
corrected docs for elliptic, gar, kalcount, kalseries, kalsurv, nbkal
	  for z$index
clarified docs for rmna, restovec, tcctomat, and tvctomat
8.9.98
removed backslash at end of Makefiles for event, gnlm, growth
moved integer declarations to the beginning in carmasub.f,
      elliptic.f, gettvc.f, survkit.f so that g77 should work
5.9.98
gar.r Corrected predictions for transformed responses
------------------------------------------------------------------------------
version 0.1
