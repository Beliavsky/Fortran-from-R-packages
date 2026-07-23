# timsac 1.3.8-6

* Changed package maintainer.


# timsac 1.3.8-4

* Fixed usage of KIND in Fortran file for warning reported in CRAN package check results.


# timsac 1.3.8-3

* Specified all variable types explicitly instead of using implicit variable types.


# timsac 1.3.8-2

* Fixed an issue with the S3 method that was reported as NOTE in the pre-test.

  Defined method for S3 classes 'decomp' and 'specmx'.


# timsac 1.3.8-1

* Collected legacy Fortran code (DCMPLX, DREAL, DIMAG, DCABS, DCONJG).


# timsac 1.3.8

* Removed C wrapper functions and registered entry points for the routines accessed by the `.Fortran` interface to call Fortran subroutines.

* Fixed bugs in `decomp()` and `mulmer()`.

* Added a ‘`NEWS.md`’ file.


# timsac 1.3.7

* Removed arguments `year` and `month`, and changed the default value of `period` in `decomp()`.
 
  If attributes of the input time series are required, set the start time and frequency using `ts()` or `tsp()`.

* Modified the plotting method for missing values in `decomp()`.


# timsac 1.3.6

* Fixed a bug in `prdctr()` for memory access issues with ASan.

* Added ‘`src/init.c`‘.

* Changed the names of C routines to be different from the function names in R.


# timsac 1.3.5

* Fixed a bug in `bsubst()` for memory access issues with Valgrind.


# timsac 1.3.4

* Fixed Fortran code according to the warning message when using gfortran with `-Wall -pedantic`.


# timsac 1.3.3-2

* Changed explanation of value tdf and examples in `decomp()` help page.


# timsac 1.3.3

* Fixed `inst/doc/index.html` according to the error message from the W3C Markup Validator Service.


# timsac 1.3.2

* Fixed a bug in `mulmar()`.
  If max.order (upper limit of AR model order) is too large, output results will be incorrect. 
  In such a case, `max.order` is changed to the appropriate value. (Reported by Jiri Baum.)

* Fixed bugs in `armafit()`, `baysea()`, `bispec()`, `decomp()`, `mlomar()`, `mulmar()` and `perars()` for memory access issues on Valgrind. (Reported by Brian Ripley.)

* Fixed a bug in `simcon.c`. (Reported by Kurt Hornik and Brian Ripley.)


# timsac 1.3.1

* Corrected spelling error in `Airpollution`. (Reported by Vivek Rao.)


# timsac 1.3.0

* Fixed array boundary dimensions (1) and (*) in Fortran source code.

* Removed `armaimp()`, `lsar2()`, `ngsmth()`, `tsmooth()`, `tvvar()`, `tvar()` and `tvspc()` functions.
  These functions have been added to the **TSSS** package.

* Deleted a document `f77program.pdf` in the doc directory.


# timsac 1.2.8

* Added `lsar2()` function.

* Defined S3methods (print) for `blocar()`, `mlocar()`, `nonst()` and `tvvar()` classes.

* Changed argument name `ord` to `subset` in `perars()`.

* Removed argument `ar.order` in `tvspc()`.

* Added a document `f77program.pdf` to the doc directory.

* Fixed bugs in `autoarmafit()` and `armafit()`.
  Added an array CXY equivalent to CXX to subroutine `SC0GRH` in autarmf.f.
	
* Removed the default value for `delta` in `tvvar()` and `tvar()` functions.

* Fixed a bug in `decomp()`.
  If argument `log` is TRUE and there is a zero or negative number in input y, print an error message.

* Fixed a bug in `mulcor()`.
  In mulcor.c, the size of allocVector is not enough.


# timsac 1.2.7

* Fixed a bug in `mulmar()` that did not work for certain types of data. (Reported by William Pleasant.)


# timsac 1.2.6

* Added C wrapper functions for calling Fortran subroutines.

* Removed arguments `tmp.file` from `autoarmafit()`, `armafit()`, `exsar()`, `markov()`, `mulmar()` and `unimar()` functions.

* Changed the argument type of `smt` from `single` to `double` to use the C wrapper function in `ngsmth()`.

* Fixed a bug in `bsubst()`.
  In the case of mtype = 3, argument `nreg` (number of regressors) of `bsubst()` is automatically set.


# timsac 1.2.5

* Fixed a bug in `markov()`.
  Added new arguments `nh` to subroutines `C0GR` and `NSUBX1` in markovf.f.

  The following warnings were found in 00install.out (GNU Fortran (GCC) 4.7.0) :
  Warning: Rank mismatch in argument `nh` at (1) (rank-1 and scalar)

* Removed arguments `param` from FUNCTION UNIF and FUNCTION DBLEXP in ngsmthf.f.
  
  The following warnings were found in 00install.out (GNU Fortran (GCC) 4.7.0) :
  
		Warning: Rank mismatch in argument 'param' at (1) (scalar and rank-1)
  
		Warning: Rank mismatch in argument 'param' at (1) (scalar and rank-1)


# timsac 1.2.4

* Removed the value `spec` from `tvar()`.

* Added `tvspc()` function.

* Fixed a bug that could cause abnormal termination of routine processing when argument `trend.order` to `tvar()` was set to 1.
  Fixed the dimension of array A in subroutine `setcar` in tvarf.f.
  

# timsac 1.2.3

* Removed arguments `ncon`, `nman` and `inw` and added new arguments `control` and `manip` to `fpec()` and `mulnos()`.

* Removed argument `niv` from `mulfrf()`.

* Changed the allowed values of argument `initd` from {0,1,2} to {1,2,3} in `ngsmth()`.

* Deleted C wrapper functions.


# timsac 1.2.2

* Fixed a bug that could cause abnormal termination of routine processing when argument `seasonal.order` to `decomp()` was zero.
  (reported by Seisho Sato.)
  Fixed the dimension of working array D in subroutine `smoth3`.

* Set the default value for `outlier` to NULL in `tvar()`.
 
* Added overview help page `timsac-package`.


# timsac 1.2.1

* Changed argument `imiss` to `miss`.

* Added `Details` description to the help page for `decomp()`. 


# timsac 1.2.0

* Added `baysea()` function.

* Added a new argument `lag` to `canarm()`.


# timsac 1.1.4

* Added a set of macros (F77_NAME(name), F77_CALL(name)) in all C source files (arma.c, auspec.c, ..., xsarma.c).

* Fixed a bug in `fpec()`.
  `aperm(arcoef, c(2,3,1))` gives an error if z$ordermin = 1. (Reported by Seisho Sato.)


# timsac 1.1.3

* Fixed a bug in `blomar()`.
  Added new arguments `F1` and `F2` to subroutine `MNONSB` in blomarf.f.

* Fixed a bug in `mlomar()`.
  Added new arguments `X` and `U` to subroutine `MNONST` in mlomarf.f.


# timsac 1.1.2

* Added `tsmooth()` function.


# timsac 1.1.1

* Changed the type of argument `smt` from `double` to `single` in `ngsmth()`.

* Changed the function name from `arma()` to `armaimp()`.


# timsac 1.1.0

* Added `arma()`, `tvvar()`, `tvar()` and `ngsmth()` functions.


# timsac 1.0.1

* Fixed a bug that occurred when `decomp()` did not take into account both seasonal components and trading day effects.

	In the following example, this problem occurs:
  
		data(Blsallfood)
    
		decomp(Blsallfood, seasonal.order=0)   # trade=FALSE (default)
  
  (Reported by Makoto Kodama.)
