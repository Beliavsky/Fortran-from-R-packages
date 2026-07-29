# smoots 1.1.4
- The NOTE "Specified C++11: please drop specification unless essential"
  has been fixed.
- a package alias has been added appropriately in the package overview
  help file.
- the plot functionality for objects of class "smoots" is no longer 
  only interactive but different plot types can now also be selected
  via a function argument.
- the long-time deprecated argument "msg" in "bootCast()", "modelCast()",
  and "rollCast()" has been removed.

# smoots 1.1.3
- a bug has been fixed that occurred when selecting p = 0 or q = 0
  for the bootstrap approach
- a bug has been fixed where using the deprecated argument msg in the function
  rollCast() caused an error

# smoots 1.1.2
- a NAMESPACE issue (stated as a NOTE in CRAN checks) regarding the package
  'progress' has been fixed.
- an ERROR that occurred when running the test files on aarch64-apple-darwin20
  (64-bit), r-patched-solaris-x86 and x86_64-pc-linux-gnu (64-bit) has been 
  fixed.

# smoots 1.1.1

- a minor bug in the automatically created plot title when applying 
  plot.smoots() to an object returned by dsmooth() was fixed.
- a notable performance increase was implemented for the algorithms
  considered in the functions msmooth(), tsmooth() and dsmooth().
- RcppArmadillo was added to Suggests to address a note on 
  r-devel-linux-x86_64-fedora-clang. 
- a few typos in the documentation were fixed.
- an error in the description of the element 'ws' within lists returned by 
  msmooth(), tsmooth(), dsmooth() and gsmooth() was fixed in the documentation.
- an acknowledgment of the support by the German DFG project GZ-FE-1500-2-1 
  was added to the DESCRIPTION file.
- [[Rcpp::export]] was removed from C++ functions only used by other C++ 
  functions.
- minor R code simplifications were implemented.
- the package now makes use of the future and future.apply packages to improve 
  the performance of the bootstraps in the functions bootCast(), modelCast() 
  and rollCast().
- due to the implementation of parallel code, the argument msg of the functions 
  bootCast(), modelCast() and rollCast() has been deprecated; instead they have 
  two new arguments pb and cores.

# smoots 1.1.0

- the README file has been adjusted.
- new functions were added for graphical testing of linearity assumptions and 
  for forecasting.
- a rescaling function was added to simplify the transfer of the obtained 
  derivative estimates on the time interval [0, 1] to the actual observation  
  time points.
- the package now uses compiled C++ code within selected functions and needs to
  import the Rcpp package for this purpose.
- due to the implemented C++ code the performance of the functions of version
  1.0.0 has been slightly improved for this version.
- the documentation has been adjusted with respect to the newly introduced 
  functions.
- typos and false or missing information in the documentation have been 
  corrected. 
- the documentation has generally been edited with a stronger focus on text 
  formatting so that actual mathematical formulae are visible in the package's 
  PDF-manual
- the S3 method plot.smoots() has been reworked and allows for more arguments 
  of the standard plot function and thus for more flexibility.
- improved feedback on incorrect input.
- some default values for arguments of functions existing since package version 
  1.0.0 have been adjusted, e.g. p = v + 1 is now the default in gsmooth() as 
  opposed to p = 1 in previous versions.

# smoots 1.0.1

- the Imports field in the DESCRIPTION file was adjusted.
- minor performance improvements were made.
- the example in the README file was adjusted.
- minor changes in the documentation were made.
- minor changes in the print method were made.
