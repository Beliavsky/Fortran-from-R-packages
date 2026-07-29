.onAttach <- function(libname, pkgname){
  packageStartupMessage("********************************************************************************")
  packageStartupMessage("                        Welcome to the package 'smoots'!")
  packageStartupMessage("********************************************************************************")
  packageStartupMessage("")
  packageStartupMessage("Please report any possible errors and bugs to dominik.schulz@uni-paderborn.de.")
  packageStartupMessage("Thank you.")
  packageStartupMessage("")
  packageStartupMessage("********************************************************************************")
}

.onUnload <- function(libpath) {
  library.dynam.unload("smoots", libpath)
}
