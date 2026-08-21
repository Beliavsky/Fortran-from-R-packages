.onAttach <- function(libname, pkgname) {
  packageStartupMessage("The argument `scale` has been replaced with `phi` in this version.")
  packageStartupMessage("The top of page 77o in Wang Louis 2003 states phi is not a scale.")
  packageStartupMessage("Please start using the `phi` argument now.")
}
