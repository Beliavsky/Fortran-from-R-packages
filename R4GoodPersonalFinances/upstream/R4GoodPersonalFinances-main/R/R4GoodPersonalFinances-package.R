#' @docType package
#' @name R4GoodPersonalFinances-package
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

.pkg_env <- new.env(parent = emptyenv())
.pkg_env$memoised <- list()

.onLoad <- function(libname, pkgname) {
  set_cache(path = file.path(tempdir(), ".cache"))
}

.onAttach <- function(libname, pkgname) {

  pkg_version <- utils::packageVersion(pkgname)
  packageStartupMessage(
    glue::glue("Welcome to {pkgname} version {pkg_version}!")
  )
  
  packageStartupMessage(
    "Cite the package: citation('R4GoodPersonalFinances')"
  )

  packageStartupMessage(
    "Package documentation: https://r4goodacademy.github.io/R4GoodPersonalFinances/"
  )

  packageStartupMessage(
    "To learn more, visit: https://www.r4good.academy/"
  )

  packageStartupMessage(
    "... and Make Optimal Financial Decisions!"
  )

}

ignore_unused_imports <- function() {
   bsicons::bs_icon
   bslib::card
   shiny::a
}
