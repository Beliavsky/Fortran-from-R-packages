#' @aliases ao-package
#' @inherit ao examples
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom checkmate assert_int
#' @importFrom checkmate assert_number
#' @importFrom checkmate test_class
#' @importFrom checkmate test_flag
#' @importFrom checkmate test_function
#' @importFrom checkmate test_list
#' @importFrom checkmate test_number
#' @importFrom cli cat_line
#' @importFrom cli cli_abort
#' @importFrom cli style_hyperlink
#' @importFrom future.apply future_lapply
#' @importFrom oeli match_arg
#' @importFrom oeli test_numeric_vector
#' @importFrom optimizeR Objective
#' @importFrom progressr progressor
#' @importFrom R6 R6Class
#' @importFrom stats rnorm
#' @importFrom utils packageVersion
## usethis namespace: end
NULL

#' @noRd

.onAttach <- function(lib, pkg) {
  doc_link <- "https://loelschlaeger.de/ao"
  issues_link <- "https://github.com/loelschlaeger/ao/issues"
  msg <- c(
    paste0("This is {ao} ", utils::packageVersion("ao")),
    ", happy alternating optimization!\n",
    "Documentation? ",
    cli::style_hyperlink(doc_link, doc_link), "\n",
    "Any issues? ",
    cli::style_hyperlink(issues_link, issues_link)
  )
  packageStartupMessage(msg)
  invisible()
}
