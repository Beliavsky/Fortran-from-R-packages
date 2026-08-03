#' List available example scripts
#'
#' @description
#' Shows all example scripts included with the PortfolioTesteR package.
#' These examples demonstrate various strategy patterns and library functions.
#'
#' @return Character vector of example filenames
#' @export
#'
#' @examples
#' # See available examples
#' list_examples()
#'
#' # Run a specific example
#' # run_example("example_momentum_basic.R")
list_examples <- function() {
  example_dir <- system.file("examples", package = "PortfolioTesteR")
  if (example_dir == "") {
    stop("Examples directory not found. Please reinstall the package.")
  }
  examples <- list.files(example_dir, pattern = "\\.R$")
  if (length(examples) == 0) {
    message("No examples found.")
  } else {
    message("Available examples:")
    for (ex in examples) {
      message("  - ", ex)
    }
  }
  invisible(examples)
}

#' Run an Example Script
#'
#' @description
#' Executes an example script bundled in the package `inst/examples/` folder.
#'
#' @param example_name Character scalar with the example filename (e.g. `"basic.R"`).
#' @param echo Logical; print code as it runs (default `TRUE`).
#'
#' @return Invisibly returns `NULL`. Runs the example for its side effects.
#' @export
#'
#' @examplesIf interactive()
#' \donttest{
#' # Example (requires a real file under inst/examples):
#' # run_example("basic.R")
#' }
run_example <- function(example_name, echo = TRUE) {
  example_path <- system.file("examples", example_name, package = "PortfolioTesteR")

  if (!file.exists(example_path)) {
    stop("Example '", example_name, "' not found. Use list_examples() to see available examples.")
  }

  message("Running example: ", example_name)
  message("=" , strrep("=", nchar(example_name) + 17))

  source(example_path, echo = echo)
  invisible(NULL)
}
