#' Locate the canonical Stan model file
#'
#' Resolves the path to \code{inst/stan/ou_nested.stan}, the single source of
#' truth for the model (single-level and 2-level nested in one file). Works both
#' for the installed package (via \code{system.file}) and during development
#' (source-tree fallbacks).
#'
#' @return Character path to the \code{.stan} file.
#'
#' @keywords internal
#' @noRd
.stan_file_path <- function() {
  p <- system.file("stan", "ou_nested.stan", package = "bayesianOU")
  if (nzchar(p) && file.exists(p)) {
    return(p)
  }

  # Development fallbacks (devtools::load_all, running tests from source, etc.)
  candidates <- c(
    file.path("inst", "stan", "ou_nested.stan"),
    file.path("..", "inst", "stan", "ou_nested.stan"),
    file.path("..", "..", "inst", "stan", "ou_nested.stan")
  )
  for (cand in candidates) {
    if (file.exists(cand)) {
      return(normalizePath(cand))
    }
  }

  stop(
    "Could not locate 'ou_nested.stan'. Reinstall the package or run from the ",
    "package root.",
    call. = FALSE
  )
}


#' Stan code for the unified nonlinear OU model
#'
#' Returns the complete Stan code for the unified Ornstein-Uhlenbeck model: a
#' single file (\code{inst/stan/ou_nested.stan}) covering the single-level mode
#' (\code{n_levels = 1}, cubic drift, stochastic volatility, Student-t
#' innovations, non-centered hierarchical priors and a TMG interaction) and the
#' 2-level nested mode (\code{n_levels = 2}, market price reverting to a latent
#' production price with its own G'-driven OU). The code is read from the
#' canonical file (single source of truth); this function does not embed a
#' duplicate copy.
#'
#' @return Character string containing the Stan model code.
#'
#' @examples
#' code <- ou_nested_stan_code()
#' cat(substr(code, 1, 300))
#'
#' @export
ou_nested_stan_code <- function() {
  paste(readLines(.stan_file_path(), warn = FALSE), collapse = "\n")
}
