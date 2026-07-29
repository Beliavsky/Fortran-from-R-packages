#' Named per-level richness configurations for the comparative experiment
#'
#' Convenience builder for the \code{level_spec} argument of
#' \code{\link{fit_ou_nested}}. Each level has four richness switches:
#' \code{cubic} (cubic drift), \code{sv} (stochastic volatility),
#' \code{student_t} (Student-t innovations) and \code{hierarchy} (cross-sector
#' partial pooling). The switches act only in the 2-level mode.
#'
#' The comparative experiment of the design contrasts these configurations by
#' parameter recovery and out-of-sample / leave-cluster-out predictive
#' performance, reporting which one the data support (no configuration is claimed
#' uniquely correct):
#' \describe{
#'   \item{\code{"canonical"} (alias \code{"n1_full_n2_lean"})}{Level 1 full
#'     (cubic + SV + Student-t + hierarchy), Level 2 lean (linear Gaussian OU with
#'     a hierarchical mean). The default.}
#'   \item{\code{"both_full"}}{Both levels full: the latent production price also
#'     gets cubic restoring, stochastic volatility and Student-t innovations.}
#'   \item{\code{"both_lean"}}{Both levels lean: linear Gaussian OU at each level,
#'     hierarchy retained.}
#'   \item{\code{"n1_lean"} (alias \code{"n1_lean_n2_lean"})}{Level 1 lean, Level 2
#'     lean.}
#' }
#' The fifth experiment arm, the single-level model, is obtained directly with
#' \code{fit_ou_nested(..., n_levels = 1)} (it takes no \code{level_spec}).
#'
#' @param config Character. One of \code{"canonical"}, \code{"both_full"},
#'   \code{"both_lean"}, \code{"n1_lean"} (and the explicit aliases above).
#'
#' @return A \code{level_spec} list with \code{level1} and \code{level2} entries,
#'   each a list of the four logical switches, ready to pass to
#'   \code{\link{fit_ou_nested}}.
#'
#' @examples
#' ou_level_spec("both_full")
#' ou_level_spec("both_lean")
#'
#' @seealso \code{\link{fit_ou_nested}}.
#' @export
ou_level_spec <- function(config = c("canonical", "n1_full_n2_lean",
                                     "both_full", "both_lean",
                                     "n1_lean", "n1_lean_n2_lean")) {
  config <- match.arg(config)
  full <- list(cubic = TRUE,  sv = TRUE,  student_t = TRUE,  hierarchy = TRUE)
  lean <- list(cubic = FALSE, sv = FALSE, student_t = FALSE, hierarchy = TRUE)
  switch(config,
    canonical        = ,
    n1_full_n2_lean  = list(level1 = full, level2 = lean),
    both_full        = list(level1 = full, level2 = full),
    both_lean        = list(level1 = lean, level2 = lean),
    n1_lean          = ,
    n1_lean_n2_lean  = list(level1 = lean, level2 = lean)
  )
}
