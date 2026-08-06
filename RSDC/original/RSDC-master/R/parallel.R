# Internal parallel apply used by the multi-start search, the parametric
# bootstrap, and the correlation-band draws.

#' Apply FUN over X, sequentially or in parallel
#'
#' Work-horse behind the \code{cores} options. Uses \code{lapply()} when
#' \code{cores = 1}, forked workers (\code{parallel::mclapply}) on Unix-alikes,
#' and a PSOCK cluster on Windows. Tasks must be independent and deterministic
#' given their inputs (all call sites warm-start local optimizers or evaluate
#' pure functions), so the result is identical for any number of cores.
#' @param X A vector/list to iterate over.
#' @param FUN Function of one element of \code{X}.
#' @param cores Integer >= 1; capped at \code{parallel::detectCores()}.
#' @return A list, as \code{lapply(X, FUN)}.
#' @importFrom parallel mclapply detectCores makeCluster stopCluster clusterEvalQ parLapply
#' @noRd
.rsdc_lapply <- function(X, FUN, cores = 1L) {
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L)
    stop("cores must be a single integer >= 1.")
  max_cores <- parallel::detectCores()
  if (is.finite(max_cores) && cores > max_cores) {
    warning("cores = ", cores, " exceeds the ", max_cores,
            " available cores; using ", max_cores, ".")
    cores <- max_cores
  }
  if (cores == 1L || length(X) < 2L) return(lapply(X, FUN))
  if (.Platform$OS.type == "unix") {
    parallel::mclapply(X, FUN, mc.cores = cores)
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(RSDC))
    parallel::parLapply(cl, X, FUN)
  }
}
