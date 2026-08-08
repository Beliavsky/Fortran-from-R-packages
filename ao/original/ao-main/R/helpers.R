#' Split estimate by target
#'
#' @description
#' This helper function splits the solution by target parameters (if provided),
#' which is used for the output.
#'
#' @author Siddhartha Chib
#'
#' @param estimate \[`numeric`\]\cr
#' A parameter vector.
#'
#' @inheritParams ao
#'
#' @return
#' A (named) \code{list}, a partition of `estimate` according to `npar`.

split_by_target <- function(estimate, target = NULL, npar) {
  if (is.null(target)) return(NULL)
  stopifnot(
    is.character(target), is.numeric(estimate),
    length(target) <= length(estimate), is.numeric(npar),
    sum(npar) == length(estimate)
  )
  out <- list()
  idx <- 1L
  for (i in seq_along(target)) {
    k <- as.integer(npar[i])
    out[[target[i]]] <- estimate[idx:(idx + k - 1L)]
    idx <- idx + k
  }
  out
}

#' Generate random partition
#'
#' @description
#' This helper function generates a random parameter partition, which is used
#' for the randomized AO procedure.
#'
#' @author Siddhartha Chib
#'
#' @param x \[`integer`\]\cr
#' The parameter indices.
#'
#' @param p \[`numeric(1)`\]\cr
#' The probability of generating a new block.
#'
#' @param min \[`integer(1)`\]\cr
#' The minimum number of blocks.
#'
#' @return
#' A \code{list}, a random partition of `x`.

generate_random_partition <- function(x, p, min) {
  if (min == length(x)) {
    return(as.list(x))
  }
  x <- sample(x, replace = FALSE)
  n <- length(x)
  y <- sample(0:1, n, replace = TRUE, prob = c(1 - p, p))
  y[1] <- 1
  ind <- which(y == 1)
  if (length(ind) < min) {
    ind <- sort(c(ind, sample(which(y == 0), size = min - length(ind))))
  }
  B <- length(ind)
  blocks <- vector("list", B)
  for (j in seq_len(B)) {
    s <- ind[j]
    e <- if (j < B) ind[j + 1] - 1 else n
    xj <- x[s:e]
    blocks[[j]] <- xj[order(xj)]
  }
  blocks
}

#' Merge optimization results
#'
#' @description
#' This helper function merges the results of multiple AO processes.
#'
#' @param results \[`list`\]\cr
#' A `list` of outputs from \code{\link[ao]{ao}}.
#'
#' @param processes \[`data.frame`\]\cr
#' A `data.frame` describing how the different processes were specified.
#'
#' @inheritParams ao
#'
#' @return
#' A \code{list}, see section "Output value" on the \code{\link[ao]{ao}} page.

merge_results <- function(
    results, minimize = TRUE, add_details = TRUE, processes = data.frame()
  ) {
  values <- vapply(results, `[[`, numeric(1), "value")
  stopping_reasons <- vapply(results, `[[`, character(1), "stopping_reason")
  optimal_process <- ifelse(
    isTRUE(minimize), which.min(values), which.max(values)
  )
  seconds_each <- vapply(results, `[[`, numeric(1), "seconds")
  has_es <- sapply(results, function(x) "estimate_split" %in% names(x)) |> any()
  if (isTRUE(add_details)) {
    details_list <- list()
    for (process in seq_along(results)) {
      details_list[[process]] <- cbind(
        process = process,
        results[[process]][["details"]]
      )
    }
    return(
      c(
        list(
          "estimate" = lapply(results, `[[`, "estimate")[[optimal_process]],
          "estimates" = lapply(results, `[[`, "estimate")
        ),
        if (isTRUE(has_es)) list(
          "estimate_split" = lapply(
            results, `[[`, "estimate_split"
          )[[optimal_process]]
        ),
        list(
          "value" = values[optimal_process],
          "values" = as.list(values),
          "details" = do.call("rbind", details_list),
          "seconds" = sum(seconds_each),
          "seconds_each" = as.list(seconds_each),
          "stopping_reason" = stopping_reasons[optimal_process],
          "stopping_reasons" = as.list(stopping_reasons),
          "processes" = processes
        )
      )
    )
  } else {
    return(
      c(
        list(
          "estimate" = lapply(results, `[[`, "estimate")[[optimal_process]]
        ),
        if (isTRUE(has_es)) list(
          "estimate_split" = lapply(
            results, `[[`, "estimate_split"
          )[[optimal_process]]
        ),
        list(
          "value" = values[optimal_process],
          "seconds" = sum(seconds_each),
          "stopping_reason" = stopping_reasons[optimal_process]
        )
      )
    )
  }
}

