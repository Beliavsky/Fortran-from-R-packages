# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Distribution Calculus                                                     ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculus               ----
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Distributional Calculus
#' @name calculus
#'
#' @description
#' Arithmetic operators and functions for probability distribution objects.
#' These methods define how standard operations like \code{+}, \code{-},
#' \code{*}, and \code{/} behave when applied to random variables, returning the
#' resulting distribution based on known properties of common distribution
#' families.
#'
#' @param x,e1,e2 objects of subclass `Distribution`.
#' @param na.rm logical. Should missing values be removed?
#' @param ... extra arguments.
#'
#' @return All calculations return `Distribution` objects (specifically, objects
#' of a class that is a subclass of `Distribution`), according to the property
#' at hand.
#'
#' @export
#'
#' @examples
#' # -----------------------------------------------------
#' # Distribution Calculus Example
#' # -----------------------------------------------------
#'
#' # Normal location - scale transformation
#' x <- Norm(2, 3)
#' y <- 3 * x + 1 # Norm(7, 9)
#'
#' # Addition of two independent Normal random variables
#' x1 <- Norm(1, 3)
#' x2 <- Norm(2, 4)
#' x3 <- x1 + x2 # Norm(3, 5)
NULL

#' @rdname calculus
setMethod("+", signature = c(e1 = "Norm", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1@mean + e2@mean,
                 sd = sqrt(e1@sd ^ 2 + e2@sd ^ 2))
          })

#' @rdname calculus
setMethod("+", signature = c(e1 = "numeric", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1 + e2@mean, sd = e2@sd)
          })

#' @rdname calculus
setMethod("+", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            e2 + e1
          })

#' @rdname calculus
setMethod("-", signature = c(e1 = "Norm", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1@mean - e2@mean,
                 sd = sqrt(e1@sd ^ 2 + e2@sd ^ 2))
          })

#' @rdname calculus
setMethod("-", signature = c(e1 = "numeric", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1 - e2@mean, sd = e2@sd)
          })

#' @rdname calculus
setMethod("-", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            Norm(mean = e1@mean - e2, sd = e1@sd)
          })

#' @rdname calculus
setMethod("*", signature = c(e1 = "numeric", e2 = "Norm"),
          function(e1, e2) {
            Norm(mean = e1 * e2@mean, sd = e1 * e2@sd)
          })

#' @rdname calculus
setMethod("*", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            e2 * e1
          })

#' @rdname calculus
setMethod("/", signature = c(e1 = "Norm", e2 = "numeric"),
          function(e1, e2) {
            (1 / e2) * e1
          })

#' @rdname calculus
setMethod("sum", signature = c(x = "Norm", na.rm = "logical"),
          function(x, ..., na.rm = FALSE) {
            d <- list(x, ...)
            m <- unlist(lapply(d, FUN = function(x) { x@mean }))
            s <- unlist(lapply(d, FUN = function(x) { x@sd }))
            Norm(mean = sum(m), sd = sqrt(sum(s ^ 2)))
          })

#' @rdname calculus
setMethod("exp", signature = c(x = "Norm"),
          function(x) {
            Lnorm(meanlog = x@mean, sdlog = x@sd)
          })
