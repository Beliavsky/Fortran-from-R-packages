setClass("FLXMCsparse",
         contains = "FLXM")

as.data.frame.simple_triplet_matrix <- function(x, ...) {
    as.data.frame.model.matrix(x, ...)
}
              
setMethod("FLXgetModelmatrix", signature(model = "FLXMCsparse"),
          function(model, data, formula, lhs=TRUE, ...) {

  formula <- RemoveGrouping(formula)
  if (length(grep("\\|", deparse(model@formula))))
      stop("no grouping variable allowed in the model")
  if(is.null(model@formula))
    model@formula <- formula

  model@fullformula <- update(terms(formula, data = data), model@formula)
  fullformula <- terms(model@fullformula, data = data)
  model@terms <- attr(fullformula, "terms")

  if (lhs) {
      env <- environment(fullformula)
      vars <- attr(fullformula, "variables")
      varnames <- vapply(vars, function(x)
          paste(deparse(x, backtick = FALSE), collapse = " "), " ")[-1L]
      variables <- eval(vars, data, env)
      resp <- attr(fullformula, "response")
      response <- variables[[resp]]
      model@y <- model@preproc.y(response)
      model@x <- matrix(nrow = nrow(model@y), ncol = 0)
  } else {
      model@x <- matrix(nrow = nrow(as.data.frame(data)), ncol = 0)
  }
  model
})
