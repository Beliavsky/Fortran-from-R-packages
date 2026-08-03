### print vector/matrix while limiting the number of rows/columns printed

printRowColLimits <- function(x,
                              max.rows=getOption("max.rows", 20),
                              max.cols=getOption("max.cols", 7),
                              ...  # other arguments to 'print.matrix'
                              ) {
   x1 <- x
   msg <- NULL
   if(is.null(dim(x))) {
      x1 <- matrix(x, nrow=1)
      colnames(x1) <- names(x)
      x <- x1
   }
   ## we have a matrix (higher-D arrays not supported)
   if(ncol(x) > max.cols) {
      x1 <- x[, seq(length=max.cols), drop=FALSE]
      msg <- paste(msg, "reached getOption(\"max.cols\") -- omitted",
                   ncol(x) - max.cols, "columns\n")
   }
   print(head(x1, max.rows), ...)
   if(nrow(x) > max.rows) {
      msg <- paste(msg, "reached getOption(\"max.rows\") -- omitted",
                   nrow(x) - max.rows, "rows\n")
   }
   cat(msg)
}
