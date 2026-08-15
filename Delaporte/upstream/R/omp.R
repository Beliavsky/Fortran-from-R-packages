# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# Function to get current active number of threads
getDelapThreads <- function() {
  .Call(gOMPT_C)
}

setDelapThreads <- function(n) {
  n <- as.integer(floor(n))
  if (n <= 0L) {
    stop("Number of threads must be > 0.")
  }
  ncpus <- get("DLPCPU", envir = DelaporteEnv)
  if (n > ncpus && !is.na(ncpus)) {
    message("Capping at system maximum of ", ncpus, ".")
    n <- ncpus
  }
  invisible(.Call(sOMPT_C, n))
}
