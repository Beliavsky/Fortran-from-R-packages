# HEADER ####################################################
# This is file spam/R/random.R.                             #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



spam_random <- function(nrow = 1L, ncol = nrow, density = 0.5, distribution = NULL,
                        digits = NULL, sym = FALSE, spd = FALSE, verbose = FALSE, ...) {

  if (!is.numeric(nrow) || length(nrow) != 1)
    stop("nrow must be a single numeric value.")
  if (!is.numeric(ncol) || length(ncol) != 1)
    stop("ncol must be a single numeric value.")
  if (!is.numeric(density) || length(density) != 1 || density < 0 || density > 1)
    stop("density must be a single numeric value, equal or between 0 and 1.")
  if (!(is.null(distribution) || (is.function(distribution) && any(names(formals(distribution)) == "n"))))
    stop("distribution must be a function which generates random deviates and must have an argument \"n\".")
  if (!(is.null(digits) || (is.numeric(digits) && length(digits) == 1 && digits >= 0)))
    stop("digits must be a single numeric value larger or equal to 0.")
  if (!is.logical(verbose) || length(verbose) != 1)
    stop("verbose must be a single locigal value.")
  if (!is.logical(sym) || length(sym) != 1)
    stop("sym must be a single locigal value.")
  if (!is.logical(spd) || length(spd) != 1)
    stop("spd must be a single locigal value.")
  if (spd && ((ncol != nrow) || density == 0))
    stop("if spd == TRUE, then it must hold that ncol == nrow and density > 0")

  if (density == 1)
    npercol <- rep.int(ncol, nrow)
  else if (density == 0)
    return(spam(0, nrow, ncol))
  else
    npercol <- as.integer(rbinom(n = nrow, size = ncol, prob = density))

  rowp <- c(1L, cumsum(npercol) + 1L)
  n <- sum(npercol)
  if (n == 0) {
    if (sym)
      if (is.null(distribution))
        return(diag.spam(nrow = nrow, ncol = ncol))
    else
      return(diag.spam(x = distribution(min(nrow, ncol), ...),
                       nrow = nrow, ncol = ncol))
    else
      return(spam(0, nrow = nrow, ncol = ncol))
  }

  coli <- integer(n)
  seqcol <- seq_len(ncol)
  for (i in seq_len(nrow)) {
    if (rowp[i] == rowp[i+1L])
      next
    coli[rowp[i]:(rowp[i+1L]-1L)] <- sort(sample(seqcol, npercol[i], replace = FALSE))
  }

  if (is.null(distribution))
    entries <- rep.int(1, n)
  else
    entries <- distribution(n, ...)
  if (!is.null(digits))
    entries <- round(entries, digits)

  rspam <- .newSpam(entries,
                    colindices = coli,
                    rowpointers = rowp,
                    dimension = c(nrow, ncol),
                    force64 = getOption("spam.force64"))

  if (spd) {
    sym <- TRUE
    diag.spam(rspam) <- diag.spam(rspam) + 1 }

  if (sym) {
    s_t <- t.spam(rspam)
    rspam <- s_t[lower.tri(s_t, diag = TRUE)] + rspam[upper.tri(rspam, diag = TRUE)]
  }

  if (spd) {
    st <- rspam
    st@entries <- abs(st@entries)
    st_dia <- diag.spam(st)
    st_rs <- rowSums(st)
    s_diag <- diag(rspam)
    ind <- 2*st_dia <= st_rs
    diag.spam(rspam)[ind] <- ifelse(s_diag >=0, st_rs+1, -st_rs-1)[ind]
  }

  if(verbose)
    message("Density is ", round(n/prod(rspam@dimension), 7), ", specified is ", density,
      " (nnz=",n,").")

  return(rspam)
}

