# HEADER ####################################################
# This is file spam/R/eigen.R.                              #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



setMode <- function(sMode, symmetric, silent = FALSE){

  if (symmetric) {
    modes <- c("LM", "SM", "LA", "SA") #,"BE")
    if (!(sMode %in% modes)) {
      sMode <- 'LM'
      if (!silent) { warning(paste("the control option \"mode\" is set to", sMode), call. = FALSE) } }
  } else {
    modes <- c("LR", "SR", "LI", "SI", "SM", "LM")
    if (!(sMode %in% modes)) {
      sMode <- 'LR'
      if(!silent) { warning(paste("the control option \"mode\" is set to", sMode), call. = FALSE) } }
  }

  if (sMode == 'LM') { imode <- as.integer(1) }
  if (sMode == 'SM') { imode <- as.integer(2) }
  if (sMode == 'LR') { imode <- as.integer(3) }
  if (sMode == 'SR') { imode <- as.integer(4) }
  if (sMode == 'LI') { imode <- as.integer(5) }
  if (sMode == 'SI') { imode <- as.integer(6) }
  if (sMode == 'LA') { imode <- as.integer(7) }
  if (sMode == 'SA') { imode <- as.integer(8) }
  if (sMode == 'BE') { imode <- as.integer(9) }

  return(imode)
}

mk_cmplxentries <- function(z)
{
  if (dim(z)[2] != 2)
    stop("wrong format from fortran return: dn_eigen_f")

  cmplx <- NULL
  cmplx <- sapply(1:length(z[ , 1]), function(x) { complex(real = z[x, 1], imaginary = z[x, 2]) })

  if (is.null(cmplx))
    stop("error while aggregating fortran return from two vectors to a complex")

  return(cmplx)
}


getEigenval <- function(values, mode, dim, nEig, symmetric){

  if (symmetric) {
    result <- rep(NA, dim)

    orderedInd <- order(values[1:nEig], decreasing = TRUE)
    values <- values[orderedInd]

    # sort works also for complex values
    if (identical(mode %% 2, 0)) { # even modes are SA, SM, SR, SI
      result[(dim - nEig + 1):dim] <- values
    } else {                       # odd modes are LA, LA, LR, LI
      result[1:nEig] <- values }

  } else {
    result <- matrix(NA, nrow = dim, ncol = 2)

    orderedInd <- order(values[1:nEig, 1], decreasing = TRUE)
    values <- values[orderedInd, ]

    if (identical(mode %% 2, 0)) {
      result[(dim - nEig + 1):dim, ] <- values
    } else {
      result[1:nEig, ] <- values }

    result <- mk_cmplxentries(result)
  }

  return(list("eigenvalues" = result , "order" = orderedInd))
}


getEigenvec <- function(v, sym, dimen, nEig, orderind, eigenvalues, cmplxeps){

  if (is.null(v))
    stop("fortran returned NULL eigenvectors")

  if (sym) {

    v <- matrix(v, nrow = dimen, ncol = nEig, byrow = FALSE)
    v <- v[ , orderind, drop = FALSE]

  } else {

    v <- matrix(v[1:(dimen*nEig*2)], nrow = dimen, ncol = nEig*2, byrow = FALSE)

    v_real <- v[ , seq(1, nEig*2, by = 2)]
    v_imag <- v[ , seq(2, nEig*2, by = 2)]

    v_real <- v_real[ , orderind, drop = FALSE]
    v_imag <- v_imag[ , orderind, drop = FALSE]

    rm(v)

    eigenvalues <- stats::na.omit(eigenvalues)
    v <- matrix(NA, nrow = dimen, ncol = nEig)
    v <- sapply(1:nEig, function(x) {
      if (abs(Im(eigenvalues[x])) > cmplxeps) {
          mk_cmplxentries(cbind(v_real[ , x], v_imag[ , x]))
      } else {
        mk_cmplxentries(cbind(v_real[ , x], base::rep.int(0, times = dimen)))
      }
    })

  }

  return(v)
}


eigen_approx <- function(x,
                         nev,
                         ncv,
                         nitr,
                         mode,
                         only.values  = FALSE,
                         verbose      = FALSE,
                         f_routine){

  nz <- x@rowpointers[x@dimension[1]+1]-1  
    if((nz==0) | (norm(x@entries) <= getOption("spam.eps"))){#trap zero matrix
        stop("'x' is an emtpy matrix")
    }
    
  # check & parse arguments
  if (x@dimension[1] <= nev)
    stop("nev: the number of eigenvalues to calculate must be smaller than the matrix dimensions", call. = TRUE)

  if (f_routine != "ds_eigen_f" && f_routine != "dn_eigen_f")
    stop("non valid fortran routine is specified", call. = TRUE)



  f_mode <- setMode(sMode = mode, symmetric = ifelse(identical(f_routine, "ds_eigen_f"), TRUE, FALSE))

  fortran_object <- result <- list(NULL)

  if(getOption("spam.force64") || .format.spam(x)$package != "spam") {
    SS <- .format64()
    f_routine <- paste0(f_routine, "64")
  } else {
    SS <- .format32 }

  # Fortran call: symmetric matrices
  if (identical(f_routine, "ds_eigen_f") || identical(f_routine, "ds_eigen_f64")) {

    # define upperbound for matrices in spam64, since ARPACK routines rely on BLAS/LAPACK functions
    # which have no integer/logical's of kind = 8
    if (max(x@dimension[1], ncv*(ncv + 8))*3 >= 2^31-1) {
      stop("the dimension and of the input matrix and/or the required number of eigenvalues
            is too large for the current eigenvalue decomposition implementation.\n
            Reducing the argument \"nev\" and\"ncv\" in the contorl options may helps.")
    }

    fortran_object <- .C64 (f_routine,
                            SIGNATURE = c("integer", "integer", "integer", "integer",
                                          "integer", "integer", "double",
                                          SS$signature, SS$signature, "double",
                                          "double", "integer"),
                            maxnev    = nev,
                            ncv       = ncv,
                            maxitr    = nitr,
                            n         = x@dimension[1],
                            iwhich    = f_mode,
                            na        = x@dimension[1],
                            a         = x@entries,
                            ja        = x@colindices,
                            ia        = x@rowpointers,
                            v         = vector_dc("double", x@dimension[1]*ncv),
                            d         = vector_dc("double", nev),
                            iparam    = integer_dc(8),
                            INTENT    = c("r", "r", "r", "r", "r", "r", "r", "r", "r",
                                           "rw", "rw", "rw"),
                            NAOK      = getOption("spam.NAOK"),
                            PACKAGE   = SS$package)

    if (is.null(fortran_object)) {
      stop("error while calling fortran routine, check (control) arguments", call. = TRUE) }

    result <-   list ("nEigenVal"     = nev,
                      "Mode"          = f_mode,
                      "Eigenvectors"  = if (!only.values) { fortran_object$v[1:(x@dimension[1]*nev)] } else { NULL },
                      "Eigenvalues"   = fortran_object$d,
                      "nconv"         = fortran_object$iparam[5])
  }

  # Fortran call: nonsymmetric matrices
  if (identical(f_routine, "dn_eigen_f") || identical(f_routine, "dn_eigen_f64")) {

    if (max(x@dimension[1], ncv^2+6*ncv)*3 >= 2^31-1) {
      stop("the dimension and of the input matrix and/or the required number of eigenvalues
            is too large for the current eigenvalue decomposition implementation.\n
            Reducing the argument \"nev\" and\"ncv\" in the contorl options may helps.")
    }


    fortran_object <- .C64 (f_routine,
                            SIGNATURE = c("integer", "integer", "integer", "integer",
                                          "integer", "integer", "double",
                                          SS$signature, SS$signature, "double",
                                          "double", "double", "integer"),
                            maxnev    = nev,
                            ncv       = ncv,
                            maxitr    = nitr,
                            n         = x@dimension[1],
                            iwhich    = f_mode,
                            na        = x@dimension[1],
                            a         = x@entries,
                            ja        = x@colindices,
                            ia        = x@rowpointers,
                            v         = vector_dc("double", x@dimension[1]*ncv),
                            dr        = vector_dc("double", nev+1),
                            di        = vector_dc("double", nev+1),
                            iparam    = integer_dc(8),
                            INTENT    = c("r", "r", "r", "r", "r", "r", "r", "r", "r",
                                          "rw", "rw", "rw", "rw"),
                            NAOK      = getOption("spam.NAOK"),
                            PACKAGE   = SS$package)

    if (is.null(fortran_object)) {
      stop("error while calling fortran routine, check (control) arguments", call. = TRUE) }

    result <-   list ("nEigenVal"     = nev,
                      "Mode"          = f_mode,
                      "Eigenvectors"  = if (!only.values) { fortran_object$v[1:(x@dimension[1]*nev*2)] } else { NULL },
                      "Eigenvalues"   = cbind(fortran_object$dr, fortran_object$di),
                      "nconv"         = fortran_object$iparam[5])
  }

  if (verbose) {
    cat("\n used options/arguments:")
    if (identical(f_routine, "dn_eigen_f") || identical(f_routine, "dn_eigen_f64")) {
      issym <- "non symmetric"
    } else {
      issym <- "symmetric" }
    cat(paste("\n      FORTRAN routine:", issym, "matrices"))
    cat(paste("\n      nitr:", nitr))
    cat(paste("\n      ncv:", ncv, "\n"))

    cat("\n summary FORTRAN return:")
    cat(paste("\n     ", nev, "eigenvalues requested and", result$nconv, "converged\n"))
  }

  if (result$nconv < nev)
    warning(paste("\n only", result$nconv, "instead of", nev ,"eigenvalues converged, try to increase 'control = list(nitr = .., ncv = ..)'"), call. = TRUE)

  if (is.null(result)) {
    stop("unpredicted error, check control options of the eigen.spam function.", call. = TRUE)
  }

  return(result)
}


eigen.spam <- function (x, nev = 10, symmetric, only.values = FALSE, control = list()){

  con <- list(mode     = 'LM',
              nitr     = NULL,
              ncv      = NULL,
              spamflag = FALSE,
              verbose  = FALSE,
              cmplxeps = .Machine$double.eps)

  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "), call. = TRUE)

  ifelse(!con$verbose, vFlag <- FALSE, vFlag <- TRUE)

  nz <- x@rowpointers[x@dimension[1]+1]-1  
  if((nz==0) | (norm(x@entries) <= getOption("spam.eps"))){#trap zero matrix
      stop("'x' is an emtpy matrix")
  }

  # arpack routines cant handle 'small' matrices
  minDimARPACK <- 100

  resContainer <- list(NULL)

  if (!is.spam(x)) {
    try(x <- as.spam(x)) }
  if (missing(symmetric)){
    symmetric <- isSymmetric.spam(x)
  } else {
    if (!identical(symmetric, isSymmetric.spam(x))) {
      warning(paste("the input matrix is", ifelse(isSymmetric.spam(x), "symmetric", "not symmetric")), call. = TRUE)
    }
  }

  ## dispatching between base::eigen and ARPACK
  if ((!con$spamflag && prod(dim(x)) <= getOption("spam.inefficiencywarning")) || (con$spamflag && prod(dim(x)) <= minDimARPACK)) {
    warning("\n The eigenvalues are calculated with the function 'eigen' from the base package.", call. = TRUE)

    if (con$spamflag) { warning(paste("\n Even 'spamflag = TRUE', since the matrix dimension is smaller than", minDimARPACK), call. = TRUE) }

    resContainer <- base::eigen(x           = x,
                                symmetric   = symmetric,
                                only.values = only.values)

    resEig       <- resContainer$values
    resVec       <- resContainer$vectors

  } else {
  # --- using ARPACK -------------------------------------------------------------- #

    if (!identical(x@dimension[1], x@dimension[2]))
      stop("non-square matrix in 'eigen'")

    if (nev >= x@dimension[1] || nev <= 0)
      stop ("the number of asked eigenvalues is higher or equal the dimension of the input matrix")

    if (symmetric) { ncvMaxMin <- 100 } else {  ncvMaxMin <- 100 }

    if (con$ncv > x@dimension[1] || con$ncv < nev || is.null(con$ncv)) {
      con$ncv <- min(x@dimension[1] + 1, max(2 * nev + 1, ncvMaxMin)) }

    if (is.null(con$nitr)) { con$nitr <- con$ncv + 1000 }

    # calculate eigenvalues and vectors
    resContainer <- eigen_approx(x            = x,
                                 nev          = nev,
                                 ncv          = con$ncv,
                                 nitr         = con$nitr,
                                 mode         = con$mode,
                                 only.values  = only.values,
                                 verbose      = vFlag,
                                 f_routine    = ifelse(symmetric, "ds_eigen_f", "dn_eigen_f"))

    if (is.null(resContainer)) {
      stop("\n FORTRAN return is empty") }

    tmpresEig    <-  getEigenval(values       = resContainer$Eigenvalues,
                                 mode         = setMode(con$mode, symmetric, silent = TRUE),
                                 dim          = x@dimension[1],
                                 nEig         = if (resContainer$nconv > nev) { nev } else { resContainer$nconv },
                                 symmetric    = symmetric)
    resEig <- tmpresEig$eigenvalues

    resVec <- NULL
    if(!only.values) {
      resVec     <-  getEigenvec(v            = resContainer$Eigenvectors,
                                 sym          = symmetric,
                                 dimen        = x@dimension[1],
                                 nEig         = if (resContainer$nconv > nev) { nev } else { resContainer$nconv },
                                 orderind     = tmpresEig$order,
                                 eigenvalues  = resEig,
                                 cmplxeps     = con$cmplxeps) }

  # ------------------------------------------------------------------------------- #
  }

  return (list("values" = resEig, "vectors" = resVec))
}


