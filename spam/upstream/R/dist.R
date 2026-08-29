# HEADER ####################################################
# This is file spam/R/dist.R.                               #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################


### in base:
# dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p=2)
#
### in fields
# rdist( x1, x2)
# rdist.earth(x1, x2, miles = TRUE, R = NULL)
# fields.rdist.near( x1, x2, delta, max.points= NULL)
#
### in sp
# spDistsN1(pts, pt, longlat=FALSE)
#
### in amap     ### nbproc  integer, Number of subprocess for parallelization
# Dist(x, method = "euclidean", nbproc = 1, diag = FALSE, upper = FALSE)
#
### in argosfilter
# distance(lat1, lat2, lon1, lon2)  # gc between two pts in km
# distanceTrack(lat,lon)            # gc between pts  in km
#
### in proxy
# dist(x, y = NULL, method = NULL, ..., diag = FALSE, upper = FALSE,
#     pairwise = FALSE, by_rows = TRUE, convert_similarities = TRUE,
#     auto_convert_data_frames = TRUE)
#
### in RFOC
#  GreatDist(LON1, LAT1, LON2, LAT2, EARTHRAD= 6371)
#

nearest.dist <- function( x, y=NULL, method = "euclidean",
                         delta = 1,
                         upper = if(is.null(y)) FALSE else NULL,
                         p = 2, miles=TRUE, R=NULL, fortran = FALSE
#                         eps =  NULL, diag = NULL
                         )
{
  # see help for exact parameter meaning

  # We always include all small distances. Hence, this function
  #   works differently than any other spam functions. An additional
  #   call to an as.spam would eliminate the small values.
#  if (!is.null(diag)) warning("Argument "diag" is deprecated")
#  if (!is.null(eps))  warning("Argument "eps" is deprecated")

  if (!is.na(pmatch(method, "euclidian")))     method <- "euclidean"
  METHODS <- c("euclidean", "maximum", "minkowski", "greatcircle")
  method <- pmatch(method, METHODS)  # result is integer

  if (is.na(method))     stop("invalid distance method")

  if (method == 4) {
    if (is.null(R))
      p <- ifelse( miles,3963.34,6378.388)
    else {
      if (R <= 0)           stop("'R' should be postiive")
      p <- R
    }
    if (abs(delta)>180.1)  stop("'delta' should be smaller than 180 degrees.")
  }

  if (is.null(upper))
    part <- 0L
  else
    part <- ifelse(upper, 1L ,-1L)
  if (is.data.frame(x))  x <- as.matrix(x)
  if (is.list(x)) stop("'x' should be an array or matrix")
  # as.matrix( list() ) does not work
  if (!is.matrix(x)) x <- as.matrix(x)
  nd <- dim(x)[2]
  n1 <- dim(x)[1]

  if (!is.null(y)) {
    # we specify x and y:
    if (is.data.frame(y))  y <- as.matrix(y)
    if (is.list(x)) stop("'x' should be an array or matrix")
    if (!is.matrix(y)) y <- as.matrix(y)
    if (nd!=dim(y)[2]) stop("'x' and 'y' should have the same number of columns.")
    n2 <- dim(y)[1]
    sizemax <- n2*n1
    # there is an as.double just in case that mi (and n1 below) is > 2^16
  } else {
    # x = y, i.e. proper distance matrix
    if (n1==1)         stop("More than a single point in 'x' is required.")
    if (method == 4) {
      p <- -p  # we save one variable...
    }
    y <- x
    n2 <- n1
    sizemax <- n1*(n1+1)/2
  }

  # which code should be used
  force64 <- getOption("spam.force64")
  if (force64 || fortran || sizemax > 2147483647 - 1 ){
    if (!is.null(y)) {
      mi <- min(n1,n2)
      ma <- max(n1,n2)
      nnz <- min(max(getOption("spam.nearestdistnnz")[1],
                     ma*getOption("spam.nearestdistnnz")[2]),
             (as.double(mi)*(mi+1)+(ma-mi)^2)/ ifelse( is.null(upper), 1, 2))
      } else {
        nnz  <- min(max(getOption("spam.nearestdistnnz")[1],
                        ma*getOption("spam.nearestdistnnz")[2]),
                  (as.double(n1)*(n1+1))/ ifelse( is.null(upper), 1, 2))
      }
   # EXPLICIT-STORAGE-FORMAT
    if(max(n1,n2,nnz) > 2147483647 - 1 || force64)
      SS <- .format64()
    else
      SS <- .format32

    repeat {
      d <- NULL # Free the memory allocated by a previous attempt
      d <- .C64("closestdist",
                ##              subroutine closestdist( ncol, x,nrowx, y, nrowy,
                ##  &    part, p, method,
                ##  &    eta, colindices, rowpointers, entries, nnz, iflag)
                SIGNATURE=c(SS$signature, "double", SS$signature, "double", SS$signature,
                            SS$signature, "double", SS$signature, "double", SS$signature,
                            SS$signature, "double", SS$signature, SS$signature),

                nd,
                x,
                n1, #w
                y,

                n2, #w
                part, #arg 6
                p[1],
                method,

                abs( delta[1]),
                colindices=vector(SS$type, nnz),
                rowpointers=vector(SS$type, n1+1), #arg 11
                entries=vector("double",nnz),

                nnz=nnz,
                iflag=0,

                INTENT = c("r", "r", "r", "r",
                           "r", "r", "rw", "r",
                           "r", "w", "w", "w",
                           "rw", "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE=SS$package)

     if (d$iflag==0) break else {
       nnz <-  nnz*getOption("spam.nearestdistincreasefactor")*n1/(d$iflag-1)

        # EXPLICIT-STORAGE-FORMAT
        if(max(n1,n2,nnz) > 2147483647 - 1 || force64)
         SS <- .format64()
        else
          SS <- .format32

      # if(nnz > 2147483647) stop("Distance matrix is too dense (2)")

        madens <- d$iflag
        warning(paste("You ask for a 'dense' sparse distance matrix, I require one more iteration.",
                            "\nTo avoid the iteration, increase 'nearestdistnnz' option to something like\n",
                            "'options(spam.nearestdistnnz=c(",d$nnz,",500))'\n(constructed ",madens,
                            " lines out of ",n1,").\n",sep=""), call. = TRUE)

      }
    }

   length(d$entries) <- d$nnz
   length(d$colindices) <- d$nnz
  } else {
    d <- closestdistCpp(x,y,part,p,method,eta=delta)
  }

  if(length(d$entries) == 0) {
    return(.newSpam(
      dimension=c(n1,n2),
      force64 = force64
    ))
  }

  return(.newSpam(
    entries=d$entries,
    colindices=d$colindices,
    rowpointers=d$rowpointers,
    dimension=c(n1,n2),
    force64 = force64
  ))
}

# in fields:
# rdist <- function (x1, x2)

spam_rdist <- function(x1, x2, delta = 1)
     nearest.dist(x1, y=x2,   delta = delta,  upper = NULL)

# in fields:
# rdist.earth <- function (x1, x2, miles = TRUE, R = NULL)
spam_rdist.earth <- function(x1, x2, delta=1, miles = TRUE, R = NULL)
    nearest.dist( x1, y=x2, method = "greatcircle",
                         delta = delta, miles=miles, R=R,  upper = NULL)

