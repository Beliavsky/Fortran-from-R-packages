insidePoly <- function(x,y,pgon,sanityCheck=FALSE,
                       tolerance=sqrt(.Machine$double.eps)) {
# This code is taken from the spatstat.utils function
# inside.xypolygon() (the "interpreted" method), with minor
# modifications, particularly in respect of the "sanity check".
#
# Also, in respect of determing boundary points, arrangements have
# been made to allow conversion of` "==" to almost.equal() ("%~%";
# see below).  This permits better determination of boundary points.

# Check that pgon is of the appropriate form.
if(inherits(pgon,"list")) {
    ok <- all(c("x","y") %in% names(pgon))
    if(!ok) {
        stop("Argument \"pgon\" must be a list with components \"x\" and \"y\".\n")
    }
    xp <- pgon$x
    yp <- pgon$y
    ok1 <- is.numeric(xp)
    ok2 <- is.numeric(yp)
    ok3 <- length(xp) == length(yp)
    if(!(ok1 & ok2 & ok3)) {
        stop("The components of \"pgon\" are not of the correct form.\n")
    }
} else {
    stop("Argument \"pgon\" must be a list determining a polygon.\n")
}

if(tolerance > 0) {
    almost.equal <- function(x, y)  {
        abs(x - y) < tolerance
    }
    assign("tolerance",tolerance,envir=environment(almost.equal))
    `%~%`  <- function(x,y) almost.equal(x,y)
} else {
    `%~%` <- `==`
}

full.npts <- npts <- length(x)
nedges <- length(xp)   # sic
vx <- x%in%xp
vy <- y%in%yp
vv <- vx & vy
retain <- !vv

# Remove vertices from subsequent consideration; replace them later
if(vertices.present <- !all(retain)) {
  x <- x[retain]
  y <- y[retain]
  npts <- sum(retain)
}

score <- numeric(npts)
on.boundary <- rep.int(FALSE, npts)

if(anyretain <- any(retain)) {
    for(i in 1:nedges) {
        x0 <- xp[i]
        y0 <- yp[i]
        x1 <- if(i == nedges) xp[1] else xp[i+1]
        y1 <- if(i == nedges) yp[1] else yp[i+1]
        dx <- x1 - x0
        dy <- y1 - y0
        if(dx < 0) {
# upper edge
            xcriterion <- (x - x0) * (x - x1)
            consider <- (xcriterion <= 0)
            if(any(consider)) {
                ycriterion <- y[consider] * dx - x[consider] * dy +  x0 * dy - y0 * dx
# closed inequality
                contrib <- (ycriterion >= 0) *
                   ifelse(xcriterion[consider] %~% 0, 1/2, 1)
# positive edge sign
                score[consider] <- score[consider] + contrib
# detect whether any point lies on this segment
                on.boundary[consider] <- on.boundary[consider] | (ycriterion %~% 0)
            }
        } else if(dx > 0) {
# lower edge
            xcriterion <- (x - x0) * (x - x1)
            consider <- (xcriterion <= 0)
            if(any(consider)) {
                ycriterion <- y[consider] * dx - x[consider] * dy + x0 * dy - y0 * dx
# open inequality
                contrib <- (ycriterion < 0) * ifelse(xcriterion[consider] %~% 0, 1/2, 1)
# negative edge sign
                score[consider] <- score[consider] - contrib
# detect whether any point lies on this segment
                on.boundary[consider] <- on.boundary[consider] | (ycriterion %~% 0)
            }
        } else {
# vertical edge
            consider <- (x %~% x0)
            if(any(consider)) {
# zero score
# detect whether any point lies on this segment
                yconsider <- y[consider]
                ycriterion <- (yconsider - y0) * (yconsider - y1)
                on.boundary[consider] <- on.boundary[consider] | (ycriterion <= 0)
            }
       }
    }
}

# replace any polygon vertices that were temporarily removed
if(vertices.present) {
  full.score <- numeric(full.npts)
  full.on.boundary <- rep.int(FALSE, full.npts)
  if(anyretain) {
    full.score[retain] <- score
    full.on.boundary[retain] <- on.boundary
  }
  full.score[vv] <- 1
  full.on.boundary[vv] <- TRUE
  score       <- full.score
  on.boundary <- full.on.boundary
  npts        <- full.npts
  
}
  
# any point recognised as lying on the boundary gets score 1.
score[on.boundary] <- 1

# Sanity Clause.
if(sanityCheck) {
    if(!all((score == 0) | (score == 1)))
        warning("Some \"scores\" are neither equal to 0 nor to 1.\n")
}

# Aw' done!
score <- as.logical(score)
attr(score,"on.boundary") <- on.boundary
score
}
