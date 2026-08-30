deldir <- local({

digOutz <- function(z1,znm,x) {
    if(inherits(z1,"try-error")) { # z not found, so the call can't
                                   # have been z=NULL!
        if(inherits(x,"matrix")) {
            if(znm %in% colnames(x)) return(x[,znm])
        }
        if(znm %in% names(x)) return(x[[znm]])
        whinge <- paste0("Object z = ",znm," not found.\n")
        stop(whinge)
    }
    if(is.null(z1)) return(NULL)
# Here z1 is explicit vector, or a text string, equal to znm.
# If the former return it; if the latter look for the
# object named by that text string.
    if(!isTRUE(identical(z1,znm))) return(z1)
    z <- try(get(znm,pos=1),silent=TRUE)
    if(inherits(z,"try-error")) {
        if(inherits(x,"matrix")) {
            if(znm %in% colnames(x)) return(x[,znm])
        }
        if(znm %in% names(x)) return(x[[znm]])
        whinge <- paste0("Object z = ",znm," not found.\n")
        stop(whinge)
    }
    return(z)
}

digOutxy <- function(x,y,znm) {
if(inherits(x,c("data.frame","matrix"))) {
    whatsit <- if(inherits(x,"data.frame")) "data frame" else "matrix"
    jj   <- 1:ncol(x)
    vnms <- if(inherits(x,"data.frame")) names(x) else colnames(x)
    if(is.null(vnms)) vnms <- paste0("V",jj)
    jt <- jj[!(vnms %in% c("y",znm))]
    if(length(jt)==0) {
        whinge <- paste0("Argument \"x\" is a ",whatsit," but does not appear\n",
                         "  to contain the x-coordinates.\n")
        stop(whinge)
    }
    nomx <- min(jt)
    jx <- match("x",vnms,nomatch=nomx)
    vnms[jx] <- "x"
    y1 <- try(y,silent=TRUE)
    if(inherits(y1,"try-error") || is.null(y1)) {
        jy <- match("y",vnms)
        if(is.na(jy)) {
            jt <- jj[!(vnms %in% c("x",znm))]
            if(length(jt)==0) {
                whinge <- paste0("Argument \"x\" is a ",whatsit," but does not appear\n",
                                 "  to contain the y-coordinates, nor is \"y\" to be\n",
                                 "  found in the global environment.\n")
                stop(whinge)
            }
            jy <- min(jt)
        }
        return(list(x=x[,jx],y=x[,jy]))
    }
    return(list(x=x[,jx],y=y1))
}
if(inherits(x,"list")) {
    if(!("x" %in% names(x))) {
        whinge <- paste0("When argument \"x\" is a generic list, it must have\n",
                         "  a component named \"x\". \n")
        stop(whinge)
    }
    y1 <- try(y,silent=TRUE)
    if(inherits(y1,"try-error") || is.null(y)) {
        y1 <- x[["y"]]
    }
    if(!is.null(y1)) return(list(x=x[["x"]],y=y1))
    stop("Argument \"y\" not found.\n")
}
list(x=x,y=y)
}

function(x,y=NULL,z=NULL,rw=NULL,eps=1e-9,sort=TRUE,plot=FALSE,
                   round=TRUE,digits=6,id=NULL,...) {
# Function deldir to compute the Delaunay Triangulation (and hence
# the Dirichlet Tesselation) of a planar point set according to the
# second (iterative) algorithm of Lee and Schacter, International
# Journal of Computer and Information Sciences, Vol. 9, No. 3, 1980,
# pages 219 to 242.
#
# ORIGINALLY PROGRAMMED BY: Rolf Turner in 1987/88, while with the
# Division of Mathematics and Statistics, CSIRO, Sydney, Australia.
# Re-programmed by Rolf Turner to adapt the implementation from a
# stand-alone Fortran program to an S function, while visiting the
# University of Western Australia, May 1995.  Further revised
# December 1996.
# 

# The triangulation is made to be with respect to the whole plane by
# `suspending' it from `ideal' points
# (-R,-R), (R,-R) (R,R), and (-R,R), where R --> infinity.

# It is also enclosed in a finite rectangle (whose boundaries truncate any
# infinite Dirichlet tiles) with corners (xmin,ymin) etc.  This rectangle
# is referred to elsewhere as `the' rectangular window.

#cat("Fuck fuck fuck!!!\n")

# Organise the x, y and possibly z arguments.

if(inherits(x,"ppp")) {
# If the first argument is an object of class "ppp", extract the x
# and y coordinates from this object. If this object is "marked"
# and if the marks are atomic (a vector or a factor) and z is NULL,
# then set z equal to the marks.
    y1 <- try(y,silent=TRUE)
    if(!is.null(y1))
         warning("Since \"x\" is of class \"ppp\", argument \"y\" is ignored.\n")
    if(is.null(z)) {
        marx <- x$marks
        ok   <- !is.null(marx) & is.atomic(marx)
        if(ok) z <- marx
    }
    if(is.null(rw)) rw <- c(x$window$xrange, x$window$yrange)
    y <- x$y
    x <- x$x
} else {
    z1  <- try(z,silent=TRUE)
    if(inherits(z1,"character") & length(z1)==1) {
        znm <- z1
    } else {
        znm <- deparse(substitute(z))
    }
    xyTemp <- digOutxy(x,y,znm)
    z      <- digOutz(z1,znm,x)
    x      <- xyTemp$x
    y      <- xyTemp$y
}
haveZ <- !is.null(z)
# Check that x and y are numeric.
if(!is.numeric(x))
    stop("The x-coordinates must be numeric.\n")
if(!is.numeric(y))
    stop("The y-coordinates must be numeric.\n")

# Check that lengths match.
n <- length(x)
if(n!=length(y)) stop("Lengths of \"x\" and \"y\" do not match.\n")
if(haveZ) {
    if(n!=length(z))
        stop("Length of \"z\" does not match lengths of \"x\" and \"y\".\n")
}

# Check on the "id" argument.
haveId <- !is.null(id)
if(haveId) {
    if(any(duplicated(id)))
        stop("Argument \"id\", if supplied, must contain no duplicate values.\n")
    if(n!=length(id))
        stop("Length of \"id\" does not match lengths of \"x\" and \"y\".\n")
    id <- as.character(id)
}

# If a data window is specified, turn it into a 4-tuple (if necessary).
if(!is.null(rw)) {
   if(inherits(rw,"owin")) {
       xr <- rw$xrange
       yr <- rw$yrange
       rw <- c(xr,yr)

# Apparently --- according to Michael Chirico 24/09/2015 --- the
# following will accommodate the bounding box of a collection of
# polygons as structured in the "sp" package.
    } else if(is.matrix(rw)) {
        rw <- as.vector(t(rw))
    }
# Check that rw is an appropriate 4-tuple.
    ok <- length(rw)==4 && (rw[1] < rw[2] & rw[3] < rw[4])
    if(!ok)
        stop("The rectangula window rw is not of the appropriate form.\n")
}

# If a data window now exists, get its corner coordinates
# and truncate the data by this window.
if(!is.null(rw)) {
    xmin <- rw[1]
    xmax <- rw[2]
    ymin <- rw[3]
    ymax <- rw[4]
    ind.orig <- 1:n
    drop     <- ind.orig[x<xmin|x>xmax|y<ymin|y>ymax]
    if(length(drop)>0) {
        x  <- x[-drop]
        y  <- y[-drop]
        if(haveZ) z <- z[-drop]
        if(haveId) id <- id[-drop]
        ind.orig <- ind.orig[-drop]
    }
}
nn <- length(x) # Could be different from "n" if the data were
                # clipped to rw.

# If the rectangular window is (still) not specified, form its corners
# from the minimum and maximum of the data +/- 10%:
if(is.null(rw)) {
    if(length(x)==0)
        stop("No points nor any rectangular window specified.\n")
    xmin <- min(x)
    xmax <- max(x)
    ymin <- min(y)
    ymax <- max(y)
    xdff <- xmax-xmin
    ydff <- ymax-ymin
# Check for the setting in which all points lie on a vertical
# or horizontal straight line.
    if(isTRUE(all.equal(xdff,0))) {
        whinge <- paste0("The x-range of the points is zero, whence a",
                         " rectangular window\n  cannot be inferred from",
                         " the data.  You must specify the rectangular\n",
                         "  window explicitly.\n")
        stop(whinge)
    }
    if(isTRUE(all.equal(ydff,0))) {
        whinge <- paste0("The y-range of the points is zero, whence a",
                         " rectangular window\n  cannot be inferred from",
                         " the data.  You must specify the rectangular\n",
                         "  window explicitly.\n")
        stop(whinge)
    }
    xmin <- xmin-0.1*xdff
    xmax <- xmax+0.1*xdff
    ymin <- ymin-0.1*ydff
    ymax <- ymax+0.1*ydff
    rw   <- c(xmin,xmax,ymin,ymax)
    ind.orig <- 1:n
}

# Eliminate duplicate points:
iii <- duplicatedxy(x,y)
if(any(iii)) {
    kkk <- !iii
    nn  <- sum(kkk)
    if(haveZ) {
        jjj <- duplicated(data.frame(x=x,y=y,z=z))
        if(sum(jjj) < sum(iii)) {
            whinge <- paste("There were different tags corresponding to\n",
                            "duplicated points.\n",sep="")
            warning(whinge)
        }
        z   <- z[kkk]
    }
    x   <- x[kkk]
    y   <- y[kkk]
    ind.orig <- ind.orig[!iii]
}

# Check there are sufficiently many points to triangulate/tessellate.
if(nn <= 1) {
    whinge <- paste("There is at most one point inside the given\n",
                    " rectangular window. Thus there are insufficiently\n",
                    " many points to triangulate/tessellate.\n")
    stop(whinge)
}

# Sort the coordinates into "bins".  There are approximately
# sqrt(nn) such bins.  The vector "ind" (index) keeps track of the
# re-ordering; if ind[i] == j then i is the index of a point in
# the original sequence of points and j is the index of the same
# point in the bin sorted sequence.  The vector "rind" (reverse
# index) does the opposite; if rind[i] == j then i is the position
# of a point in the bin sorted sequence and j is its position in
# the original sequence.  Thus ind[rind[k]] = k and rind[ind[k]] = k
# for all k.  So xs[ind] (where xs is the bin sorted sequence of
# x's) is equal to x, the original sequence of x's.  Likewise ys[ind]
# (where ys is the bin sorted sequence of y's) is equal to y, the
# original sequence of y's. Conversely x[rind] = xs and y[rind] = ys.
#
# Added 20/03/2017:  I think I rather made a meal of this.  The
# vector of indices "rind" is just order(ind); ind is a permutation
# of 1, 2, ..., nn.  Thus rind[ind] = ind[rind] = 1, 2, ..., nn.
# However there's no real harm done, so I won't change the shaganappi
# code at this stage.
if(sort) {
    xy   <- binsrtR(x,y,rw)
    x    <- xy$x
    y    <- xy$y
    ind  <- xy$ind
    rind <- xy$rind
} else {
    ind  <- 1:nn
    rind <- 1:nn
}

# Make space for the total number of points as well as 4
# ideal points and 4 extra corner points which get used (only) by
# subroutines dirseg and dirout in the ``output'' process (returning
# a description of the triangulation after it has been calculated):
ntot <- nn + 4               # ntot includes the 4 ideal points but
                             # but NOT the 4 extra corners
x <- c(rep(0,4),x,rep(0,4))
y <- c(rep(0,4),y,rep(0,4))

# Set up fixed dimensioning constants:
ntdel  <- 4*nn
ntdir  <- 3*nn

# Set up dimensioning constants which might need to be increased:
madj <- max(20,ceiling(3*sqrt(ntot)))
tadj <- (madj+1)*(ntot+4)
ndel <- madj*(madj+1)/2
tdel <- 6*ndel
ndir <- ndel
tdir <- 10*ndir

# Call the master subroutine to do the work:
repeat {
    tmp <- .Fortran(
            'master',
            x=as.double(x),
            y=as.double(y),
            rw=as.double(rw),
            nn=as.integer(nn),
            ntot=as.integer(ntot),
            nadj=integer(tadj),
            madj=as.integer(madj),
            eps=as.double(eps),
            delsgs=double(tdel),
            ndel=as.integer(ndel),
            delsum=double(ntdel),
            dirsgs=double(tdir),
            ndir=as.integer(ndir),
            dirsum=double(ntdir),
            incAdj=integer(1),
            incSeg=integer(1),
            PACKAGE='deldir'
        )

# Check for problems with insufficient storage:
    incAdj <- tmp$incAdj
    incSeg <- tmp$incSeg
    if(incAdj==0 & incSeg==0) break
    if(incAdj==1) {
        nmadj <- ceiling(1.2*madj)
        wrds <- paste('Increasing madj from',madj,'to',nmadj,
                      'and trying again.')
        message(wrds)
        madj <- nmadj
        tadj <- (madj+1)*(ntot+4)
        ndel <- max(ndel,madj*(madj+1)/2)
        tdel <- 6*ndel
        ndir <- ndel
        tdir <- 10*ndir
    }
    if(incSeg==1) {
        nndel <- ceiling(1.2*ndel)
        wrds <-paste('Increasing ndel and ndir from',ndel,
                     'to',nndel,'and trying again.')
        message(wrds)
        ndel <- nndel
        tdel <- 6*ndel
        ndir <- ndel
        tdir <- 10*ndir
    }
}

# Collect up the results for return:
ndel             <- tmp$ndel
delsgs           <- if(round) {
                        round(t(as.matrix(matrix(tmp$delsgs,nrow=6)[,1:ndel])),digits)
                    } else {
                        t(as.matrix(matrix(tmp$delsgs,nrow=6)[,1:ndel]))
                    }
delsgs           <- as.data.frame(delsgs)
names(delsgs)    <- c('x1','y1','x2','y2','ind1','ind2')
delsum           <- matrix(tmp$delsum,ncol=4)
del.area         <- sum(delsum[,4])
delsum           <- if(round) {
                        round(cbind(delsum,delsum[,4]/del.area),digits)
                    } else {
                        cbind(delsum,delsum[,4]/del.area)
                    }
del.area         <- if(round) round(del.area,digits) else del.area
ndir             <- tmp$ndir
dirsgs           <- if(round) {
                        round(t(as.matrix(matrix(tmp$dirsgs,nrow=10)[,1:ndir])),digits)
                    } else {
                        t(as.matrix(matrix(tmp$dirsgs,nrow=10)[,1:ndir]))
                    }
dirsgs           <- as.data.frame(dirsgs)
dirsum           <- matrix(tmp$dirsum,ncol=3)
dir.area         <- sum(dirsum[,3])
dirsum           <- if(round) {
                        round(cbind(dirsum,dirsum[,3]/dir.area),digits)
                    } else {
                        cbind(dirsum,dirsum[,3]/dir.area)
                    }
dir.area         <- if(round) round(dir.area,digits) else dir.area
names(dirsgs)    <- c('x1','y1','x2','y2','ind1','ind2','bp1','bp2',
                      'thirdv1','thirdv2')
mode(dirsgs$bp1) <- 'logical'
mode(dirsgs$bp2) <- 'logical'
allsum           <- as.data.frame(cbind(delsum,dirsum))
names(allsum)    <- c('x','y','n.tri','del.area','del.wts',
                              'n.tside','nbpt','dir.area','dir.wts')

# The foregoing results are in terms of the indices of the bin sorted coordinates.
# Put things in terms of the indices of the original coordinates.
delsgs$ind1 <- rind[delsgs$ind1]
delsgs$ind2 <- rind[delsgs$ind2]
dirsgs$ind1 <- rind[dirsgs$ind1]
dirsgs$ind2 <- rind[dirsgs$ind2]
dirsgs$thirdv1  <- with(dirsgs,ifelse(thirdv1<0,thirdv1,rind[abs(thirdv1)]))
dirsgs$thirdv2  <- with(dirsgs,ifelse(thirdv2<0,thirdv2,rind[abs(thirdv2)]))

# If "id" was supplied, change the entries of "ind1", "ind2", "thirdv1",
# and "thirdv2" to be the corresponding entries of "id".
if(haveId) {
    delsgs$ind1 <- id[delsgs$ind1]
    delsgs$ind2 <- id[delsgs$ind2]
    dirsgs$ind1 <- id[dirsgs$ind1]
    dirsgs$ind2 <- id[dirsgs$ind2]
    dirsgs$thirdv1  <- with(dirsgs,ifelse(thirdv1<0,as.character(thirdv1),
                                          id[abs(thirdv1)]))
    dirsgs$thirdv2  <- with(dirsgs,ifelse(thirdv2<0,as.character(thirdv2),
                                          id[abs(thirdv2)]))
}

# The points in "allsum" appear in bin-sorted order; rearrange
# the rows of "allsum" so that the points appear in the original order.
allsum          <- allsum[ind,]

# The following is a furphy --- it just makes the rownames into
# 1, 2, ..., n.  At this point the rownames of "allsum" were
# (1:n)[ind].  So we're getting (1:n)[ind][rind] = ind[rind]
# = 1:n !!!
# rownames(allsum) <- rownames(allsum)[rind]
# So we could just set rownames(allsum) <- 1:nrow(allsum) and
# get the same effect.  However that does not take account of
# *duplicated* points.  So it is better to use ind.orig.  Note that
# the resulting rowname corresponding to a point is the index (in
# the original sequence of points) of the *first* in its sequence
# of duplicated points.

rownames(allsum) <- ind.orig

# Arrange for the tags and id to be in the summary, given
# that they were supplied.
dfz <- if(haveZ) {
    data.frame(z=z)
} else {
    as.data.frame(matrix(nrow=nn,ncol=0))
}
dfid <- if(haveId) {
    data.frame(id=id)
} else {
    as.data.frame(matrix(nrow=nn,ncol=0))
}

allsum <- cbind(allsum[,1:2],dfid,dfz,allsum[,3:9])
rw     <- if(round) round(rw,digits) else rw

# Aw' done!!!
rslt <- list(delsgs=delsgs,dirsgs=dirsgs,summary=allsum,n.data=nn,
             del.area=del.area,dir.area=dir.area,rw=rw,
             ind.orig=ind.orig)
attr(rslt,"round") <- round
attr(rslt,"digits") <- if(round) digits else NA
class(rslt) <- 'deldir'
if(plot) {
    plot(rslt,...)
    return(invisible(rslt))
} else return(rslt)
}
})
