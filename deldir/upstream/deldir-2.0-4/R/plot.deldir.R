plot.deldir <- local({

fixColours <- function(cmpnt_col) {
    cmpnt_nms <- c("tri","tess","points","labels","rect")
    if(is.null(cmpnt_col)) {
        cmpnt_col <- rep(1,5)
        names(cmpnt_col) <- cmpnt_nms
    } else {
        cmpnt_col <- unlist(cmpnt_col)
        if(length(cmpnt_col) > 5) cmpnt_col <- cmpnt_col[1:5]
        if(!is.null(names(cmpnt_col))) {
            if(!all(names(cmpnt_col) %in% cmpnt_nms)) {
                stop("Argument \"cmpnt_col\" has incorrect names.\n")
            }
            ctmp <- rep(NA,5)
            names(ctmp) <- cmpnt_nms
            ctmp[names(cmpnt_col)] <- cmpnt_col
            cmpnt_col <- ctmp
        } else {
            cmpnt_col <- rep(cmpnt_col,length.out=5)
            names(cmpnt_col) <- cmpnt_nms
        }
        if(any(is.na(cmpnt_col))) {
            mde <- mode(cmpnt_col)
            switch(EXPR=mde,
               character={cmpnt_col[is.na(cmpnt_col)] <- palette()[1]},
               numeric={cmpnt_col[is.na(cmpnt_col)] <- 1},
               stop("Argument \"cmpnt_col\" is of an inappropriate mode.\n")
            )
        }
    }
    cmpnt_col
}

fixLines <- function(cmpnt_lty) {
    lty_nms <- c("tri","tess")
    if(is.null(cmpnt_lty)) {
        cmpnt_lty <- 1:2
        names(cmpnt_lty) <- lty_nms
    } else {
        cmpnt_lty <- unlist(cmpnt_lty)
        if(length(cmpnt_lty) > 2) cmpnt_lty <- cmpnt_lty[1:2]
        if(mode(cmpnt_lty) == "numeric") {
            if(!all(cmpnt_lty %in% 1:6)) {
                whinge <- paste("Numeric values of \"cmpnt_lty\" must",
                                "be integers between 1 and 6.\n")
                stop(whinge)
            }
        } else if(mode(cmpnt_lty) == "character") {
            linetypes <- c("solid","dashed","dotted","dotdash",
                           "longdash","twodash") 
            if(!all(cmpnt_lty %in% linetypes)) {
                whinge <- paste0("Text string values of \"cmpnt_lty\" must ",
                                "be one of the strings\n",paste(linetypes,collapse=", "),
                                ".\n")
                stop(whinge)
            }
        } else {
            whinge <- paste0("Argument \"cmpnt_lty\" must be of mode either",
                            " \"numeric\" or \"character\".\n")
            stop(whinge)
        }
        if(!is.null(names(cmpnt_lty))) {
            if(!all(names(cmpnt_lty) %in% lty_nms)) {
                stop("Argument \"cmpnt_lty\" has incorrect names.\n")
            }
            ltmp <- rep(NA,2)
            names(ltmp) <- lty_nms
            ltmp[names(cmpnt_lty)] <- cmpnt_lty
            cmpnt_lty <- ltmp
            dflt <- if(mode(cmpnt_lty) == "character") "solid" else 1
            if(any(is.na(cmpnt_lty))) {
                cmpnt_lty[is.na(cmpnt_lty)] <- dflt
            }
        } else {
            cmpnt_lty <- rep(cmpnt_lty,length.out=2)
            names(cmpnt_lty) <- lty_nms
        }
    }
    cmpnt_lty
}

function(x,add=FALSE,wlines=c('both','triang','tess'),
                        showpoints=TRUE,labelPts=FALSE,cex=1,lex=1,
                        cmpnt_col=NULL,cmpnt_lty=NULL,pch=1,
                        xlim=NULL,ylim=NULL,axes=FALSE,
                        xlab=if(axes) 'x' else '',
                        ylab=if(axes) 'y' else'',
                        showrect=FALSE,asp=1,...)
{
#
# Function plot.deldir to produce a plot of the Delaunay triangulation
# and Dirichlet tesselation of a point set, as produced by the
# function deldir().
#

# Check for use of the defunct argument name "number".
ccc <- match.call()
i   <- match("number",names(ccc))
if(!is.na(i)) {
    if("labelPts" %in% names(ccc)) {
        whinge <- paste0("Both \"labelPts\" and the defunct argument",
                         " \"number\" have been\n  specified.  Do not use",
                         " the defunct argument \"number\".  Use\n",
                         "  \"labelPts\" only.\n")
        stop(whinge)
    }
    whinge <- paste0("The argument name \"number\" is defunct. Please",
                     " use \"labelPts\"\n  instead.\n")
    warning(whinge)
    names(ccc)[i] <- "labelPts"
    return(eval(ccc))
}

# Carry on.

# Check that x is of the appropriate class.
if(!inherits(x, "deldir")) 
        stop("Argument \"x\" is not of class deldir.\n")

wlines    <- match.arg(wlines)
cmpnt_col <- fixColours(cmpnt_col)
cmpnt_lty <- fixLines(cmpnt_lty)

plot.del <- switch(wlines,both=TRUE,triang=TRUE,tess=FALSE)
plot.dir <- switch(wlines,both=TRUE,triang=FALSE,tess=TRUE)

delsgs <- x$delsgs
dirsgs <- x$dirsgs
n      <- x$n.data
rw     <- x$rw

if(plot.del) {
    x1<-delsgs[,1]
    y1<-delsgs[,2]
    x2<-delsgs[,3]
    y2<-delsgs[,4]
} else {
    x1 <- y1 <- x2 <- y2 <- numeric(0)
}

if(plot.dir) {
    u1<-dirsgs[,1]
    v1<-dirsgs[,2]
    u2<-dirsgs[,3]
    v2<-dirsgs[,4]
} else {
    u1 <- v1 <- u2 <- v2 <- numeric(0)
}

X<-x$summary[,"x"]
Y<-x$summary[,"y"]

if(!add) {
    if(is.null(xlim)) xlim <- rw[1:2]
    if(is.null(ylim)) ylim <- rw[3:4]
        x.all <- c(x1,x2,u1,u2,X)
        y.all <- c(y1,y2,v1,v2,Y)
        pty <- list(...)$pty
        if(!is.null(pty)) {
            OP <- par(pty=pty)
            on.exit(par(OP))
        }
    plot(x.all,y.all,type='n',xlim=xlim,ylim=ylim,
         xlab=xlab,ylab=ylab,axes=axes,asp=asp)
}

if(plot.del) {
    dotargs <- list(...)
    if(is.null(dotargs$col)) dotargs$col <- cmpnt_col[1]
    if(is.null(dotargs$lty)) dotargs$lty <- cmpnt_lty[1]
    arhgs <- c(list(x1,y1,x2,y2),dotargs)
    do.call(segments,arhgs)
}
if(plot.dir) {
    dotargs <- list(...)
    if(is.null(dotargs$col)) dotargs$col <- cmpnt_col[2]
    if(is.null(dotargs$lty)) dotargs$lty <- cmpnt_lty[2]
    arhgs <- c(list(u1,v1,u2,v2),dotargs)
    do.call(segments,arhgs)
}
if(showpoints) {
    dotargs <- list(...)
    dotargs$pch <- NULL
    dotargs$cex <- NULL
    do.call(points,c(list(x=X,y=Y,pch=pch,col=cmpnt_col[3],cex=cex),dotargs))
}
if(labelPts) {
    xoff <- 0.02*diff(range(X))
    yoff <- 0.02*diff(range(Y))
    dotargs <- list(...)
    dotargs$ces <- NULL
    dotargs$col <- NULL
    lbls <- x$summary[["id"]]
    if(is.null(lbls)) lbls <- 1:nrow(x$summary)
    do.call(text,c(list(x=X+xoff,y=Y+yoff,labels=lbls,cex=lex,
            col=cmpnt_col[4]),dotargs))
}
if(showrect) do.call(rect,c(as.list(x$rw)[c(1,3,2,4)],list(border=cmpnt_col[5])))
invisible()
}
})
