# updated in february 2013 by Raimon
pairwisePlot <- function(X,Y,...) UseMethod("pairwisePlot",X)

# updated in february 2013 by Raimon                                 
pairwisePlot.default <- function(X,Y=X,...,xlab=deparse(substitute(X)),ylab=deparse(substitute(Y)),nm=c(length(Y),length(X)),panel=plot,
                                 add.line=FALSE, line.col=2,add.robust=FALSE,rob.col=4) {
  ylab
  xlab
  if( is.data.frame(X) ) {
    X<-X # data.frames are lists already
  } else if( is.list(X) ) {
    if( !missing(xlab) || is.null(names(X)) )
      names(X) <- if( length(xlab) < length(X) ) paste(xlab,1:length(X),sep="") else xlab
  } else {
    if( length(dim(X))== 0 ) {
      X <- t(t(X))
      colnames(X) <- xlab
    }
    if( is.null(colnames(X)) )
      colnames(X) <- if( length(xlab) < ncol(X) ) paste(xlab,1:ncol(X),sep="") else xlab
    cn <- 1:ncol(X)
    names(cn)<-colnames(X)
    X <- lapply(cn,function(i) X[,i])
  }
                                        #else {
                                        #stop("Unhandled type of Y in pairwisePlot:",class(Y))
                                        #}
  if( is.data.frame(Y) ){
    Y<-Y # data.frames are lists already
  } else if( is.list(Y) ) {
    if( !missing(ylab) || is.null(names(Y)) )
      names(Y) <- if( length(ylab) < length(Y) ) paste(ylab,1:length(Y),sep="") else ylab
  } else {
    if( length(dim(Y))== 0 ) {
      Y <- t(t(Y))
      colnames(Y) <- ylab
    }
    if( is.null(colnames(Y)) )
      colnames(Y) <- if( length(ylab) < ncol(Y) ) paste(ylab,1:ncol(Y),sep="") else ylab
    cn <- 1:ncol(Y)
    names(cn)<-colnames(Y)
    Y <- lapply(cn,function(i) Y[,i])
  }
                                        #else {
                                        # stop("Unhandled type of Y in pairwisePlot:",class(Y))
                                        #}
  if( !is.null(nm) ) {
    opar <- par(mfrow=nm)
    on.exit(par(opar))
  }
  withformula = length( grep("formula", methods(panel) ) )>0
  for(j in 1:length(Y) )
    for(i in 1:length(X)) {
      # if there is a formula interface for the desired panel, use it
      if(withformula){
        panel(Y[[j]]~X[[i]],
            xlab=names(X)[i],
            ylab=names(Y)[j],...)
      }else{
        panel(X[[i]], Y[[j]],
              xlab=names(X)[i],
              ylab=names(Y)[j],...)
      }  
      if( !any(is.factor(Y[[j]]), is.factor(X[[i]]) ) ){
       if(add.line) abline( lm(Y[[j]]~X[[i]]), col=line.col )
       if(add.robust) abline(lmrob(Y[[j]]~X[[i]])$coefficients, col=rob.col, lwd=2)
      } 
    }
}

# updated in february 2013 by Raimon
pwlrPlot = function(x, y, ..., add.line=FALSE, line.col=2, add.robust=FALSE, rob.col=4){
  if(is.null(y)){
    stop("argument y cannot be empty")
  }
  if(is.null(x)){
    stop("argument x cannot be empty")
  }
  tk = "acomp" %in% c(class(x), class(y))
  if(!tk){
    stop("one of x or y must be an acomp composition!")
  }
  X = list(x,y)
  # check which of the two is not compositional, and it is a vector
  funcond = function(z) any(is.data.frame(z) & ncol(z)==1, 
                            is.matrix(z) & ncol(z)==1, 
                            (is.factor(z) || is.numeric(z)) && class(z)!="acomp" )
  tk = sapply(X, funcond)
  covarisfactor = is.factor(X[tk][[1]])
  if(!any(tk)){
    stop("one of x or y must be a unique covariable, either a one-column data.frame or matrix, a factor or a vector")
  }
  if( any(add.line, add.robust) & covarisfactor){
    add.line = FALSE
    add.robust = FALSE
    warning("regression line meaningless with factor covariable; add.line and add.robust set to FALSE")
  }
  # apply a pwlr to the compositional part
  D = ncol(X[!tk][[1]])
  partnames = colnames(X[!tk][[1]])
  covariablename = c(deparse(substitute(x)), deparse(substitute(y)))[tk]
  pwlr = function(i,j) log( X[!tk][[1]][,i]/X[!tk][[1]][,j] )
  X[!tk][[1]] <- mapply(i=rep(1:D,times=D),j=rep(1:D, each=D), pwlr)
  # ensure that the covariable is represented as a vector, then repeat as many times as D*D
  aux =  unlist(X[tk][[1]], recursive=FALSE)
  dim(aux) = NULL
  X[tk][[1]] = rep(aux, times=D*D)
  dim(X[tk][[1]]) = dim(X[!tk][[1]])
  # set graphical parameters now, store prior back as default
  opar <- par()
  on.exit(par(opar))
  par(mfrow=c(D,D), mar=c(1,1,0,0),oma=c(3,3,2,2))
  # plot!
  for(i in 1:D){
    for(j in 1:D){
      k = (i-1)*D+j
      if(i!=j){
        plot( X[[2]][,k]~X[[1]][,k], xaxt=ifelse(i==D,"s","n"), yaxt=ifelse(j==1,"s","n"),... )
        if(add.line) abline(lm( X[[2]][,k]~X[[1]][,k] ), col=line.col, lwd=2)
        if(add.robust) abline(lmrob(X[[2]][,k]~X[[1]][,k] )$coefficients, col=rob.col, lwd=2)
        if(!is.factor(X[[1]][,k]) ) axis(side=3, labels=(i==1))
        if(i==1 & is.factor(X[[1]][,k]) ) axis(side=3, at=1:length(levels( X[[1]][,k])), labels=abbreviate(levels( X[[1]][,k] )))
        axis(side=4, labels=(j==D))
      }else{
        plot( 0, 0, xaxt="n", yaxt="n", bty="n", type="n", xlab="", ylab=""  )
        text(0,0, labels=partnames[i],cex=1.5)
      }
    }
  }
  titols = rep("pairwise logratio", 2)
  titols[tk] = covariablename
  mtext(side=1:2, titols, cex=1.25, outer=TRUE, line=1.5)
}



#### pairs plots -----
pairs.acomp <-  function(x, labels, 
      panel = vp.lrdensityplot, ...,
      horInd = 1:ncol(x), verInd = 1:ncol(x),
      lower.panel = panel, upper.panel = panel,
      diag.panel = NULL, text.panel = textPanel,
      label.pos = 0.5 + has.diag/3, line.main = 3,
      cex.labels = NULL, font.labels = 1,
      row1attop = TRUE, gap = 1, log = ""){
  if (doText <- missing(text.panel) || is.function(text.panel)) 
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    xpd <- NA
    if (side%%2L == 1L && xl[j]) 
      xpd <- FALSE
    if (side%%2L == 0L && yl[i]) 
      xpd <- FALSE
    if (side%%2L == 1L) 
      Axis(x, side = side, xpd = xpd, ...)
    else Axis(y, side = side, xpd = xpd, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 3L) 
    stop("only 1 or 2 columns in the argument to 'pairs'")
  if (!all(horInd >= 1L & horInd <= nc)) 
    stop("invalid argument 'horInd'")
  if (!all(verInd >= 1L & verInd <= nc)) 
    stop("invalid argument 'verInd'")
  if (doText) {
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      doText <- FALSE
  }
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  main <- if ("main" %in% nmdots) 
    dots$main
  if (is.null(oma)) 
    oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
  opar <- par(mfrow = c(length(horInd), length(verInd)), mar = rep.int(gap/2, 
                                                                       4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  xl <- yl <- logical(nc)
  if (is.numeric(log)) 
    xl[log] <- yl[log] <- TRUE
  else {
    xl[] <- grepl("x", log)
    yl[] <- grepl("y", log)
  }
  for (i in if (row1attop) 
    verInd
    else rev(verInd)) for (j in horInd) {
      l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", 
                                                 ""))
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ..., log = l)
      xyl = sapply(1:ncol(x), function(i) -c(alr(x, ivar = i)))
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        box()
        if (i == 1 && (!(j%%2L) || !has.upper || !has.lower)) 
          localAxis(1L + 2L * row1attop, xyl[, j], xyl[, i], 
                    ...)
        if (i == nc && (j%%2L || !has.upper || !has.lower)) 
          localAxis(3L - 2L * row1attop, xyl[, j], xyl[, i], 
                    ...)
        if (j == 1 && (!(i%%2L) || !has.upper || !has.lower)) 
          localAxis(2L, xyl[, j], xyl[, i], ...)
        if (j == nc && (i%%2L || !has.upper || !has.lower)) 
          localAxis(4L, xyl[, j], xyl[, i], ...)
        mfg <- par("mfg")
        if (i == j) {
          if (has.diag) 
            localDiagPanel(as.vector(x[, i]), ...)
          if (doText) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            xlp <- if (xl[i]) 
              10^0.5
            else 0.5
            ylp <- if (yl[j]) 
              10^label.pos
            else label.pos
            text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                       font = font.labels)
          }
        }
        else if (i < j) 
          localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                         i]), ...)
        else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                            i]), ...)
        if (any(par("mfg") != mfg)) 
          stop("the 'panel' function made a new plot")
      }
      else par(new = FALSE)
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
          font = font.main)
  }
  invisible(NULL)
}



pairs.rcomp <- function(x, labels, 
               panel = vp.diffdensityplot,
               ...,
               horInd = 1:ncol(x), verInd = 1:ncol(x),
               lower.panel = panel, upper.panel = panel,
               diag.panel = NULL, text.panel = textPanel,
               label.pos = 0.5 + has.diag/3, line.main = 3,
               cex.labels = NULL, font.labels = 1,
               row1attop = TRUE, gap = 1, log = ""){
  if (doText <- missing(text.panel) || is.function(text.panel)) 
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    xpd <- NA
    if (side%%2L == 1L && xl[j]) 
      xpd <- FALSE
    if (side%%2L == 0L && yl[i]) 
      xpd <- FALSE
    if (side%%2L == 1L) 
      Axis(x, side = side, xpd = xpd, ...)
    else Axis(y, side = side, xpd = xpd, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 3L) 
    stop("only 1 or 2 columns in the argument to 'pairs'")
  if (!all(horInd >= 1L & horInd <= nc)) 
    stop("invalid argument 'horInd'")
  if (!all(verInd >= 1L & verInd <= nc)) 
    stop("invalid argument 'verInd'")
  if (doText) {
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      doText <- FALSE
  }
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  main <- if ("main" %in% nmdots) 
    dots$main
  if (is.null(oma)) 
    oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
  opar <- par(mfrow = c(length(horInd), length(verInd)), mar = rep.int(gap/2, 
                                                                       4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  xl <- yl <- logical(nc)
  if (is.numeric(log)) 
    xl[log] <- yl[log] <- TRUE
  else {
    xl[] <- grepl("x", log)
    yl[] <- grepl("y", log)
  }
  for (i in if (row1attop) 
    verInd
    else rev(verInd)) for (j in horInd) {
      l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", 
                                                 ""))
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ..., log = l)
      xyl = sapply(1:ncol(x), function(i) -c(alr(x, ivar = i)))
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        box()
        if (i == 1 && (!(j%%2L) || !has.upper || !has.lower)) 
          localAxis(1L + 2L * row1attop, xyl[, j], xyl[, i], 
                    ...)
        if (i == nc && (j%%2L || !has.upper || !has.lower)) 
          localAxis(3L - 2L * row1attop, xyl[, j], xyl[, i], 
                    ...)
        if (j == 1 && (!(i%%2L) || !has.upper || !has.lower)) 
          localAxis(2L, xyl[, j], xyl[, i], ...)
        if (j == nc && (i%%2L || !has.upper || !has.lower)) 
          localAxis(4L, xyl[, j], xyl[, i], ...)
        mfg <- par("mfg")
        if (i == j) {
          if (has.diag) 
            localDiagPanel(as.vector(x[, i]), ...)
          if (doText) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            xlp <- if (xl[i]) 
              10^0.5
            else 0.5
            ylp <- if (yl[j]) 
              10^label.pos
            else label.pos
            text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                       font = font.labels)
          }
        }
        else if (i < j) 
          localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                         i]), ...)
        else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                            i]), ...)
        if (any(par("mfg") != mfg)) 
          stop("the 'panel' function made a new plot")
      }
      else par(new = FALSE)
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
          font = font.main)
  }
  invisible(NULL)
}


#### compositional panel functions ---------------

vp.lrdensityplot = function (x, y, col=2,..., alpha = NULL){
  yh = hist(log(y/x), plot=FALSE)
  yd = density(log(y/x))
  usr <- par("usr")
  usr[1:2] <- range(yd$x)
  usr[3:4] <- range(yd$y)
  par(usr = usr)
  reject <- FALSE
  if (!is.null(alpha) && shapiro.test(y)$p < alpha) 
    reject <- TRUE
  plot(yh, main="", freq=F, add=T, col=ifelse(reject,NA,col))
  lines(yd, main="", xlab="", ylab="", xaxt="n", yaxt="n", col=ifelse(reject,col,1), lwd=2)
}


vp.lrboxplot = function (x, y, ... ){
  lxy = log(y/x)
  usr <- par("usr")
  #if(is.null(z)){
    boxplot(lxy, add=T,  ...)
  #}else{
  #  boxplot(lxy, add=T, ...)
  #}
}


vp.diffdensityplot <- function (x, y, col=2,..., alpha = NULL){
  yh = hist(y-x, plot=FALSE)
  yd = density(y-x)
  usr <- par("usr")
  usr[1:2] <- range(yd$x)
  usr[3:4] <- range(yd$y)
  par(usr = usr)
  if (!is.null(alpha) && is.factor(x)) 
    alpha <- alpha/nlevels(x)
  reject <- FALSE
  if (!is.null(alpha) && shapiro.test(y)$p < alpha) 
    reject <- TRUE
  plot(yh, main="", freq=F, add=T, col=ifelse(reject,NA,col))
  lines(yd, main="", xlab="", ylab="", xaxt="n", yaxt="n", col=ifelse(reject,col,1), lwd=2)
}



vp.kde2dplot <- function (x, y, grid = TRUE, legpos = "bottomright", colpalette = heat.colors, ...) {
  args = list(...)
  aux = MASS::kde2d(x, y, n = 50)
  aux$z = sqrt(aux$z)
  bks = hist(aux$z, plot = F, breaks = 20)$breaks
  cols = c(NA, colpalette(length(bks) - 2))
  args$add = c(args$add, TRUE)[1]
  args$xlab = c(args$xlab, "")[1]
  args$ylab = c(args$ylab, "")[1]
  args$breaks = bks
  args$col = cols
  args$x = aux
  do.call("image", args)
  xgrid = seq(from = floor(min(x)), to = ceiling(max(x)), by = 1)
  ygrid = seq(from = floor(min(y)), to = ceiling(max(y)), by = 1)
  abline(lm(y ~ x), col = 2, lwd = 2)
  if (grid) 
    abline(v = xgrid, h = ygrid, col = "#999999")
  legend(legpos, legend = round(cor(x, y), digits = 3), bg = "#999999")
}
