plot3D <- function(x,...) UseMethod("plot3D")

plot3D.rmult <- function(x,parts=1:3,...,center=FALSE,scale=FALSE,add=FALSE,axes=!add,cex=2,vlabs=colnames(x),size=cex,bbox=FALSE,col=1) {
  requireNamespace("rgl")
  ddd = TRUE
  X<-x
  if( ! add ) 
      gsi.reset3D()
    x <- X[,parts, drop=ddd]
    if( length(scale) == 1 ) {
      if( is.logical(scale) )
        if( scale )
          scale <- diag(c(1/sqrt(diag(var(x))),1))
        else
          scale <- diag(rep(1,4))
    } else {
      scale <- diag(scale,1)
    }
    scale[4,4] <- max(diag(scale)[-4]) 
    rgl::points3d(x[,1, drop=ddd],x[,2, drop=ddd],x[,3, drop=ddd],size=size,...,col=col)
    rgl::par3d(userMatrix=scale)
    gsi.filtercall("axis3D")
                                        #  axis3D(axis.origin=axis.origin,axis.scale=axis.scale,
                                        #         axis.col=axis.col,vlabs=vlabs,
                                        #         vlabs.col=vlabs.col,bbox=bbox,
                                        #         axis.lwd=axis.lwd,axis.len=axis.len,
                                        #         axis.angle=axis.angle,oth=orth,axes=axes)
    invisible(rmult(x))
  }

gsi.reset3D <- function(userMatrix=diag(rep(1,4))) {
  requireNamespace("rgl")
  rgl::clear3d()
  rgl::par3d(userMatrix=userMatrix)
}


gsi.filtercall <- function(fkt,...,overwrite=NULL,transmit=character(0),default=list(...),prefix=if(is.character(fkt)) fkt else NULL) {
  #oo <- function(x) {cat(deparse(substitute(x)),"=\n");print(x)}
  #oo <- function(x) {}
  requireNamespace("rgl")
  prefix
  if(is.character(fkt)) fkt <- get(fkt,mode="function")
  #oo(fkt)
  call=sys.call(sys.parent())
  #oo(call[[1]])
  f <- get(as.character(call[[1]]),mode="function")
  match <- as.list(match.call(f,call=call))[-1]
  #oo(match)
  nn <- names(match)
  nn <- nn[!is.null(nn)& !is.na(nn)]
  #oo(match[nn])
  args <- default
  #oo(args)
  args[nn] <- match[nn]
  #oo(args)
  #oo(default)
  if( !is.null(prefix) ) {
    e <- as.list(match)[grep(paste("^",prefix,"\\.",sep=""),names(match))]
    #oo(e)
    if( length(e) > 0 ) {
      nc <- nchar(prefix)
      names(e) <- substring(names(e),nc+2)
      #oo(e)
      args[names(e)]<-e           
    }
  }
  #oo(args)
  args[names(overwrite)]<-overwrite
  #oo(args)

  nn <- unique(c(as.character(names(formals(fkt))),transmit,names(overwrite),names(default)))
  #oo(nn)
  nn <- nn[!is.null(nn)&!is.na(nn)]
  nn <- nn[nn!="..." & nn!=""]
  nn <- nn[nn %in% names(args)]
  #oo(nn)
  args <- args[nn]
  #oo(args)
  do.call(fkt,args[nn],envir=parent.frame(2)) # this line is the source of error!!!
}

plot3D.default <- function(x,...,add=FALSE,bbox=TRUE,axes=FALSE,cex=1,size=cex,col=1) {
  requireNamespace("rgl")
  ddd = TRUE
  X<-x
  if( ! add )
    gsi.reset3D()
  radius <- max(max(X[,1, drop=ddd])-min(X[,1, drop=ddd]),max(X[,2, drop=ddd])-min(X[,2, drop=ddd]),max(X[,3, drop=ddd])-min(X[,3, drop=ddd]))/100*cex
  rgl::spheres3d(X[,1, drop=ddd],X[,2, drop=ddd],X[,3, drop=ddd],...,radius=radius,col=col)
  if( bbox) rgl::bbox3d(color = c("white", "black"))
  if( any(axes) ) arrows3D(diag(c(0,0,0)),diag(c(1,1,1))) 
  invisible(X[,1:3, drop=ddd])
}



plot3D.acomp <- function(x,parts=1:min(ncol(X),4),...,lwd=2,axis.col="gray",add=FALSE,cex=2,vlabs=colnames(x),vlabs.col=axis.col,center=FALSE,scale=FALSE,log=FALSE,bbox=FALSE,axes=TRUE,size=cex,col=1) {
  requireNamespace("rgl")
  ddd = TRUE
  X<-x
  out = NULL
  if( length(parts) == 3 ) {
   if( log ) {
     x <- clr(scale(acomp(X,parts=parts),center=center,scale=scale))
     if( ! add ) {
       gsi.reset3D()
       if( axes )
         arrows3D(diag(c(0,0,0)),diag(c(1,1,1)),labs=vlabs,col=axis.col)
     }
      rgl::points3d(x[,1, drop=ddd],x[,2, drop=ddd],x[,3, drop=ddd],size=size,...,col=col)
      out = rmult(x[,1:3, drop=ddd])
    } else {
      x <- scale(acomp(X,parts=parts),center=center,scale=scale)
      if( ! add ) {
        gsi.reset3D()
        corners <- rbind(diag(rep(1,3)),c(0,0,0))
        cl <- corners[c(1,2,3,4,1,3,2,4),]
        if( axes ) 
          rgl::lines3d(cl[,1],cl[,2],cl[,3],col=axis.col,size=lwd)
        if( !is.null(vlabs) )
          rgl::texts3d(corners[,1],corners[,2],corners[,3],c(vlabs,"0"),col=vlabs.col)
      }
      rgl::points3d(x[,1, drop=ddd],x[,2, drop=ddd],x[,3, drop=ddd],size=size,...,col=col)
      out = rmult(x[,1:3, drop=ddd])
    }
   rgl::view3d(45,35.4)
 } else if( length(parts)==4 ) {
   x <- clo(X,parts=parts)
   if( log ) {
     if( ! add ) {
       gsi.reset3D()
       corners <- normalize(ilr(diag(rep(0.5,4))+0.1))
       if( axes )
         arrows3D(corners*0,corners,col=axis.col,size=lwd,labs=vlabs)
     }
     ilrx <- ilr(scale(acomp(x),center=center,scale=scale))
     rgl::points3d(ilrx[,1, drop=ddd],ilrx[,2, drop=ddd],ilrx[,3, drop=ddd],size=size,...,col=col)
     out = rmult(ilrx[,1:3, drop=ddd])
   } else {
     if( ! add ) {
       gsi.reset3D()
       corners <- diag(rep(1,4))
       cornerlines <- corners[c(1,2,3,4,1,3,2,4),]
       cl <- ipt(cornerlines)
       if( axes )
         rgl::lines3d(cl[,1, drop=ddd],cl[,2, drop=ddd],cl[,3, drop=ddd],col=axis.col,size=lwd)
     }
     iptx <- ipt(scale(acomp(x),center=center,scale=scale))
     rgl::points3d(iptx[,1, drop=ddd],iptx[,2, drop=ddd],iptx[,3, drop=ddd],size=size,...,col=col)
     out = rmult(iptx[,1:3, drop=ddd])
     if( !is.null(vlabs) ) {
       cc <- ipt(corners)
       rgl::texts3d(cc[,1, drop=ddd],cc[,2, drop=ddd],cc[,3, drop=ddd],c(vlabs),col=vlabs.col)
     }
   }
 } else
  stop("Wrong number of parts")
  if( bbox )
    rgl::bbox3d()
  invisible(out)
}
  
plot3D.rcomp <- plot3D.acomp



axis3D <- function(axis.origin=c(0,0,0),axis.scale=1,axis.col="gray",vlabs=c("x","y","z"),vlabs.col=axis.col,bbox=FALSE,axis.lwd=2,axis.len=mean(axis.scale)/10,axis.angle=30,orth=c(1,0.0001,0.000001),axes=TRUE,...) {
  requireNamespace("rgl")
  if( axes ) {
    M <- rbind(axis.origin,axis.origin,axis.origin)
    if( length(axis.scale)== 1) axis.scale <- rep(axis.scale,3)
    arrows3D(M,M+diag(axis.scale),length=axis.len,lwd=axis.lwd,
             angle=axis.angle,
             code=2,col=axis.col,orth=orth,labs=vlabs,...)
  }
  if( bbox ) {
    rgl::bbox3d(xlab=vlabs[1],ylab=vlabs[2],zlab=vlabs[3])
  }
}

plot3D.rplus <- function(x,parts=1:3,...,vlabs=NULL,add=FALSE,bbox=FALSE,cex=1,size=cex,axes=TRUE,col=1) {
  requireNamespace("rgl")
  X<-x
  if(! add )
    gsi.reset3D()
  XX <- oneOrDataset(iit(rplus(X,parts=parts)))
  out = plot3D(rmult(XX),...,bbox=bbox,col=col)
  if( missing(vlabs) )
    vlabs <- colnames(XX)
  if( axes ) arrows3D(diag(c(0,0,0)),diag(mean(rplus(XX))),labs=vlabs)
  invisible(out)
}

plot3D.aplus <- function(x,parts=1:3,...,vlabs=NULL,add=FALSE,log=TRUE,bbox=FALSE,axes=TRUE,col=1) {
  requireNamespace("rgl")
  X<-x
  if( ! log ) {
    out = plot3D(rplus(X),parts=parts,...,col=col,vlabs=vlabs,add=add,bbox=bbox)
    return()
  }
  if(! add )
    gsi.reset3D()
  XX <- oneOrDataset(ilt(aplus(X,parts=parts)))
  if( missing(vlabs) )
    vlabs <- colnames(XX)
  out = plot3D(rmult(XX),...,col=col,bbox=bbox)
  mm <- mean(rmult(XX))
  if( axes )
    arrows3D(diag(mm+c(0,0,0)),diag(mm+c(1,1,1)),labs=vlabs)
 # bbox3d()
  invisible(out)
}


arrows3D <- function(...) UseMethod("arrows3D")

#par("fg") and par("col") # when the plotting functions put their reference arrows, the call to par() generates a blank plot, replaced by col=1
arrows3D.default <- function(x0,x1,...,length=0.25,
                     angle=30,code=2,col="black",lty=NULL,lwd=2,orth=c(1,0.0001,0.0000001),labs=NULL,size=lwd) {
  requireNamespace("rgl")
  X <- oneOrDataset(x0,x1)
  Y <- oneOrDataset(x1,x0)
  if( ! is.null(labs))
    rgl::texts3d(Y[,1],Y[,2],Y[,3],labs,...,col=col)
  XX <- rbind(X,Y)
  XX[(1:nrow(X)-1)*2+1,]<- X
  XX[(1:nrow(X)-1)*2+2,]<- Y
  rgl::segments3d(XX[,1],XX[,2],XX[,3],col=col,size=size,...)
  out = rmult(XX[,1:3])
  if( code > 0 ) {
    XY <- rmult(X)-rmult(Y)
    XL <- norm(XY)
    XD <- XY / XL
    ZD <- normalize(rmult(cbind(XD[,2]*orth[3]-XD[,3]*orth[2],
                                XD[,3]*orth[1]-XD[,1]*orth[3],
                                XD[,1]*orth[2]-XD[,2]*orth[1])))
    s <- sin(angle/180*pi)
    c <- cos(angle/180*pi)
    if( code  %% 2 == 1 ) {
      XX[(1:nrow(X)-1)*2+1,]<- X
      XX[(1:nrow(X)-1)*2+2,]<- X - c*length*XY + ZD*(XL*length*s)
      rgl::segments3d(XX[,1],XX[,2],XX[,3],col=col,size=size,...)
      XX[(1:nrow(X)-1)*2+1,]<- X
      XX[(1:nrow(X)-1)*2+2,]<- X - c*length*XY - ZD*(XL*length*s)
      rgl::segments3d(XX[,1],XX[,2],XX[,3],col=col,size=size,...)
    }
    if( floor(code/2)  %% 2 == 1 ) {
      XX[(1:nrow(X)-1)*2+1,]<- Y
      XX[(1:nrow(X)-1)*2+2,]<- Y + c*length*XY + ZD*(XL*length*s)
      rgl::segments3d(XX[,1],XX[,2],XX[,3],col=col,size=size,...)
      XX[(1:nrow(X)-1)*2+1,]<- Y
      XX[(1:nrow(X)-1)*2+2,]<- Y + c*length*XY - ZD*(XL*length*s)
      rgl::segments3d(XX[,1],XX[,2],XX[,3],col=col,size=size,...)
    }
  }
  invisible(out)
}

biplot3D <- function(x,...) UseMethod("biplot3D")

biplot3D.default <- function(x,y,var.axes=TRUE,
                             col=c("green","red"),cex=c(2,2),
            xlabs = NULL, ylabs = NULL, expand = 1, arrow.len = 0.1,
#            xlim  = NULL, ylim  = NULL, 
#            main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                             ...,add=FALSE){
  requireNamespace("rgl")
  if( !add ) {
    gsi.reset3D();
  }
  x <- oneOrDataset(x)
  rgl::points3d(x[,1],x[,2],x[,3],size=cex[1],col=col[[1]],...)
  if( ! is.null(xlabs) ) {
    rgl::texts3d(x[,1],x[,2],x[,3],xlabs,...)
  }
  if( var.axes ) {
    y <- oneOrDataset(expand*y)
    arrows3D(c(0,0,0),y,length=arrow.len,
             lwd=cex[[1%%length(cex)+1]],
             col=col[[1%%length(col)+1]],...)
    
  } else {
    y <- oneOrDataset(expand*y)
    rgl::points3d(y[,1],y[,2],y[,3],
             size=cex[[1%%length(cex)+1]],
             col=col[[1%%length(col)+1]],...)
  }
  out = rmult(y[,1:3])
  if( ! is.null(ylabs) ) {
    rgl::texts3d(y[,1],y[,2],y[,3],ylabs,...)
  }
  invisible(out)
}

biplot3D.princomp <- function(x,choices=1:3,scale=1,...,comp.col=1,comp.labs=paste("Comp.",1:3),scale.scores=lambda[choices]^(1-scale),scale.var=scale.comp,scale.comp=sqrt(lambda[choices]),scale.disp=1/scale.comp) {
  requireNamespace("rgl")
  lambda <- x$sdev^2
  if( length(scale.var) == 1 ) scale.var <- rep(scale.var,3)
  if( length(scale.scores) == 1 ) scale.scores <-  rep(scale.scores,3)
  if( length(scale.comp) == 1 ) scale.comp <-  rep(scale.comp,3)
  if( length(scale.disp) == 1 ) scale.disp <-  rep(scale.disp,3)
#  biplot3D.default( unclass(x$scores[,choices]) %*% diag(lambda[choices]^(1-scale)),
#  unclass(x$loadings[,choices])%*% diag(lambda[choices]^scale),...)
  out = biplot3D.default( unclass(x$scores[,choices]) %*% diag(scale.scores),
  unclass(x$loadings[,choices])%*% diag(scale.var),...,ylabs=rownames(x$loadings),xlabs=rownames(x$scores))#
  comp <- diag(scale.comp)
  arrows3D(c(0,0,0),comp,col=comp.col)
  if( ! is.null(comp.labs) )
    rgl::texts3d(comp[,1],comp[,2],comp[,3],comp.labs,col=comp.col)
  rgl::par3d(userMatrix=diag(c(scale.disp/max(scale.disp),1)))
  invisible(out)
}
