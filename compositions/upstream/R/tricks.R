#### wrappers around standard functions ----------

# test:
# data(Hydrochem)
# X1 = Hydrochem$H
# Y = acomp(Hydrochem[, 7:12])
# model = lm(alr(Y)~log(X1))
# anova(model)
anova = function(...){
  o = getStickyClassOption()
  setStickyClassOption(FALSE)
  on.exit(setStickyClassOption(o))
  #o = options("compositions")$compositions
  #fun = o$compositions$wrappedFunctions$anova
  res = stats::anova(...)
  return(res)
}


#summary.manova = function(object, ...){
#  for(nm in names(which(sapply(object, class)=="rmult"))){
#    object[[nm]] = unclass(object[[nm]])
#  }
#  stats::manova(object, ...)
#}


as.matrix.rmult = function(x, ...) oneOrDataset(unclass(x), B=x)


#### subsetting ----------

getStickyClassOption <- function() return(getOption("compositions")$stickyClass)

setStickyClassOption <- function(value){
  if(!is.logical(value) | is.na(value)) stop("setStickyClassOption: value must be logical")
  o = options("compositions")
  o$stickyClass <- value
  options(compositions=o)
}

gsi.subsetrow <- function(x, i, drop=FALSE){
  if(missing(i)) return(x)
  .stickyClass = getStickyClassOption()
  if(is.null(.stickyClass)) .stickyClass=FALSE
  y = unclass(x)[i,,drop=drop]
  if(.stickyClass) return(gsi.mystructure(y,class=class(x)))
  return(y)
}

gsi.LengthOne <- function(j){
  #!missing(j) && is.null(j) && length(j)==1
  if(missing(j)) return(FALSE)
  if(all(is.null(j))) return(FALSE)
  length(j)==1
}


"[.acomp" <- "[.rcomp"  <- 
  function(x, i, j, drop=gsi.LengthOne(j)){
  if(length(dim(x))==0) return(gsi.subsetvector(x,i))
  y = gsi.subsetrow(x,i, drop=missing(j))
  if(missing(j)){
    return(y)
  }else{
    y = unclass(y)[,j,drop=drop]
  }
  if(drop) return(y)
  if(ncol(y)>1) 
    return(gsi.mystructure(clo(y), class=class(x)))
  return(y)
}


"[.ccomp" <- "[.aplus" <- "[.rplus" <- 
  function(x, i, j, drop=gsi.LengthOne(j)){
  if(length(dim(x))==0) return(gsi.subsetvector(x,i))
  y = gsi.subsetrow(x,i, drop=missing(j))
  if(missing(j)){
    return(y)
  }else{
    y = unclass(y)[,j,drop=drop]
  }
  if(drop) return(y)
  return(gsi.mystructure(y, class=class(x)))
}


"[.rmult" <- function(x, i, j, drop=gsi.LengthOne(j)){
  if(length(dim(x))==0) return(gsi.subsetvector(x,i))
  .orig = gsi.orig(x)
  .V = gsi.getV(x)
  y = gsi.subsetrow(x,i, drop=missing(j))
  if(!is.null(.orig) & !missing(i)).orig=.orig[i,]
  if(!missing(j)){
    y = unclass(y)[,j,drop=drop]
    if(!is.null(.V)) .V = .V[,j]
  }
  .stickyClass = getStickyClassOption()
  if(is.null(.stickyClass) | drop) .stickyClass=FALSE
  if(.stickyClass) return(rmult(y, orig=.orig, V=.V))
  y
}


gsi.subsetvector <- function(x,i){
  if(missing(i)) return(x)
  .stickyClass = getStickyClassOption()
  if(is.null(.stickyClass)) .stickyClass=FALSE
  if(length(i)==1) .stickyClass=FALSE
  y = unclass(x)[i]
  if(.stickyClass & is.rmult(x)){
     .orig = gsi.orig(x)
     .V = gsi.getV(x)
      if(!is.null(.V) & ncol(.V)>1) .V = .V[,i] # if ncol(V)==1 we are in a (2,1) situation, and x is actually a column vector
     return(rmult(y, orig=.orig, V=.V))
  }
  if(.stickyClass & !(is.acomp(x)|is.rcomp(x)) ) return(gsi.mystructure(y,class=class(x)))
  if(.stickyClass & length(y)>1 ) return(gsi.mystructure(clo(y),class=class(x)))
  return(y)
}





"$.rmult" <- "$.rcomp" <- "$.rplus" <- "$.ccomp" <- "$.acomp" <- "$.aplus" <- function(x,name){
  if(name %in% colnames(x))
    return( unclass(x)[,name] )
  if(name=="pwlr") return(pwlr(x))
  if(name=="pwlrInv") return(pwlrInv(x))
  if(name=="alr") return(alr(x))
  if(name=="alrInv") return(alrInv(x))
  if(name=="clr") return(clr(x))
  if(name=="clrInv") return(clrInv(x))
  if(name=="ilr") return(ilr(x))
  if(name=="ilrInv") return(ilrInv(x))
  if(name=="apt") return(apt(x))
  if(name=="aptInv") return(aptInv(x))
  if(name=="cpt") return(cpt(x))
  if(name=="cptInv") return(cptInv(x))
  if(name=="ipt") return(ipt(x))
  if(name=="iptInv") return(iptInv(x))
  if(name=="cdt") return(cdt(x))
  if(name=="cdtInv") return(cdtInv(x))
  if(name=="idt") return(idt(x))
  if(name=="idtInv") return(idtInv(x))
  if(name=="iit") return(iit(x))
  if(name=="iitInv") return(iitInv(x))
  if(name=="ilt") return(ilt(x))
  if(name=="iltInv") return(iltInv(x))
  stop("[.xxx: 'name' given neither a column name nor one of the known data representations")
}
  

# "$.rmult" <- "$.rcomp" <- "$.rplus" <- "$.ccomp" <- "$.acomp" <- "$.aplus" <- function(x,name){
#   if(name %in% colnames(x))
#     return( unclass(x)[,name] )
#   if(name=="unclass") return(unclass(x))
#   if(name=="raw") return(unclass(x))
#   if(name=="clo") return(clo(x))
#   if(name=="pwlr") return(pwlr(x))
#   if(name=="pwlrInv") return(pwlrInv(x))
#   if(name=="alr") return(alr(x))
#   if(name=="alrInv") return(alrInv(x))
#   if(name=="clr") return(clr(x))
#   if(name=="clrInv") return(clrInv(x))
#   if(name=="ilr") return(ilr(x))
#   if(name=="ilrInv") return(ilrInv(x))
#   if(name=="apt") return(apt(x))
#   if(name=="aptInv") return(aptInv(x))
#   if(name=="cpt") return(cpt(x))
#   if(name=="cptInv") return(cptInv(x))
#   if(name=="ipt") return(ipt(x))
#   if(name=="iptInv") return(iptInv(x))
#   if(name=="cdt") return(cdt(x))
#   if(name=="cdtInv") return(cdtInv(x))
#   if(name=="idt") return(idt(x))
#   if(name=="idtInv") return(idtInv(x))
#   if(name=="iit") return(iit(x))
#   if(name=="iitInv") return(iitInv(x))
#   if(name=="ilt") return(ilt(x))
#   if(name=="iltInv") return(iltInv(x))
#   
# }




#### S4 behavour and inheritance as data.frame ---------
gsi.registerS4compositions <- function(from){
  to = as.data.frame(from)
  attr(to, "origClass")=class(from)
  return(to)
}

gsi.asMatrix <- function(from) as.matrix(x=gsi.registerS4compositions(from))

setOldClass(c("aplus", "data.frame"))
setAs("aplus", "data.frame", def=gsi.registerS4compositions)
setAs("aplus", "structure", def=gsi.asMatrix)

setOldClass(c("rplus", "data.frame"))
setAs("rplus", "data.frame", def=gsi.registerS4compositions)
setAs("rplus", "structure", def=gsi.asMatrix)

setOldClass(c("acomp", "data.frame"))
setAs("acomp", "data.frame", def=gsi.registerS4compositions)
setAs("acomp", "structure", def=gsi.asMatrix)

setOldClass(c("rcomp", "data.frame"))
setAs("rcomp", "data.frame", def=gsi.registerS4compositions)
setAs("rcomp", "structure", def=gsi.asMatrix)

setOldClass(c("ccomp", "data.frame"))
setAs("ccomp", "data.frame", def=gsi.registerS4compositions)
setAs("ccomp", "structure", def=gsi.asMatrix)

setOldClass(c("rmult", "data.frame"))
setAs("rmult", "data.frame", def=function(from){
  to = as.data.frame(from)
  attr(to, "origClass")="rmult"
  attr(to, "orig")= gsi.orig(from)
  attr(to, "V")= gsi.getV(from)
  return(to)
})
setAs("rmult", "structure", def=gsi.asMatrix)

setClassUnion(name="compositional", members=c("acomp", "rcomp", "ccomp"))
setClassUnion(name="amounts", members=c("aplus", "rplus", "ccomp"))


cdt.data.frame <- function(x,...){
  cl = attr(x, "origClass")
  if(is.null(cl)) return(x)
  class(x) = cl
  cdt(x)
}
idt.data.frame <- function(x,...){
  cl = attr(x, "origClass")
  if(is.null(cl)) return(x)
  class(x) = cl
  idt(x)
}

cdtInv.data.frame <- function(x, orig=gsi.orig(x),...){
  cl = c(attr(orig, "origClass"), class(orig))[1] 
  if(is.null(cl)) return(x)
  class(orig) = cl
  if("data.frame" %in% cl) return(x)
  cdtInv(x, orig)
}    
idtInv.data.frame <- function(x, orig=gsi.orig(x),...){
  cl = c(attr(orig, "origClass"), class(orig))[1] 
  if(is.null(cl)) return(x)
  class(orig) = cl
  if("data.frame" %in% cl) return(x)
  idtInv(x, orig)
}    


gsi.ORsequentially <- function(...){
  ll = list(...)
  for(l in ll){
    if(length(l)>0)if(l) return(TRUE)
  }
  return(FALSE)
}

gsi.ANDsequentially <- function(...){
  ll = list(...)
  for(l in ll){
    if(length(l)==0) return(FALSE)
    if(!l) return(FALSE)
  }  
  return(TRUE)
}



gsi.mystructure <- function(x, ...){
  if(is.null(x)) return(x)
  else return(structure(x, ...))
}



