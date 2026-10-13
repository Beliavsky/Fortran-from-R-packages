# Genuine monobin 0.2.4 exported functions, not reimplemented R algorithms.
suppressPackageStartupMessages(library(monobin))
stopifnot(as.character(packageVersion('monobin')) == '0.2.4')
# Hmisc 5.3-0 was used for the checked-in reference fixtures.

outdir <- 'test/fixtures'
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
fmt <- function(x) ifelse(is.na(x), 'NaN', ifelse(is.infinite(x),ifelse(x>0,'Inf','-Inf'),sprintf('%.17g',x)))
methods <- c('pct','iso','woe','cum','sts','ndr','mdt')
# Fixed construction exercises repeated predictors, oscillation, differing bin
# prevalences, min-count corrections, rare events, and both target families.
n <- 600L
i <- seq_len(n)
x0 <- rep(seq_len(30),each=20)
y0 <- as.numeric(((i*17)%%100)/100 < (.06+.025*x0))
yc <- .1+.021*x0 + ((i*13)%%31)/100
sets <- list(
 binary=list(x=x0,y=y0),
 continuous=list(x=x0,y=yc,target=2L),
 decreasing=list(x=x0,y=rev(y0),direction=-1L),
 increasing=list(x=x0,y=y0,direction=1L),
 rare=list(x=x0,y=as.numeric((i%%83==0)&(x0>8))),
 options=list(x=x0,y=y0,g=9L,minobs=.12,minrate=.1,gap=.65,p=.2),
 specials=list(x=c(NaN,-999,Inf,-Inf,x0,5),y=c(0,1,0,1,y0,NaN),special=-999),
 separate=list(x=c(-999,Inf,-Inf,NaN,x0),y=c(1,0,1,0,y0),special=-999,separate=TRUE),
 noisy=list(x=(i*19)%%71,y=as.numeric((i*7)%%11<3)),
 zero_cells=list(x=x0,y=as.numeric(x0>15)),
 negative_cont=list(x=x0,y=-yc,target=2L,woe=FALSE),
 reverse_cont=list(x=x0,y=rev(yc),target=2L,direction=-1L),
 forced_opposite=list(x=x0,y=y0,direction=-1L),
 separated=list(x=x0,y=as.numeric(((i*17)%%100)/100 < pmin(1,.2+.04*x0))),
 negative_woe=list(x=x0,y=-yc,target=2L,woe=TRUE),
 mixed_cont=list(x=x0,y=yc-.5,target=2L),
 scarce_predictors=list(x=rep(c(1,2,1000),each=200),y=y0,g=60L),
 dense_groups=list(x=x0,y=y0,g=1000L)

)
records <- list(); failures <- list(); idx <- 0L
for (nm in names(sets)) for (method in methods) {
 d <- sets[[nm]]
 val <- function(name,default) if(is.null(d[[name]])) default else d[[name]]
 target<-val('target',0L); direction<-val('direction',0L); g<-val('g',if(method=='mdt')50L else 15L)
 separate<-val('separate',FALSE); woe<-val('woe',TRUE)
 minobs<-val('minobs',.05); minrate<-val('minrate',.01); gap<-val('gap',.1); p<-val('p',.05)
 sp<-val('special',numeric())
 args<-list(x=d$x,y=d$y,sc=c(NA,NaN,Inf,-Inf,sp),sc.method=if(separate)'separately' else 'together',
 y.type=if(target==0)NA else if(target==1)'bina' else 'cont', force.trend=if(direction==0)NA else if(direction==1)'i' else 'd')
 if(method %in% c('pct','cum','mdt'))args$g<-g
 if(method=='pct') args$woe.trend<-woe
 if(method %in% c('iso','woe','sts','ndr','mdt')) {args$min.pct.obs<-minobs;args$min.avg.rate<-minrate}
 if(method=='woe')args$woe.gap<-gap
 if(method %in% c('sts','ndr'))args$p.val<-p
 ans<-tryCatch(suppressWarnings(do.call(getExportedValue('monobin',paste0(method,'.bin')),args)),error=function(e)e)
 if(inherits(ans,'error') || !is.list(ans) || is.null(ans$summary.tbl)) {
  failures[[paste(nm,method)]]<-if(inherits(ans,'error'))conditionMessage(ans) else 'upstream diagnostic return'
  next
 }
 tab<-ans$summary.tbl
 if(anyNA(tab$type)) {failures[[paste(nm,method)]]<-'upstream missing type column on early cumulative return';next}
 # Fortran special groups are first-seen rather than formatted-label collation.
 special_rows<-which(tab$type=='special cases')
 if(length(special_rows)>1){
  first_seen<-unique(as.character(d$x[!is.na(d$y)&(d$x%in%args$sc)]))
  special_rows<-special_rows[match(first_seen,as.character(tab$bin[special_rows]))]
 }
 tab<-tab[c(special_rows,which(tab$type=='complete cases')),]
 member<-match(as.character(ans$x.trans),as.character(tab$bin))
 if(anyNA(member)) {failures[[paste(nm,method)]]<-'upstream x.trans contains unmatched labels';next}
 fullmember<-rep(0L,length(d$x));fullmember[!is.na(d$y)]<-member
 idx<-idx+1L; filename<-sprintf('%03d_%s_%s.dat',idx,nm,method)
 con<-file(file.path(outdir,filename),'w')
 writeLines(paste(length(d$x),method,target,direction,g,as.integer(separate),as.integer(woe),
 fmt(minobs),fmt(minrate),fmt(gap),fmt(p),length(sp)),con)
 if(length(sp))writeLines(paste(fmt(sp),collapse=' '),con)
 writeLines(paste(fmt(d$x),fmt(d$y)),con)
 writeLines(paste(nrow(tab),length(special_rows)),con)
 p.ref<-if('p.val'%in%names(tab))tab$p.val else rep(NA_real_,nrow(tab))
 for(j in seq_len(nrow(tab)))writeLines(paste(tab$no[j],paste(fmt(c(unlist(tab[j,c('y.sum','y.avg','x.avg','x.min','x.max','woe','iv.b')]),p.ref[j])),collapse=' ')),con)
 writeLines(paste(fullmember,collapse=' '),con)
 close(con)
 records[[idx]]<-c(nm,method,filename,nrow(tab))
}
writeLines(vapply(records,`[`,'',3),file.path(outdir,'index.txt'))
writeLines(c(paste('Fixtures:',length(records)),capture.output(sessionInfo()),'Upstream failures:',capture.output(failures)),
 'test/REFERENCE-SESSION.txt')
print(do.call(rbind,records));print(failures)
