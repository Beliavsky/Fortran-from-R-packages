# Recreate deterministic pinned posterior 1.7.0 references, from package root.
stopifnot(as.character(packageVersion("posterior")) == "1.7.0")
suppressPackageStartupMessages(library(posterior))
number <- function(x) ifelse(is.na(x),"NaN",ifelse(is.infinite(x),ifelse(x>0,"Inf","-Inf"),sprintf("%.17g",x)))
con <- file("test/reference_pareto.dat","w")
put <- function(x) writeLines(paste(number(x),collapse=" "),con)
p <- ((1:100)-.5)/100
set.seed(271828)
cases <- list(bounded=qgeneralized_pareto(p,k=-.4),exponential=qgeneralized_pareto(p,k=0),
              heavy=qgeneralized_pareto(p,k=.65),gaussian=rnorm(100),
              ties=round(qgeneralized_pareto(p,k=.2),1),unsorted=sample(qgeneralized_pareto(p,k=.3)))
put(length(cases))
for(nm in names(cases)) {
 x <- cases[[nm]]
 n <- length(x)
 weights <- seq(.1,2,length.out=n)
 if (nm == "unsorted") weights[seq(1,n,2)] <- 0
 excess <- x-min(x)+.01
 fit <- gpdfit(excess)
 fit_no_prior <- gpdfit(excess,wip=FALSE)
 fit_weighted <- gpdfit(excess,weights=weights)
 result <- suppressWarnings(lapply(c("both","left","right"),function(tail) pareto_smooth(x,tail=tail,r_eff=.73,return_k=TRUE,verbose=FALSE)))
 result[[4]] <- suppressWarnings(pareto_smooth(x,r_eff=.73,return_k=TRUE,verbose=FALSE,are_log_weights=TRUE))
 diags <- suppressWarnings(pareto_diags(x,r_eff=.73))
 writeLines(nm,con)
 put(n)
 put(x)
 put(weights)
 put(c(fit$k,fit$sigma,fit_no_prior$k,fit_no_prior$sigma,fit_weighted$k,fit_weighted$sigma,
       sapply(result,function(r)r$diagnostics$khat),unlist(diags)))
 for(r in result) put(r$x)
 # PIT comparisons have no ties, testing left, central and right GPD refinement.
 targets <- c(min(x)-.1, median(x)+.001, max(x)+.1)
 put(targets)
 put(sapply(targets,function(y) pit(matrix(x,ncol=1),y)))
 put(sapply(targets,function(y) pit(matrix(x,ncol=1),y,weights=matrix(weights,ncol=1))))
 put(sapply(targets,function(y) pareto_pit(matrix(x,ncol=1),y)))
 put(sapply(targets,function(y) pareto_pit(matrix(x,ncol=1),y,weights=matrix(log(weights),ncol=1),log=TRUE)))
}
close(con)
con <- file("test/reference_uniformity.dat","w")
cases <- list(c(.001,.025,.15,.19,.25,.43,.55,.58,.66,.81,.91,.94,.993),
              c(.01,.02,.03,.04,.07,.16,.26,.31,.38,.4,.44,.48,.6),
              c(.04,.04,.2,.2,.4,.4,.6,.6,.8,.8,.97,.97))
put(length(cases)*3)
for(x in cases) for(method in c("POT","PRIT","PIET")) {
 result <- uniformity_test(x,method)
 writeLines(method,con)
 put(length(x))
 put(x)
 put(c(result$pvalue,result$pointwise))
}
close(con)
