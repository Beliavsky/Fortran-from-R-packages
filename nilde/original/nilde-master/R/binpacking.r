
gen.bin.inner <- function(tab, Sums, NumVag, Train, remSacks, Cap, vagon, flag){
 #OptVag <- length(inputa)
 #TrainList <- vector("list", length(inputa))
  OptVag <- bin.globals$OptVag
      if (sum(tab) > 0) {
      if (flag == 0) {
        if ( Cap == 0 ) { 
           Train <- paste0(Train," ",vagon) 
           NumVag1 <- NumVag + 1
           vagon <- character(0) 
           Cap <- input.n
        } else {NumVag1 <- NumVag} 
     if ( Cap == input.n ) {
        if (NumVag1 > OptVag) {
             return()
        }
      }
    if(length(Sums)>0) { smalls <- Sums <= Cap } else smalls <- 0
    i <- 1
    while (sum(smalls) > 0 & i <= sum(smalls)) {
        Ind <- which(smalls)[i]
        c.sol <- as.vector(tab[, Ind, drop=F], mode="logical")
            remSacks1 <- remSacks[!c.sol]
            vagon1 <- paste0(vagon, ifelse(identical(vagon, character(0)), "", ","), paste0(dimnames(tab)[[1]][c.sol],collapse=","))
            Cap1 <- Cap - Sums[i]
            NonZero <- (colSums(tab[!c.sol, , drop=F]) > 0) 
            Complementary <- colSums(tab[c.sol, , drop=F]) == 0 
            tab1 <- tab[!c.sol, which(NonZero & Complementary), drop=F] 
           Sums1 <- colSums(tab1*input.a[remSacks1])
           gen.bin.inner(tab=tab1, Sums=Sums1, NumVag=ifelse(Cap1==0, (NumVag+1), NumVag), Train=Train, remSacks=remSacks1, Cap=Cap1, vagon=vagon1, flag=ifelse(Cap1==0, 0, 1))
        i <- i + 1
    } 
    } else {
    if(length(Sums)>0) { smalls <- Sums <= Cap } else smalls <- 0
    j <- 1
    while (sum(smalls) > 0 & j <= sum(smalls)) {
        Ind <- which(smalls)[j]
        c.sol <- as.vector(tab[, Ind, drop=F], mode="logical")
        remSacks1 <- remSacks[!c.sol]
        vagon1 <- paste0(vagon, ifelse(identical(vagon, character(0)), "", ","), paste0(dimnames(tab)[[1]][c.sol],collapse=","))
        Cap1 <- Cap - Sums[Ind]
        NonZero <- (colSums(tab[!c.sol, , drop=F]) > 0) 
        Complementary <- colSums(tab[c.sol, , drop=F]) == 0 
        tab1 <- tab[!c.sol, which(NonZero & Complementary), drop=F] 
        Sums1 <- colSums(tab1*input.a[remSacks1]) 
        gen.bin.inner(tab=tab1, Sums=Sums1, NumVag=ifelse(Cap1==0, (NumVag+1), NumVag), Train=Train, remSacks=remSacks1, Cap=Cap1, vagon=vagon1, flag=ifelse(Cap1==0, 0, 1))
        j <- j + 1
    } 
     if (sum(smalls) == 0) {
         Train1 <- paste0(Train," ",vagon) 
         vagon1 <- character(0)
         NumVag1 <- NumVag + 1
         Cap1 <- input.n
         gen.bin.inner(tab=tab, Sums=Sums, NumVag=NumVag1, Train=Train1, remSacks=remSacks, Cap=Cap1, vagon=vagon1, flag=0)
     } 
    }
} else {
if (!identical(vagon, character(0))) {
        Train1 <- paste0(Train," ",vagon) 
        NumVag1 <- NumVag + 1
      } else {Train1 <- Train}
        Train1 <- substr(Train1, 2,nchar(Train1))
        . <- gsub("s", "", Train1, fixed=T)
        . <-strsplit(., split=" ", fixed=T)[[1]] 
        . <-strsplit(., split=",", fixed=T)     
        . <- lapply(., as.numeric)             
        . <- lapply(., function(x) x[order(x)])
        vmin <- sapply(., min)                
        . <- .[order(vmin)]            
        Len <- length(.)              
        vagIneff <- sapply(., function(x) input.n - sum(input.a[x]))
        trIneff <- sum(vagIneff)          
        code <- sapply(., paste0, collapse=",")
        code <- paste0(code, collapse=" ")
        if (Len <= OptVag){
            #  OptVag <<- Len
             # assign("OptVag", Len, envir = .GlobalEnv)
              bin.globals$OptVag <- Len
            if(is.null(bin.globals$TrainList[[Len]])){
              #TrainList[[Len]]$code <<- c(TrainList[[Len]]$code, code)
              #TrainList[[Len]]$trIneff <<- c(TrainList[[Len]]$trIneff, trIneff)
              # TrainList[[Len]]$vagIneff <<- c(TrainList[[Len]]$vagIneff, vagIneff)   
              bin.globals$TrainList[[Len]]$code <- c(bin.globals$TrainList[[Len]]$code, code)
              bin.globals$TrainList[[Len]]$trIneff <- c(bin.globals$TrainList[[Len]]$trIneff, trIneff)
              bin.globals$TrainList[[Len]]$vagIneff <- c(bin.globals$TrainList[[Len]]$vagIneff, vagIneff)   
          } else {
            if(is.na(match(code, bin.globals$TrainList[[Len]]$code))) {
              #TrainList[[Len]]$code <<- c(TrainList[[Len]]$code, code)
              #TrainList[[Len]]$trIneff <<- c(TrainList[[Len]]$trIneff, trIneff)
              #TrainList[[Len]]$vagIneff <<- c(TrainList[[Len]]$vagIneff, vagIneff)
              bin.globals$TrainList[[Len]]$code <- c(bin.globals$TrainList[[Len]]$code, code)
              bin.globals$TrainList[[Len]]$trIneff <- c(bin.globals$TrainList[[Len]]$trIneff, trIneff)
              bin.globals$TrainList[[Len]]$vagIneff <- c(bin.globals$TrainList[[Len]]$vagIneff, vagIneff)
                       }       
              }
        }
}
return()
}


## wrapper function for gen.bin.inner()...

bin.packing <- function(input.a, input.n,bin.globals){
## input.a - weights of items
## input.n - capacity of a bin
  ## prepare inputs for gen.bin.innner function ...
   a2 <- nlde(a=input.a, n=input.n, option=2)
   tab <-  a2$solution[, -1] 
   Sums <- colSums(tab*input.a); names(Sums) <- NULL # row of efficiencies
   MaxSacks <- length(input.a)
   remSacks <- 1:length(input.a) # indexes of remaining sacks in initial vector
   Cap <- input.n
   NumVag <- 0
 # bin.globals <- new.env() # new environment for global variables
 # bin.globals$OptVag <- length(input.a) # Note: global, initial min # of vagons
  # bin.globals$TrainList <- vector("list", length(input.a) ) # Note: global
   Train <- character(0)
   vagon <- character(0)
   flag <- 1 # 0, curr bin is complete, need to start new one 
   g <- gen.bin.inner(tab=tab, Sums=Sums,NumVag=NumVag,Train=Train,remSacks=remSacks,
              Cap=Cap, vagon=vagon, flag=flag)
   #return(bin.globals$TrainList)
   min.bins <- NULL ## min number of bins required
   for (jj in 1:length(bin.globals$TrainList)){
        if (!is.null(bin.globals$TrainList[[jj]])){
           min.bins <- jj
           break
       }
   }
   list(solution=bin.globals$TrainList[[min.bins]]$code, min.bins=min.bins, bin.ineff=bin.globals$TrainList[[min.bins]]$vagIneff, 
        total.ineff=bin.globals$TrainList[[min.bins]]$trIneff)
}




