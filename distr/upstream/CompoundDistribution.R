
setMethod("NumbOfSummandsDistr", signature="CompoundDistribution",
           function(object) object@NumbOfSummandsDistr)
setMethod("SummandsDistr", signature="CompoundDistribution",
           function(object) object@SummandsDistr)

CompoundDistribution<- function( NumbOfSummandsDistr, SummandsDistr, .withSim = FALSE,
                                 withSimplify = FALSE){

  Symmetry <- NoSymmetry()

  if(!is(NumbOfSummandsDistr,"DiscreteDistribution"))
    stop("Argument 'NumbOfSummandsDistr' must be of class 'DiscreteDistribution'")

  supp <- support(NumbOfSummandsDistr)
  if(!(all(.isInteger(supp))&&all(supp >=0)))
    stop("Support of 'NumbOfSummandsDistr' must be non neg. integers")

  if(!is(SummandsDistr,"UnivDistrListOrDistribution"))
    stop("Argument 'SummandsDistr' must be of class 'UnivDistrListOrDistribution'")

##20200918 can be deleted:  supp <- support(NumbOfSummandsDistr)

  supp <- as(supp,"integer")
  suppNot0 <- supp[supp!=0L]

  ## new 20200918 triggered by mail by Vlada Milchevskaya vmilchev@uni-koeln.de
  ## special treatment of case support is of length 1
  if(length(supp)==1L){
     if(supp[1]==0L) return(Dirac(0))
     return(convpow(SummandsDistr,supp[1]))
  }

  is0 <- 0 %in% supp
  lI <- vector("list", length(supp))
  if(is0) lI[[1]] <- Dirac(0)
  ##  bugfix :: bug detected by Wolfgang Kreitmeier <wkreitmeier@gmx.de> 29.07.2016
  if(length(suppNot0)){
     if(is(SummandsDistr,"UnivariateDistribution")){
#        dsuppNot0 <- c(suppNot0,diff(suppNot0))
#        S <- 0
        for (i in 1:length(suppNot0)){
#             x0 <- convpow(SummandsDistr,suppNot0[i])
             S <- convpow(SummandsDistr,suppNot0[i])
#             S <- S + x0
             lI[[i+is0]] <- S
        }
      Symmetry <- Symmetry(SummandsDistr)
     }else{
       supp <- min(supp):max(supp)
       if( (length(supp)!=length(SummandsDistr)) &&
           !(is0 && length(supp)==1+length(SummandsDistr)))
          stop("Lengths of support of 'NumbOfSummandDistr' and list in 'SummandDistr' do not match")
       if(is0 && length(supp)==length(SummandsDistr))
          SummandsDistr <- SummandsDistr[2:length(SummandsDistr)]
       S <- 0
       Symm1 <- Symmetry(SummandsDistr[[1]])
       SymmL <- is(Symm1, "SphericalSymmetry")
       SymmC <- if(SymmL) SymmCenter(Symm1) else NULL
       for(i in 1:(length(supp)-is0)){
           if(SymmL && i>1){
              SymmI <- Symmetry(SummandsDistr[[i]])
              SymmL <- is(SymmI, "SphericalSymmetry")
              if(SymmL)
                 SymmL <- .isEqual(SymmCenter(SymmI),SymmC)
           }
           S <- S + SummandsDistr[[i]]
           lI[[i+is0]] <- S
       }
       if(SymmL) Symmetry <- SphericalSymmetry(SymmC)
     }
  UV <- do.call("UnivarMixingDistribution",
                 args = c(list(mixCoeff = d(NumbOfSummandsDistr)(supp),
                               withSimplify = FALSE),
                               lI)
                 )
  obj <- new("CompoundDistribution",
              NumbOfSummandsDistr = NumbOfSummandsDistr,
              SummandsDistr = SummandsDistr,
              p = UV@p, r = UV@r, d = UV@d, q = UV@q,
              mixCoeff = UV@mixCoeff, mixDistr = UV@mixDistr,
              .withSim = .withSim, .withArith = TRUE,
              Symmetry = Symmetry)

   if(withSimplify) return(simplifyD(obj))
   else return(obj)

  }
}

setMethod("+", c("CompoundDistribution","numeric"),
          function(e1, e2) simplifyD(e1)+e2)
setMethod("*", c("CompoundDistribution","numeric"),
          function(e1, e2) simplifyD(e1)*e2)


setMethod("Math", "AcDcLcDistribution",
          function(x){
            callGeneric(simplifyD(.ULC.cast(x)))
          })

setAs("CompoundDistribution", "UnivarLebDecDistribution",
       function(from)simplifyD(from))







####################################
