
### should move checkMaxControl to a separate file but how to do it?

setClassUnion("functionOrNULL", c("function", "NULL"))
setClassUnion("integerOrNULL", c("integer", "NULL"))
setClassUnion("numericOrNULL", c("numeric", "NULL"))

checkMaxControl <- function(object) {
   ## check validity of MaxControl objects
   if(!inherits(object, "MaxControl")) {
      stop("'MaxControl' object required.  Currently '",
           paste(class(object), sep=", "), "'")
   }
   ##
   errors <- character(0)
   ## Check length of componenents
   for(s in slotNames(object)) {
      if(s == "sann_cand") {
         if(length(slot(object, s)) > 1) {
            errors <- c(errors,
                        paste("'", s, "' must be either 'NULL' or ",
                              "a function of length 1, not of length ",
                              length(slot(object, s)), sep=""))
         }
      }
      else if(s %in% c("SG_batchSize", "SG_clip", "SG_patience")) {
                           # integerOrNULL
         if(length(slot(object, s)) > 1) {
            errors <- c(errors,
                        paste("'", s, "' must be either 'NULL' or ",
                              "of length 1, not of length ",
                              length(slot(object, s)), sep=""))
         }
      }
      else if(length(slot(object, s)) != 1) {
                           # length 1
         errors <- c(errors,
                     paste("'", s, "' must be of length 1, not ",
                           length(slot(object, s)), sep=""))
      }
   }
   ## check missings
   for(s in slotNames(object)) {
      if(is.vector(slot(object, s)) && any(is.na(slot(object, s)))) {
                           # is.na only works for vectors
         errors <- c(errors,
                     paste0("NA in '", s, "'")
                     )
         return(errors)
                           # return errors here as otherwise NA-s will interfere the
                           # block of if-s below
      }
   }
   ##
   if(slot(object, "steptol") < 0) {
      errors <- c(errors, paste("'steptol' must be non-negative, not",
                                slot(object, "steptol")))
   }
   if(slot(object, "lambdatol") < 0) {
      errors <- c(errors, paste("'lambdatol' must be non-negative, not",
                                slot(object, "lambdatol")))
   }
   ## qac valid values--only check if length 1
   if(length(slot(object, "qac")) == 1 &&
      !pmatch(slot(object, "qac"), c("stephalving", "marquardt"))) {
      errors <- c(errors, paste("'qac' must be 'stephalving' or 'marquadt', not",
                                slot(object, "qac")))
   }
   if(slot(object, "qrtol") < 0) {
      errors <- c(errors, paste("'qrtol' must be non-negative, not",
                                slot(object, "qrtol")))
   }
   if(slot(object, "marquardt_lambda0") < 0) {
      errors <- c(errors, paste("'lambda0' must be non-negative, not",
                                slot(object, "lambda0")))
   }
   if(slot(object, "marquardt_lambdaStep") <= 1) {
      errors <- c(errors, paste("'lambdaStep' must be > 1, not",
                                slot(object, "lambdaStep")))
   }
   if(slot(object, "marquardt_maxLambda") < 0) {
      errors <- c(errors, paste("'maxLambda' must be non-negative, not",
                                slot(object, "maxLambda")))
   }
   ## NM
   if(slot(object, "nm_alpha") < 0) {
      errors <- c(errors, paste("Nelder-Mead reflection factor 'alpha' ",
                                "must be non-negative, not", slot(object, "nm_alpha")))
   }
   if(slot(object, "nm_beta") < 0) {
      errors <- c(errors, paste("Nelder-Mead contraction factor 'beta' ",
                                "must be non-negative, not", slot(object, "nm_beta")))
   }
   if(slot(object, "nm_gamma") < 0) {
      errors <- c(errors, paste("Nelder-Mead expansion factor 'gamma' ",
                                "must be non-negative, not", slot(object, "nm_gamma")))
   }
   ## SANN
   if(!inherits(slot(object, "sann_cand"), c("function", "NULL"))) { #
      errors <- c(errors, paste("'SANN_cand' must be either NULL or a function, not",
                                slot(object, "SANN_cand")))
   }
   if(slot(object, "sann_tmax") < 1) {
      errors <- c(errors, paste("SANN number of calculations at each temperature ",
                                "'tmax' ",
                                "must be positive, not", slot(object, "sann_tmax")))
   }
   ## SGA
   if(slot(object, "SGA_momentum") < 0 || slot(object, "SGA_momentum") > 1) {
      errors <- c(errors, paste("SGA momentum parameter must be in [0,1], not",
                                slot(object, "SGA_momentum")))
   }
   ## Adam
   if(slot(object, "Adam_momentum1") < 0 || slot(object, "Adam_momentum1") > 1) {
      errors <- c(errors, paste("Adam momentum1 parameter must be in [0,1], not",
                                slot(object, "Adam_momentum1")))
   }
   if(slot(object, "Adam_momentum2") < 0 || slot(object, "Adam_momentum2") > 1) {
      errors <- c(errors, paste("Adam momentum2 parameter must be in [0,1], not",
                                slot(object, "Adam_momentum2")))
   }
   ## SG general
   if(slot(object, "SG_learningRate") <= 0) {
      errors <- c(errors, paste("learning rate for SGA must be positive, not",
                                slot(object, "SG_learningRate")))
   }
   if(length(slot(object, "SG_batchSize")) > 0 && slot(object, "SG_batchSize") <= 0L) {
      errors <- c(errors, paste("SGA batch size must be positive, not",
                                slot(object, "SG_batchSize")))
   }
   if(length(slot(object, "SG_clip")) > 0 && slot(object, "SG_clip") <= 0L) {
      errors <- c(errors, paste("SGA gradient clip norm threshold must be positive, not",
                                slot(object, "SG_clip")))
   }
   if(length(slot(object, "SG_patience")) > 0 && slot(object, "SG_patience") <= 0L) {
      errors <- c(errors, paste("SG patience must be positive (or NULL), not",
                                slot(object, "SG_patience")))
   }
   if(slot(object, "SG_patienceStep") <= 0L) {
      errors <- c(errors, paste("SG patience step must be positive, not",
                                slot(object, "SG_patienceStep")))
   }
   ## general
   if(slot(object, "iterlim") < 0) {
      errors <- c(errors, paste("'iterlim' must be non-negative, not",
                                slot(object, "iterlim")))
   }
   if(slot(object, "max.rows") < 0) {
      errors <- c(errors, paste("'max.rows' must be non-negative, not",
                                slot(object, "max.rows")))
   }
   if(slot(object, "max.cols") < 0) {
      errors <- c(errors, paste("'max.cols' must be non-negative, not",
                                slot(object, "max.cols")))
   }
   if(length(errors) > 0)
      return(errors)
   return(TRUE)
}

### MaxControls contains all control parameters for max* family
setClass("MaxControl",
         slots=representation(
             tol="numeric",
             reltol="numeric",
             gradtol="numeric",
             steptol="numeric",
                           #
             lambdatol="numeric",
             qrtol="numeric",
             ## Qadratic Approximation Control
             qac="character",
         marquardt_lambda0="numeric",
         marquardt_lambdaStep="numeric",
         marquardt_maxLambda="numeric",
         ## Optim Nelder-Mead:
         nm_alpha="numeric",
         nm_beta="numeric",
         nm_gamma="numeric",
         ## SANN
         sann_cand="functionOrNULL",
         sann_temp="numeric",
         sann_tmax="integer",
         sann_randomSeed="integer",
         ## SGA
         SGA_momentum = "numeric",
         ## Adam
         Adam_momentum1 = "numeric",
         Adam_momentum2 = "numeric",
         ## SG general
         SG_patience = "integerOrNULL",  # NULL: don't care about patience
         SG_patienceStep = "integer",  # check patience at every epoch
         SG_learningRate="numeric",
         SG_batchSize = "integerOrNULL",  # NULL: full batch
         SG_clip="numericOrNULL",  # NULL: do not clip
         ##
         iterlim="integer",
         max.rows="integer",
         max.cols="integer",
         printLevel="integer",
         storeValues="logical", storeParameters="logical"
         ),
         ##
         prototype=prototype(
             tol=1e-8,
             reltol=sqrt(.Machine$double.eps),
             gradtol=1e-6,
             steptol=1e-10,
                           #
             lambdatol=1e-6,
                           #
             qac="stephalving",
             qrtol=1e-10,
            marquardt_lambda0=1e-2,
         marquardt_lambdaStep=2,
         marquardt_maxLambda=1e12,
         ## Optim Nelder-Mead
         nm_alpha=1,
         nm_beta=0.5,
         nm_gamma=2,
         ## SANN
         sann_cand=NULL,
         sann_temp=10,
         sann_tmax=10L,
         sann_randomSeed=123L,
         ## SGA
         SGA_momentum = 0,
         ## Adam
         Adam_momentum1 = 0.9,
         Adam_momentum2 = 0.999,
         ##
         SG_learningRate=0.1,
         SG_batchSize=NULL,
         SG_clip=NULL,
         SG_patience = NULL,
         SG_patienceStep = 1L,
         ##
         iterlim=150L,
         max.rows=as.integer(getOption("max.rows", 20L)),
         max.cols=as.integer(getOption("max.cols", 7L)),
         printLevel=0L,
         storeValues=FALSE, storeParameters=FALSE),
         ##
         validity=checkMaxControl
         )
