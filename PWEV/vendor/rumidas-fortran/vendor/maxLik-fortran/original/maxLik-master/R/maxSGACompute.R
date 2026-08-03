maxSGACompute <- function(fn, grad, hess,
                           start, 
                           nObs,
                           finalHessian=FALSE,
                           bhhhHessian = FALSE,
                           fixed=NULL,
                           control=maxControl(),
                           optimizer="SGA",  # type of optimizer: SGA, Adam
                           ...) {
   ## Stochastic Gradient Ascent: implements
   ## * SGA with momentum
   ## * Adam 
   ## Parameters:
   ## fn          - the function to be maximized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ## start       - initial parameter vector (eventually w/names)
   ## control       MaxControl object:
   ##     The stopping criteria
   ##     tol         - maximum allowed absolute difference between sequential values
   ##     reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ##     gradtol     - maximum allowed norm of gradient vector
   ## 
   ##     iterlim     - maximum # of iterations
   ##     
   ## finalHessian  include final Hessian?  As computing final hessian does not carry any extra penalty for NR method, this option is
   ##               mostly for compatibility reasons with other maxXXX functions.
   ##               TRUE/something else  include
   ##               FALSE                do not include
   ## fixed       - a logical vector -- which parameters are taken as fixed.
   ##               Other paramters are treated as variable (free).
   ## ...           additional argument to 'fn'.  This may include
   ##               'fnOrig', 'gradOrig', 'hessOrig' if called fromm
   ##               'maxNR'.
   ##
   ## RESULTS:
   ## an object of class 'maxim'
   ##      
   ## -------------------------------------------------
   maximType <- "Stochastic Gradient Ascent"
   iterlim <- slot(control, "iterlim")
   nParam <- length(start)
   start1 <- start
   storeParameters <- slot(control, "storeParameters")
   storeValues <- slot(control, "storeValues")
   learningRate <- slot(control, "SG_learningRate")
   clip <- slot(control, "SG_clip")
   max.rows <- slot(control, "max.rows")
   max.cols <- slot(control, "max.cols")
   patience <- slot(control, "SG_patience")
   patienceStep <- slot(control, "SG_patienceStep")
   printLevel <- slot(control, "printLevel")
   batchSize <- slot(control, "SG_batchSize")
   if(optimizer == "Adam") {
      maximType <- "Stochastic Gradient Ascent/Adam"
      Adam.momentum1 <- slot(control, "Adam_momentum1")
      Adam.momentum2 <- slot(control, "Adam_momentum2")
      Adam.delta <- 1e-8  # maybe make it a parameter in the future
      Adam.s <- 0
      Adam.r <- 0
      Adam.time <- 0
   } else if(optimizer == "SGA") {
      momentum <- slot(control, "SGA_momentum")
      v <- 0  # velocity that retains the momentum
   } else {
      stop(paste("unknown optimizer", optimizer))
   }
   ## ---------- How many batches
   if(is.null(batchSize)) {
      nBatches <- 1
      index <- seq(from=1, to=nObs, by=nBatches)
   } else {
      nBatches <- max(1L, nObs %/% batchSize)
                           # ensure that we get at least one batch if batchSize set too large
      shuffledIndex <- sample(nObs, nObs)
      index <- shuffledIndex[seq(from=1, to=nObs, by=nBatches)]
   }
   ##
   f1 <- NULL
                           # mark that we haven't computed the fcn value
   if(printLevel > 0) {
      f1 <- fn(start, fixed = fixed, sumObs = TRUE, index=index, ...)
      cat("Initial function value:", f1, "\n")
      if( isTRUE( attr( f1, "gradBoth" ) ) ) {
         warning( "the gradient is provided both as attribute 'gradient' and",
                 " as argument 'grad': ignoring argument 'grad'" )
      }
      if( isTRUE( attr( f1, "hessBoth" ) ) ) {
         warning( "the Hessian is provided both as attribute 'hessian' and",
                 " as argument 'hess': ignoring argument 'hess'" )
      }
   }
   if(!is.null(patience)) {
      if(is.null(f1)) {
         f1 <- fn(start, fixed = fixed, sumObs = TRUE, index=index, ...)
      }
      fBest <- f1  # remember the previous best value
      paramBest <- start
      patienceCount <- 0
                           # how many times have we hit a worse outcome
   }
   G1 <- grad(start, fixed = fixed, sumObs = TRUE, index=index, ...)
                           # have to compute fn as we cannot get gradient otherwise
   if(any(is.na(G1[!fixed]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(G1[!fixed]))) {
      stop("Infinite initial gradient")
   }
   if(length(G1) != nParam) {
      stop( "length of gradient (", length(G1),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   if(length(clip) > 0) {
      if((norm2 <- sum(G1*G1)) > clip)
         G1 <- G1/sqrt(norm2)*sqrt(clip)
   }
   if(storeValues) {
      valueStore <- rep(NA_real_, iterlim + 1)
      if(is.null(f1)) {
         f1 <- fn(start, fixed = fixed, sumObs = TRUE, index=index, ...)
      }
      valueStore[1] <- f1
   }
   if(storeParameters) {
      parameterStore <- matrix(NA_real_, iterlim + 1, nParam)
      dimnames(parameterStore) <- list(epoch=c("start", 1:iterlim), parameter=names(start))
      parameterStore[1,] <- start
   }
   if(printLevel > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(f1), "\n")
      a <- cbind(start1, G1, as.integer(!fixed))
      dimnames(a) <- list(names(start1), c("parameter", "initial gradient",
                                           "free"))
      printRowColLimits(a, max.rows, max.cols)
   }
   ## ---------------- Main interation loop ------------------------
   iter <- 0L
   ## we do not need to compute the function itself here, except for
   ## printing
   repeat {
                           # repeat over epochs
      ## break here if iterlim == 0
      if( iter >= iterlim) {
         code <- 4; break
      }
      ## break here to avoid potentially costly gradient computation
      if( iter >= slot(control, "iterlim")) {
         code <- 4; break
      }
      iter <- iter + 1L
      if(printLevel > 1) {
         cat( "----- epoch", iter, "-----\n")
      }
      for(iBatch in 1:nBatches) {
                           # repeat over minibatches
         if(!is.null(batchSize)) {
            index <- shuffledIndex[seq(from=iBatch, to=nObs, by=nBatches)]
         }
         start0 <- start1
         G0 <- G1
         if(any(is.na(G0[!fixed]))) {
            stop("NA in gradient")
         }
         if(optimizer == "SGA") {
            v <- momentum*v + learningRate*G0
            start1 <- start0 + v
         } else if(optimizer == "Adam") {
            Adam.time <- Adam.time + 1
            Adam.s <- Adam.momentum1*Adam.s + (1 - Adam.momentum1)*G0
            Adam.r <- Adam.momentum2*Adam.r + (1 - Adam.momentum2)*G0*G0
            Adam.shat <- Adam.s/(1 - Adam.momentum1^Adam.time)
            Adam.rhat <- Adam.r/(1 - Adam.momentum2^Adam.time)
            v <- learningRate*Adam.shat/(sqrt(Adam.rhat) + Adam.delta)
            start1 <- start0 + v
         }
         f1 <- NULL
                           # we are at a new location, mark that we haven't computed the f1 values
         ## still iterations to go, hence compute gradient
         G1 <- grad(start1, fixed = fixed, sumObs = TRUE, index=index, ...)
         if(any(is.na(G1[!fixed])) || any(is.infinite(G1[!fixed]))) {
            cat("Iteration", iter, "\n")
            cat("Parameter:\n")
            print(headDots(start1, max.cols), quote=FALSE)
            cat("Gradient:\n")
            printRowColLimits(G1, max.rows, max.cols)
            stop("NA/Inf in gradient")
         }
         if(length(clip) > 0) {
            if((norm2 <- sum(G1*G1)) > clip)
                           # compute norm w/o cross-product as grad may not be a vector
               G1 <- G1/sqrt(norm2*clip)
         }
         ## print every batch if someone wants...
         if(printLevel > 4) {
            f1 <- fn(start1, fixed = fixed, sumObs = TRUE, index=index, ...)
            cat(" - batch", iBatch, "index", index, "learning rate", learningRate, " fcn value:",
                formatC(as.vector(f1), digits=8, format="f"),  "\n")
            a <- cbind(learningRate*G0, start1, G1, as.integer(!fixed))
            dimnames(a) <- list(names(start0), c("delta-v", "param",
                                                 "gradient", "active"))
            printRowColLimits(a, max.rows, max.cols)
         }
         if(any(is.infinite(G1))) {
            code <- 6; break;
         }
      }  # end of repeat over batches
      if(storeValues) {
         ## store last value of the epoch
         if(is.null(f1)) {
            f1 <- fn(start1, fixed = fixed, sumObs = TRUE, index=index, ...)
         }
         valueStore[iter + 1L] <- c(f1)
                           # c removes dimensions and attributes
      }
      if(storeParameters) {
         ## store last value of the epoch
         parameterStore[iter + 1L,] <- c(start1)
                           # c removes dimensions and attributes
      }
      if(slot(control, "printLevel") > 2) {
         if(is.null(f1)) {
            f1 <- fn(start1, fixed = fixed, sumObs = TRUE, index=index, ...)
         }
         cat(" learning rate", learningRate, " fcn value:",
            formatC(as.vector(f1), digits=8, format="f"),  "\n")
         a <- cbind(learningRate*G0, start1, G1, as.integer(!fixed))
         dimnames(a) <- list(names(start0), c("amount", "param",
                                             "gradient", "active"))
         printRowColLimits(a, max.rows, max.cols)
      }
      ## stopping criteria
      if( sqrt( crossprod( G1[!fixed] ) ) < slot(control, "gradtol") ) {
         code <-1; break
      }
      if(!is.null(patience) && (iter %% patienceStep == 0)) {
         if(is.null(f1)) {
            f1 <- fn(start1, fixed = fixed, sumObs = TRUE, index=index, ...)
         }
         if(f1 < fBest) {
            patienceCount <- patienceCount + 1
         } else {
            patienceCount <- 0
            fBest <- f1
            paramBest <- start1
         }
         if(patienceCount > patience) {
            code <- 10
            f1 <- fBest
            start1 <- paramBest
            break
         }
      }
   }  # main iteration loop over epochs
   if(printLevel > 0) {
      cat( "--------------\n")
      cat( maximMessage( code), "\n")
      cat( iter, " iterations\n")
      cat( "estimate:", headDots(start1, max.cols), "\n")
      if(is.null(f1)) {
         f1 <- fn(start1, fixed = fixed, sumObs = TRUE, index=index, ...)
      }
      cat( "Function value:", f1, "\n")
   }
   if(finalHessian & !bhhhHessian) {
      G1 <- grad( start1, fixed = fixed, sumObs = FALSE, index=index, ... )
   }
   if(observationGradient(G1, length(start1))) {
      gradientObs <- G1
      colnames( gradientObs ) <- names(start1)
      G1 <- colSums(as.matrix(G1 ))
   }
   else {
      gradientObs <- NULL
   }
   names( G1 ) <- names(start1)
   ## calculate (final) Hessian
   if(tolower(finalHessian) == "bhhh") {
      if(!is.null(gradientObs)) {
         hessian <- - crossprod( gradientObs )
         attr(hessian, "type") <- "BHHH"
      } else {
         hessian <- NULL
         warning("For computing the final Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   } else if( finalHessian != FALSE ) {
      hessian <- hess( start1, fixed = fixed, index=index, ... )
   } else {
       hessian <- NULL
   }
   if( !is.null( hessian ) ) {
      rownames( hessian ) <- colnames( hessian ) <- names(start1)
   }

   ## remove attributes from final value of objective (likelihood) function
   attributes( f1 )$gradient <- NULL
   attributes( f1 )$hessian <- NULL
   attributes( f1 )$gradBoth <- NULL
   attributes( f1 )$hessBoth <- NULL
   ##
   result <- list(
      maximum = unname( drop( f1 ) ),
      estimate=start1,
      gradient=drop(G1),
      hessian=hessian,
      code=code,
      message=maximMessage( code),
      fixed=fixed,
      iterations=iter,
      type=maximType,
      valueStore = if(storeValues) valueStore else NULL,
      parameterStore = if(storeParameters) parameterStore else NULL
   )
   if( exists( "gradientObs" ) ) {
      result$gradientObs <- gradientObs
   }
   result <- c(result, control=control)
                           # attach the control parameters
   ##
   class(result) <- c("maxim", class(result))
   invisible(result)
}
