#' Evaluates a Two-Dimensional Monte Carlo Model in a Loop
#'
#' \samp{evalmccut} evaluates a Two-Dimensional Monte Carlo model using a loop on the
#' uncertainty dimension. Within each loop, it calculates statistics in the variability
#' dimension and stores them for further analysis. It allows evaluation of very high
#' dimension models using time instead of memory.
#'
#' \samp{x} (or \samp{model}) should be built as three blocks, separated by \samp{\{}:
#' \enumerate{
#'   \item The first block is evaluated once (and only once) before the loop. It should
#'     evaluate all \samp{"V"}, \samp{"0"} and \samp{"U"} nodes.
#'   \item The second block, which should lead to an \samp{mc} object, is evaluated using
#'     \samp{nsu=1}. It must end with an expression like \samp{mymc <- mc(...)}.
#'   \item The third block is evaluated on the \samp{mc} object. All resulting statistics
#'     are stored. It should build a list of statistics such as \samp{summary}, \samp{plot},
#'     \samp{tornado} or any function leading to a vector, matrix or list of those.
#' }
#'
#' Steps 2 and 3 are repeated \samp{nsu} times. At each iteration, the loop index
#' (from 1 to \samp{nsu}) is assigned to the variable specified in \samp{ind}.
#'
#' @note The seed is set at the beginning of the evaluation. Complete similarity of two
#' evaluations with the same seed is not certain depending on the model structure.
#' In particular, the results will differ from \code{\link{evalmcmod}} with the same seed.
#' Do not use upper case variables in your model to avoid conflicts with internal variables.
#'
#' @param model a \samp{mcmodelcut} object or a valid call including three blocks.
#' @param nsv the number of simulations for variability.
#' @param nsu the number of simulations for uncertainty.
#' @param seed the random seed. If \samp{NULL} the seed is unchanged.
#' @param ind the variable name used in \samp{model} to refer to the uncertainty loop index.
#' @return An object of class \samp{mccut}.
#' @seealso \code{\link{evalmcmod}}, \code{\link{mcmodelcut}}.
#' @keywords methods
#' @name mccut
#' @examples
#' modEC3 <- mcmodelcut({
#' ## First block: evaluates all 0, V and U nodes
#'   {
#'     cook <- mcstoc(rempiricalD, type="V", values=c(0, 1/5, 1/50), prob=c(0.027, 0.373, 0.6))
#'     serving <- mcstoc(rgamma, type="V", shape=3.93, rate=0.0806)
#'     conc <- mcstoc(rnorm, type="U", mean=10, sd=2)
#'     r <- mcstoc(runif, type="U", min=5e-04, max=0.0015)
#'   }
#' ## Second block: evaluates VU nodes, leads to mc object
#'   {
#'     expo <- conc * cook * serving
#'     dose <- mcstoc(rpois, type="VU", lambda=expo)
#'     risk <- 1 - (1 - r)^dose
#'     res <- mc(conc, cook, serving, expo, dose, r, risk)
#'   }
#' ## Third block: statistics
#'   {
#'     list(
#'       sum = summary(res),
#'       plot = plot(res, draw=FALSE),
#'       minmax = lapply(res, range)
#'     )
#'   }
#' })
#'
#' x <- evalmccut(modEC3, nsv=101, nsu=101, seed=666)
#' summary(x)
#' @export
evalmccut <- function(model, nsv = ndvar(), nsu = ndunc(), seed = NULL,  ind = "index")
{

  # DEFINE SOME FUNCTIONS
  # DEFINEOUT build the empty list that the results will fill

  DEFINEOUT <- function(ARG){
    LAFUNC <- function(ARGJBD){
      if(inherits(ARGJBD,"tornado")) ARGJBD <- ARGJBD$value
      if(is.list(ARGJBD)) return(DEFINEOUT(ARGJBD))
      if(is.data.frame(ARGJBD)) ARGJBD <- as.matrix(ARGJBD)
      if(is.array(ARGJBD)){
         DIMS <- dim(ARGJBD)
         if(inherits(ARGJBD,"mcnode")){
            DIMS <- c(DIMS[1],nsu,DIMS[3])
            NOM <- vector(mode="list",length=3)}
         else {if(length(DIMS) > 2)                                              #Other arrays = matrix
                stop("Array > 2 dimensions are not supported in this function")
          DIMS <- c(DIMS[1],nsu,DIMS[2])
          NOM <- list(rownames(ARGJBD),NULL,colnames(ARGJBD))}
      }
      else {
        DIMS <- c(length(ARGJBD),nsu,1)
        NOM <- list(names(ARGJBD),NULL,NULL)
      }
      return(array(NA,dim=DIMS,dimnames=NOM))
      }
    OUT <- lapply(ARG,LAFUNC)
    names(OUT) <- names(ARG)
    class(OUT) <- class(ARG)
    return(OUT)
  }

  # FONCSPEC is a tricky "autocall" function to fill the results
  FONCSPEC <- function(STAT, PREC, LOOP){
    if(inherits(STAT,"tornado")) {STAT <- STAT$value}
    if(is.list(PREC)){
      PREC <- mapply(FONCSPEC,STAT,PREC,MoreArgs=list(LOOP=LOOP),SIMPLIFY=FALSE)
    }
    else PREC[,LOOP,] <- STAT
    return(PREC)
    }

OLDV <- ndvar()
OLDU <- ndunc()
RESULTAT <- try({                                                                # Debut du try
  ndvar(nsv)
  ndunc(nsu)
  if(!is.null(seed)) set.seed(seed)

  #Analyse expression
  NBEXPR <- length(model[[1]])
  if(NBEXPR != 4) stop("The expression should include three blocks")

  #Evaluate the first block and stock the U and VU mcnode
  LSBEG <- ls()
  eval(model[[1]][[2]])
  LSEND <- ls()
  NEW <- LSEND[!(LSEND %in% LSBEG)]                                             # New objects built
  NEW <- NEW[sapply(NEW,function(x) inherits(get(x),"mcnode"))]                 # New mcnode built
  TYPE <- sapply(NEW, function(x) attr(get(x),which="type"))
  OUTM <- sapply(NEW, function(x) attr(get(x),which="outm"))
    cat("'0' mcnode(s) built in the first block:",NEW[TYPE == "0"],"\n")
    cat("'V' mcnode(s) built in the first block:",NEW[TYPE == "V"],"\n")
    cat("'U' mcnode(s) built in the first block:",NEW[TYPE == "U"],"\n")
    cat("'VU' mcnode(s) built in the first block:",NEW[TYPE == "VU"],"\n")
    cat("The 'U' and 'VU' nodes will be sent column by column in the loop\n")

  QUEL <- (TYPE == "U" | TYPE == "VU")
  OUTM <- OUTM[QUEL]
  NN    <- sum(QUEL)
  TORIG <- TYPE[QUEL]                                                           # Type of Origin
  TORIG1D <- ifelse(TORIG %in% c("0","U"),"0","V")                              # Corresponding Type in 1D
  NORIG <- NEW[QUEL]                                                            # Name of the node to save
  NSAVED <- paste("SAVE",NORIG,sep="")                                          # New names

  for(i in 1:NN)  assign(NSAVED[i], get(NORIG[i]))                              # Save the nodes

  # Prepare the loop
  ndunc(1)                                                                      # nsu = 1
  for(i in 1:5) cat("---------|")
  cat("\n")
  quelcroix <- as.integer(seq(1,nsu,length=50))                                 # define the step for the progress bar
  cat(rep("*",sum(1 == quelcroix)),sep="")                                    # progress bar
  flush.console()

# First simulation
  assign(ind,1)

  # Get the first row of the stocked nodes
  NODE <- vector(mode="list",length=NN)
  for(i in 1:NN) {
    NODE[[i]] <- get(NSAVED[i])[,1,,drop=FALSE]                                 # Can not be put in an mapply function
    attr(NODE[[i]],which="type") <- TORIG[1]                                    # First simulation: keep the original "type"
    attr(NODE[[i]],which="outm") <- OUTM[i]
    class(NODE[[i]])  <- "mcnode"
    assign(NORIG[i],NODE[[i]]) }                                                 # And assign them to their original name


  # Evaluate the model
    eval(model[[1]][[3]])

    NOM <- as.character(model[[1]][[3]][[length(model[[1]][[3]])]][[2]])        # name of the mc second block
    MCOBJ <- get(NOM)
    if(!is.mc(MCOBJ)) stop("The second block result is not an mc object")

    TYPEORI <- sapply(MCOBJ,attr,which="type")                                  #Original type
    TYPENEW <- ifelse(TYPEORI %in% c("0","U"),"0","V")                          #New type

 #Evaluate the model with original type and keep the first results
    PREMS <- eval(model[[1]][[4]])

 #Evaluate the model with new types
    NODE <- mapply("attr<-",NODE,"type",TORIG1D,SIMPLIFY=FALSE)                 # Assign the new type
    for(i in 1:NN) assign(NORIG[i],NODE[[i]])                                   # And assign them to their original name

    MCOBJ[] <- mapply("attr<-",MCOBJ,"type",TYPENEW,SIMPLIFY=FALSE)
    assign(NOM,MCOBJ)

    PREMSPRIM <- eval(model[[1]][[4]])

 #Define the final structure and stock the first simu
    SORTIE <- DEFINEOUT(PREMSPRIM)
    SORTIE <- mapply(FONCSPEC,PREMSPRIM,SORTIE,MoreArgs=list(LOOP=1),SIMPLIFY=FALSE)

# Other Simulations
  for(LOOP in 2:nsu){

    assign(ind,LOOP)

    for(i in 1:NN) {
      NODE[[i]][] <- get(NSAVED[i])[,LOOP,,drop=FALSE]
      assign(NORIG[i],NODE[[i]])
    }

    eval(model[[1]][[3]])
    MCOBJ <- get(NOM)
    MCOBJ[] <- mapply("attr<-",MCOBJ,which="type",value=TYPENEW,SIMPLIFY=FALSE)
    assign(NOM,MCOBJ)

    STAT <- eval(model[[1]][[4]])
    SORTIE <- mapply(FONCSPEC,STAT,SORTIE,MoreArgs=list(LOOP=LOOP),SIMPLIFY=FALSE)

    cat(rep("*",sum(LOOP == quelcroix)),sep="")                                    # progress bar
    flush.console()
    }

  # Post Production : Affecte les classes + petites modifications pour qq fonctions connues

  for(JBD in 1:length(SORTIE)){
    if(inherits(PREMS[[JBD]],"summary.mc")) {
      LESTYPES <- sapply(PREMS[[JBD]],"attr","type")
      SORTIE[[JBD]] <- mapply("attr<-",SORTIE[[JBD]],"type",LESTYPES,SIMPLIFY=FALSE)
      class(SORTIE[[JBD]]) <- c("summary.mccut")
      }

    else if(inherits(PREMS[[JBD]],"mcnode")) {
      attr(SORTIE[[JBD]],"type") <- attr(PREMS[[JBD]],"type")
      attr(SORTIE[[JBD]],"outm") <- attr(PREMS[[JBD]],"outm")
      class(SORTIE[[JBD]]) <- c("mcnode")
      }

    else if(inherits(PREMS[[JBD]],"tornado")) {
      PREMS[[JBD]]$value <- SORTIE[[JBD]]
      SORTIE[[JBD]] <- PREMS[[JBD]]
      class(SORTIE[[JBD]]) <- c("tornado.mccut")
      }

    else if(inherits(PREMS[[JBD]],"plotmc")){
      TYPEPLOT <- sapply(PREMS[[JBD]],"attr",which="type")
      SORTIE[[JBD]] <- mapply("attr<-",SORTIE[[JBD]],"type",TYPEPLOT,SIMPLIFY=FALSE)
      class(SORTIE[[JBD]]) <- c("plot.mccut")}

  }

  class(SORTIE) <- "mccut"

SORTIE}, silent=TRUE)     # Fin du try
ndvar(OLDV)
ndunc(OLDU)
cat("\n")
if(inherits(RESULTAT,"try-error")) stop(RESULTAT,call. = FALSE)
return(invisible(RESULTAT))
}
