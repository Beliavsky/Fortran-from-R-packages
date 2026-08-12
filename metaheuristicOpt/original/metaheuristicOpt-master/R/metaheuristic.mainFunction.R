#############################################################################
#
#  This file is a part of the R package "metaheuristicOpt".
#
#  Author: Iip
#  Co-author: Muhammad Bima Adi Prabowo
#  Supervisors: Lala Septem Riza, Eddy Prasetyo Nugroho, Enjun Junaeti
#
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' A main funtion to compute the optimal solution using a selected algorithm.
#'
#' This function makes accessible all algorithm that are implemented
#' in this package. All of the algorithm use this function as interface to find
#' the optimal solution, so users do not need to call other functions.
#' In order to obtain good results, users need to adjust some parameters such as the
#' objective function, optimum type, number variable or dimension, number populations,
#' the maximal number of iterations, lower bound, upper bound, or other algorithm-dependent parameters
#' which are collected in the control parameter.
#'
#' @title metaOpt The main function to execute algorithms for getting optimal solutions
#'
#' @param FUN an objective function or cost function,
#'
#' @param optimType a string value that represents the type of optimization.
#'        There are two options for this arguments: \code{"MIN"} and \code{"MAX"}.
#'        The default value is \code{"MIN"}, referring the minimization problem.
#'        Otherwise, you can use \code{"MAX"} for maximization problem.
#'
#' @param algorithm a vector or single string value that represent the algorithm used to
#'        do optimization. There are currently twenty one implemented algorithm:
#' \itemize{
#' \item \code{"PSO"}: Particle Swarm Optimization. See \code{\link{PSO}};
#' \item \code{"ALO"}: Ant Lion Optimizer. See \code{\link{ALO}};
#' \item \code{"GWO"}: Grey Wolf Optimizer. See \code{\link{GWO}}
#' \item \code{"DA"} : Dragonfly Algorithm. See \code{\link{DA}}
#' \item \code{"FFA"}: Firefly Algorithm. See \code{\link{FFA}}
#' \item \code{"GA"} : Genetic Algorithm. See \code{\link{GA}}
#' \item \code{"GOA"}: Grasshopper Optimisation Algorithm. See \code{\link{GOA}}
#' \item \code{"HS"}: Harmony Search Algorithm. See \code{\link{HS}}
#' \item \code{"MFO"}: Moth Flame Optimizer. See \code{\link{MFO}}
#' \item \code{"SCA"}: Sine Cosine Algorithm. See \code{\link{SCA}}
#' \item \code{"WOA"}: Whale Optimization Algorithm. See \code{\link{WOA}}
#' \item \code{"CLONALG"}: Clonal Selection Algorithm. See \code{\link{CLONALG}}
#' \item \code{"DE"}: Differential Evolution Algorithm. See \code{\link{DE}}
#' \item \code{"SFL"}: Shuffled Frog Leaping Algorithm. See \code{\link{SFL}}
#' \item \code{"CSO"}: Cat Swarm Optimization Algorithm. See \code{\link{CSO}}
#' \item \code{"ABC"}: Artificial Bee Colony Algorithm. See \code{\link{ABC}}
#' \item \code{"KH"}: Krill-Herd Algorithm. See \code{\link{KH}}
#' \item \code{"CS"}: Cuckoo Search Algorithm. See \code{\link{CS}}
#' \item \code{"BA"}: Bat Algorithm. See \code{\link{BA}}
#' \item \code{"GBS"}: Gravitation Based Search Algorithm. See \code{\link{GBS}}
#' \item \code{"BHO"}: Black Hole Based Optimization Algorithm. See \code{\link{BHO}}
#' }
#'
#' @param numVar a positive integer to determine the number variables.
#'
#' @param rangeVar a matrix (\eqn{2 \times n}) containing the range of variables,
#'        where \eqn{n} is the number of variables, and first and second rows
#'        are the lower bound (minimum) and upper bound (maximum) values, respectively.
#'        If all variable have equal upper bound, you can define \code{rangeVar} as
#'        matrix (\eqn{2 \times 1}).
#'
#' @param control a list containing all arguments, depending on the algorithm to use. The following list are
#'        parameters required for each algorithm.
#'
#' \itemize{
#'
#' \item \code{PSO}:
#'
#'     \code{list(numPopulation, maxIter, Vmax, ci, cg, w)}
#'
#' \item \code{ALO}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#' \item \code{GWO}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#' \item \code{DA}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#' \item \code{FFA}:
#'
#'     \code{list(numPopulation, maxIter, B0, gamma, alphaFFA)}
#'
#' \item \code{GA}:
#'
#'     \code{list(numPopulation, maxIter, Pm, Pc)}
#'
#' \item \code{GOA}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#' \item \code{HS}:
#'
#'     \code{list(numPopulation, maxIter, PAR, HMCR, bandwith)}
#'
#' \item \code{MFO}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#' \item \code{SCA}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#' \item \code{WOA}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#' \item \code{CLONALG}:
#'
#'     \code{list(numPopulation, maxIter, selectionSize, multipicationFactor, hypermutationRate)}
#'
#' \item \code{DE}:
#'
#'     \code{list(numPopulation, maxIter, scalingVector, crossOverRate, strategy)}
#'
#' \item \code{SFL}:
#'
#'     \code{list(numPopulation, maxIter, numMemeplex, frogLeapingIteration)}
#'
#' \item \code{CSO}:
#'
#'     \code{list(numPopulation, maxIter, mixtureRatio, tracingConstant, maximumVelocity, smp, srd, cdc, spc)}
#'
#' \item \code{ABC}:
#'
#'     \code{list(numPopulation, maxIter, cycleLimit)}
#'
#' \item \code{KH}:
#'
#'     \code{list(numPopulation, maxIter, maxMotionInduced, inertiaWeightOfMotionInduced, epsilon, foragingSpeed, inertiaWeightOfForagingSpeed, maxDifussionSpeed, constantSpace, mu)}
#'
#'  \item \code{CS}:
#'
#'     \code{list(numPopulation, maxIter, abandonedFraction)}
#'
#'  \item \code{BA}:
#'
#'     \code{list(numPopulation, maxIter, maxFrequency, minFrequency, gama, alphaBA)}
#'
#'  \item \code{GBS}:
#'
#'     \code{list(numPopulation, maxIter, gravitationalConst, kbest)}
#'
#'  \item \code{BHO}:
#'
#'     \code{list(numPopulation, maxIter)}
#'
#'
#' \bold{Description of the \code{control} Parameters}
#' \itemize{
#' \item \code{numPopulation}: a positive integer to determine the number populations.
#'       The default value is 40.
#'
#' \item \code{maxIter}: a positive integer to determine the maximum number of iterations.
#'       The default value is 500.
#'
#' \item \code{Vmax}: a positive integer to determine the maximum velocity of particle.
#'       The default value is 2.
#'
#' \item \code{ci}: a positive integer to determine the individual cognitive.
#'       The default value is 1.49445.
#'
#' \item \code{cg}: a positive integer to determine the group cognitive.
#'       The default value is 1.49445.
#'
#' \item \code{w}: a positive integer to determine the inertia weight.
#'       The default value is 0.729.
#'
#' \item \code{B0}: a positive integer to determine the attractiveness firefly at r=0.
#'       The default value is 1.
#'
#' \item \code{gamma}: a positive integer to determine light absorption coefficient.
#'       The default value is 1.
#'
#' \item \code{alphaFFA}: a positive integer to determine randomization parameter.
#'       The default value is 0.2.
#'
#' \item \code{Pm}: a positive integer to determine mutation probability.
#'       The default value is 0.1.
#'
#' \item \code{Pc}: a positive integer to determine crossover probability.
#'       The default value is 0.8.
#'
#' \item \code{PAR}: a positive integer to determine Pinch Adjusting Rate.
#'       The default value is 0.3.
#'
#' \item \code{HMCR}: a positive integer to determine Harmony Memory Considering Rate.
#'       The default value is 0.95.
#'
#' \item \code{bandwith}: a positive integer to determine distance bandwith.
#'       The default value is 0.05.
#'
#' \item \code{selectionSize}: a positive integer between 0 and numVar to determine selection size.
#'       The default value is \code{as.integer(numPopulation/4)}.
#'
#' \item \code{multipicationFactor}: a positive numeric between 0 and 1 to determine number of clones.
#'       The default value is 0.5.
#'
#' \item \code{hypermutationRate}: a positive numeric between 0 and 1 to determine probabilty of variable in
#'        clone candidate solutions to be mutated, close to 1 probability is high and vice versa.
#'        The default value is 0.1.
#'
#' \item \code{scalingVector}: a positive numeric between 0 and 1 to determine scalingVector for mutation operator.
#'       The default value is 0.8.
#'
#' \item \code{crossOverRate}: a positive numeric between 0 and 1 to determine crossOver probability.
#'       The default value is 0.5.
#'
#' \item \code{strategy}: characters to determine mutation method. They are six methods to choose:
#'    \itemize{
#'    \item "classical".
#'    \item "best 1"
#'    \item "target to best"
#'    \item "best 2"
#'    \item "rand 2"
#'    \item "rand 2 dir"
#'    }
#'    The default value is "best 1".
#'
#' \item \code{numMemeplex}: a positive integer (as.integer()) between 0 and numVar to
#'        determine number of memeplex (see details).The default value is \code{as.integer(numPopulation/3)}.
#'
#' \item \code{frogLeapingIteration}: a positive integer (as.integer()) to determine number
#'        of iteration for each memeplex. The default value is \code{as.integer(10)}.
#'
#' \item \code{mixtureRatio}: a positive numeric between 0 and 1 to determine flaging proportion.
#'        higher mixtureRatio increase number of candidate solutions in seeking mode
#'        and vice versa. The default value is 0.5.
#'
#' \item \code{tracingConstant}: a positive numeric between 0 and 1 to determine tracingConstant.
#'        The default value is 0.1.
#'
#' \item \code{maximumVelocity}: a positive numeric to determine maximumVelocity while candidate solutions
#'        in tracing mode performing local search. The default value is 1.
#'
#' \item \code{smp}: a positive integer to determine number of duplication in genetic operator.
#'        The default value is \code{as.integer(20)}.
#'
#' \item \code{srd}: a positive numeric between 0 and 100 to determine mutation length in genetic operator.
#'        The default value is 20.
#'
#' \item \code{cdc}: a positive integer between 0 and numVar to determine number of variabel in
#'        candidate solutions in seeking mode to be mutated during mutation step in
#'        genetic operator. The default value is \code{as.integer(numVar)}.
#'
#' \item \code{spc}: a logical. if spc is TRUE smp = smp else smp = smp - 1. The default value is TRUE.
#'
#' \item \code{cycleLimit}: a positive integer to determine number of times allowed for
#'        candidate solution to not move. The default value is \code{as.integer(numVar * numPopulation)}.
#'
#' \item \code{maxMotionInduced}: a positive numeric between 0 and 1 to determine
#'        maximum motion induced. The default value is 0.01.
#'
#' \item \code{inertiaWeightOfMotionInduced}: a positive numeric between 0 and 1 to determine
#'        how much motion induced affect krill (candidate solution) movement. the
#'        greater the value the greater the affect of motion induced on krill movement.
#'        The default value is 0.01.
#'
#' \item \code{epsilon}: a positive numeric between 0 and 1 to determine epsilon constant.
#'        The default value is 1e-05.
#'
#' \item \code{foragingSpeed}: a positive numeric between 0 and 1 to determine foraging speed.
#'        The default value is 0.02
#'
#' \item \code{inertiaWeightOfForagingSpeed}: a positive numeric between 0 and 1 to determine
#'        how much foraging speed affect krill (candidate solution) movement. the
#'        greater the value the greater the affect of foraging speed on krill movement.
#'        The default value is 0.01.
#'
#' \item \code{maxDifussionSpeed}: a positive numeric between 0 and 1 to determine maximum
#'        difussion speed. The default value is 0.01.
#'
#' \item \code{constantSpace}: a numeric between 0 and 1 to determine how much range affect
#'        krill movement. The default value is 1.
#'
#' \item \code{mu}: a numeric between 0 and 1 to determine constant number for mutation operator.
#'       The default value is 0.1.
#'
#' \item \code{abandonedFraction}: a positive numeric between 0 and 1 to determine fraction
#'        of population to be replaced. The default value is 0.5.
#'
#' \item \code{maxFrequency}: a numeric to determine maximum frequency. The default value is 0.1.
#'
#' \item \code{minFrequency}: a numeric to determine minimum frequency. The default value is -0.1.
#'
#' \item \code{gama}: a numeric greater than equal to 1. It use to increase pulse rate. The default value is 1.
#'
#' \item \code{alphaBA}: a numeric between 0 and 1. It use to decrease loudness. The default value is 0.1.
#'
#' \item \code{gravitationalConst}: a numeric to determine gravitational constant while
#'        calculating total force. The default value is \code{max(rangeVar)}.
#'
#' \item \code{kbest}: a positive numeric between 0 and 1 to determine fraction of population
#'        with best fitness which will affect every candidate solution in population.
#'        The default value is 0.1.
#' }
#' }
#'
#' @param seed a number to determine the seed for RNG.
#'
#' @examples
#' ##################################
#' ## Optimizing the sphere function
#'
#' ## Define sphere function as an objective function
#' sphere <- function(X){
#'     return(sum(X^2))
#' }
#'
#' ## Define control variable
#' control <- list(numPopulation=40, maxIter=100, Vmax=2, ci=1.49445, cg=1.49445, w=0.729)
#'
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#'
#' ## Define control variable
#' best.variable <- metaOpt(sphere, optimType="MIN", algorithm="PSO", numVar,
#'                          rangeVar, control)
#'
#' @return \code{List} that contain list of variable, optimum value and execution time.
#'
#' @export

metaOpt <- function(FUN, optimType="MIN", algorithm="PSO", numVar, rangeVar, control=list(), seed=NULL){

	## get optimType
	optimType <- toupper(optimType)

	## get algorithm
	algorithm <- toupper(algorithm)

	## initialize result
	result <- matrix(ncol=numVar, nrow=length(algorithm))

	## initialize time elapsed
	timeElapsed <- matrix(ncol=3, nrow=length(algorithm))

	## checking consistency between variable numVar and rangeVar
	if(numVar != ncol(rangeVar) & ncol(rangeVar) != 1){
		stop("Inconsistent between number variable and number range variable")
	}


	for(i in 1:length(algorithm)){

		## PSO Algorithm
		if(algorithm[i] == "PSO"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500,
						Vmax=2, ci=1.49445, cg=1.49445, w=0.729))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter
			Vmax <- control$Vmax
			ci <- control$ci
			cg <- control$cg
			w <- control$w

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- PSO(FUN, optimType, numVar, numPopulation, maxIter, rangeVar, Vmax, ci, cg, w)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Ant Lion Optimizer Algorithm
		else if(algorithm[i] == "ALO"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- ALO(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Grey Wolf Optimizer Algorithm
		else if(algorithm[i] == "GWO"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- GWO(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Dragonfly Algorithm
		else if(algorithm[i] == "DA"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- DA(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Firefly Algorithm
		else if(algorithm[i] == "FFA"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500, B0=1, gamma=1, alpha=0.2))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter
			B0 <- control$B0
			gamma <- control$gamma
			alpha <- control$alpha

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- FFA(FUN, optimType, numVar, numPopulation, maxIter, rangeVar, B0, gamma, alpha)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Genetic Algorithm
		else if(algorithm[i] == "GA"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500, Pm=0.1, Pc=0.8))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter
			Pm <- control$Pm
			Pc <- control$Pc

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- GA(FUN, optimType, numVar, numPopulation, maxIter, rangeVar, Pm, Pc)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Grasshopper Optimisation Algorithm
		else if(algorithm[i] == "GOA"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- GOA(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Harmony Search Algorithm
		else if(algorithm[i] == "HS"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500, PAR=0.3, HMCR=0.95, bandwith=0.05))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- HS(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Moth Flame Optimizer
		else if(algorithm[i] == "MFO"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- MFO(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Sine Cosine Algorithm
		else if(algorithm[i] == "SCA"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter

			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- SCA(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Whale Optimization Algorithm
		else if(algorithm[i] == "WOA"){
			## checking missing parameters
			control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

			## get all parameter
			numPopulation <- control$numPopulation
			maxIter <- control$maxIter
			# generate result while calculating time elapsed
			set.seed(seed)
			temp<-system.time(
				result[i,] <- WOA(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
			)
			temp <- c(temp[1], temp[2], temp[3])
			timeElapsed[i,]=temp;
		}
		# Clonal Selection Algorithm
		else if(algorithm[i] == "CLONALG"){
		  ## checking missing parameters
		  control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

		  ## get all parameter
		  numPopulation <- control$numPopulation
		  maxIter <- control$maxIter
		  # generate result while calculating time elapsed
		  set.seed(seed)
		  temp<-system.time(
		    result[i,] <- CLONALG(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
		  )
		  temp <- c(temp[1], temp[2], temp[3])
		  timeElapsed[i,]=temp;
		}
	  # Artificial Bee Colony Algorithm
	  else if(algorithm[i] == "ABC"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- ABC(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Bat Algorithm
	  else if(algorithm[i] == "BA"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- BA(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Cuckoo Search
	  else if(algorithm[i] == "CS"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- CS(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Cat Swarm Optimization
	  else if(algorithm[i] == "CSO"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- CSO(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Differential Evolution
	  else if(algorithm[i] == "DE"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- DE(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Gravitational Based Search Algorithm
	  else if(algorithm[i] == "GBS"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- GBS(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Krill-Heard Algorithm
	  else if(algorithm[i] == "KH"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- KH(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Shuffled Frog Leaping Algorithm
	  else if(algorithm[i] == "SFL"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- SFL(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;
	  }
	  # Black Hole-based Algorithm
	  else if(algorithm[i] == "BHO"){
	    ## checking missing parameters
	    control <- setDefaultParametersIfMissing(control, list(numPopulation=40, maxIter=500))

	    ## get all parameter
	    numPopulation <- control$numPopulation
	    maxIter <- control$maxIter
	    # generate result while calculating time elapsed
	    set.seed(seed)
	    temp<-system.time(
	      result[i,] <- BHO(FUN, optimType, numVar, numPopulation, maxIter, rangeVar)
	    )
	    temp <- c(temp[1], temp[2], temp[3])
	    timeElapsed[i,]=temp;

		}else{
			stop("unknown Algorithm argument value")
		}
	}

	# generating optimum value foreach algorithm
	optimumValue <- c()
	for (i in 1:nrow(result)) {
		optimumValue[i] <- FUN(result[i,])
	}
	optimumValue <- as.matrix(optimumValue)

	# set name for each row
	rownames(result) <- algorithm
	rownames(optimumValue) <- algorithm
	rownames(timeElapsed) <- algorithm

	#set name for column
	colName <- c()
	for (i in 1:numVar) {
		colName[i] <- paste("var",i,sep="")
	}
	colnames(result) <- colName
	colnames(optimumValue) <- c("optimum_value")
	colnames(timeElapsed) <- c("user", "system", "elapsed")

	# build list
	allResult <- list(result=result, optimumValue=optimumValue, timeElapsed=timeElapsed)

	return(allResult)
}

## checking missing parameters
# @param control parameter values of each algorithm
# @param defaults default parameter values of each algorithm

setDefaultParametersIfMissing <- function(control, defaults) {
  for(i in names(defaults)) {
    if(is.null(control[[i]])) control[[i]] <- defaults[[i]]
  }
  control
}
