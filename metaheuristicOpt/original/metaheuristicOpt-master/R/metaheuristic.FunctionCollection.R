#############################################################################
#
#  This file is a part of the R package "metaheuristicOpt".
#
#  Author: Iip
#  Co-author: -
#  Supervisors: Lala Septem Riza, Eddy Prasetyo Nugroho
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

# this function will generate the random number with defined boundary
# @param numPopulation number population / number row
# @param dimension number variable / number column
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable

generateRandom <- function(numPopulation, dimension, lowerBound, upperBound){
	result <- matrix()
	if(length(lowerBound)==1){
		result <- matrix(runif(numPopulation*dimension, lowerBound, upperBound), nrow=numPopulation, ncol=dimension)

	}else{
		result <- matrix(nrow=numPopulation, ncol=dimension)
		for (i in 1:dimension){
			result[,i] = runif(numPopulation, lowerBound[i], upperBound[i])
		}
	}
	return(result)
}

# this function is for calculating the fitness
# @param FUN objective function
# @param optimType type optimization
# @param popu population of candidate

calcFitness <- function(FUN, optimType, popu){
	fitness <- c()
	for (i in 1:nrow(popu)) {
		fitness[i] <- optimType*FUN(popu[i,])
	}
	return(fitness)
}

# this function is for calculating the best fitness
# @param FUN objective function
# @param optimType type optimization
# @param popu population of candidate

calcBest <- function(FUN, optimType, popu){
	fitness <- calcFitness(FUN, optimType, popu)
	best <- popu[which.max(fitness),]
	return(best)
}

# this function used to check the boundary for each dimension/variable
# @param position is vector of position
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable

checkBound <- function(position, lowerBound, upperBound){
	check1 <- position > upperBound
	check2 <- position < lowerBound

	result <- (position*(!(check1+check2)))+upperBound*check1+lowerBound*check2
	return(result)
}

# function to get index resulting roulette whell selection method
# @param weight vector of double

rouletteWhell <- function(weight){
	# handle negative number
	if(any(weight<0)){
		c <- -weight[which.min(weight)]
		weight <- weight+c+1
	}
	# count the cumulative sum
	accumulation <- cumsum(weight)
	# pick a random number
	r <- runif(1)*accumulation[length(accumulation)]

	# index to pick [default value is 1]
	result <- 1
	for (i in 1:length(accumulation)){
		if(accumulation[i] > r){
			result <- i
			break
		}
	}
	return(result)
}

# example of test function sphere
sphere <- function(X){
	X <- sum(X^2)
	return(X)
}

# example of test function schwefel
schwefel <- function(X){
	X <- sum(-X*sin(sqrt(abs(X))))
	return(X)
}

# example of test function rastrigin
rastrigin <- function(X){
	X <- sum((X^2)-10*cos(2*pi*X)+10)
	return(X)
}

# example of test function for determing the centilever beam
# the real problem
centileverBeam <- function(X){
	FX <- 0.6224*sum(X)
	x3 <- X^3
	up <- c(61,37,19,7,1)
	part <- sum(up/x3) - 1
	decision <- max(part,0)
	res <- FX + (10e+17*decision)
	return(res)
}

## p-1 p-2 p-3

MAE <- function(dataHarga, X){
	value <- dataHarga[,4] ## get value
	dataHarga <- dataHarga[,1:3] ## get feature
	dataHarga <- cbind(dataHarga, rep(1, nrow(dataHarga)))
	matrix_of_x <- matrix(rep(X,nrow(dataHarga)), nrow = nrow(dataHarga), byrow = TRUE)
	result <- rowSums(dataHarga*matrix_of_x)
	# calculating MAE
	result <- mean(abs(result-value))
	return(result)
}