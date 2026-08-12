\name{bin.packing}
\alias{bin.packing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Enumeration of all existing solutions for one-dimensional bin-packing problem}
\description{The algorithm used for this function is a permutational modification of First Find (FF) algorithm described in Martello and Toth (1990). However, there are significant differences. First, the algorithm suggested by Martello and Toth does not set an objective to find all possible optimal solutions while algorithm suggested here finds them all. Apparently, these changes result in a significant increase of required computing time because we need to analyse and process all possible item  permutations. Noteworthy, the time of the optimized algorithm is still polynomial. Second, Martello  and Toth used an iterative embodiment of the algorithm while a recursive function is used here in order to reduce computing time by advantageous employment of "lazy evaluation" feature of R program and in order to optimize the code of the script.

According to its name, our algorithm is build around "generation of a bin". The objective is achieved by finding the next set of items that fit into existing bin (FF choice). The combination of added items is selected from the solutions provided by \code{nilde} function. This method allows to optimize computing time since we process more than one item at a time when we call the function. Further optimization is achieved by canceling any recursive calls when the addition of a new bin results in the number of bins exceeding the currently found local optimum (minimal number of bins achieved so far).

The algorithm segregates two types of calls from the parent function. Namely, scenarios when the current bin is complete and incomplete are treated separately. If the remaining unused capacity of the bin is zero, i.e. the bin is complete, then we check if creation of a new bin pushes the number of bins above the current optimum. If the number is still optimal, then we start a new bin and generate the list of item clusters that can fit and move on.

If the remaining unused capacity if not zero (bin is incomplete), then, first, we try to complete the existing bin by adding items that can fit and, then, if it is not possible, we close existing bin (even if its unused capacity is more than zero) and start processing item clusters that would not fit the closed bin anyway. By doing so, solutions with smaller number of bins will be generated earlier. Thus, we will be able to find globally optimal number of bins a.s.a.p.

Obviously, the algorithm stops recursive calls only in two cases: either we have distributed all the items or we exceeded the optimal number of bins by adding next new bin. In both cases, we proceed with processing the next possible combination of items, i.e. process the next leaf on our decision-making tree.

As for robustness checks, we have tested several versions of our recursive algorithm. For example, Martello and Toth demonstrate that on sample of any complexity First Find Decreasing (FFD) algorithm leads to a significant decrease in computing time as compared to FF. Therefore, we tested FFD modification of the algorithm. However, in our case, since we are looking for all optimal solutions, implementation of FFD algorithm has yielded no results. Furthermore, we experimented with the format and list of variables transferred recursively. Specifically, a version of the algorithm that transfers only logical vector of scenarios to be processed resulted in increase of computing time.

The function demonstrates the best computing time for all the sampled scenarios of item weights and bin capacity. However, there are some limitations to be addressed. For example, if the initial set includes multiple items with the same weight but different IDs, then the output of \code{GenVagonE} will need to be filtered from seemingly different solutions. Yet, the filtering is not computationally demanding and definitely polynomial in terms of time.

Note: majority of input variables are pre-computed in advance, separately, see example.

}
\usage{
 bin.packing(input.a, input.n, bin.globals)

}
\arguments{
   \item{input.a}{a vector of items weights.}
   \item{input.n}{capacity of a bin.}   
   \item{bin.globals}{an environment for global variables.}   
}

\value{
%  \item{TrainList}{a list where each node summarizes information about obtained solutions for trains of length \code{N}. So, the %first element of the list is trains with only one bin, second - with two bins, and so on. Empty node means that there are no %solution of this length are found.}
  
%TrainList$code - solutions of the bin-packing problem. Each number is a position of an item in the input string (Input string %specifies item weights. Items are numbered 1, 2, 3, etc.)  Items included into a bin are separated by commas. Bins are separated %by space character. Trains are enclosed by double quotes.

%TrainList$trInef - "inefficiency", i.e. unused space of each train (each train is a solution of a problem)

%TrainList$vagInef - "inefficiency", i.e. unused space of each bin respectively.

\item{min.bins}{minimum number of bins required.}
\item{solution}{solutions of the bin-packing problem. Each number is a position of an item in the input string (Input string specifies item weights. Items are numbered 1, 2, 3, etc.)  Items included into one bin are separated by commas. Bins are separated by space character. Different solutions are enclosed by double quotes.}
\item{bin.ineff}{bin "inefficiency", i.e. unused space of each bin respectively, for every solution.}
\item{total.ineff}{total "inefficiency", i.e. unused space of every solution (sum of bin "inefficiencies" per solution).}

}

\author{
Rashid Makarov 
}
\references{
  Martello, S. and Toth, P. (1990) Knapsack Problems: Algorithms and Computer Implementations, Wiley, Chichester, 1990.

   Voinov, V., Makarov, R., Voinov, Y. (2019) An exact polynomial in time solution of the one-dimensional bin-packing problem. In: Christos H Skiadas (ed.) Proceedings of the ASMDA 2019, published by ISAST (Int. Society for the Advancement of Science and Technology), December 2019, pp. 787-798.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{nilde-package}}, \code{\link{get.partitions}}, \code{\link{get.subsetsum}}, \code{\link{nlde}}    
}

\examples{
library(nilde)
input.a <- c(70, 60, 50, 40, 30, 20, 10) # weights of items
input.n <- 100 # capacity of a bin
bin.globals <- new.env() # a new environment for global variables
bin.globals$OptVag <- length(input.a) # initial min # of bins
bin.globals$TrainList <- vector("list",length(input.a)) # output with solutions
g <- bin.packing(input.a, input.n,bin.globals)
g$min.bins # minimum number of bins
g$solution # solutions
}

\keyword{optimize}% at least one, from doc/KEYWORDS


