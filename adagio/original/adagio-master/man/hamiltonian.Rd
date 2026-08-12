\name{hamiltonian}
\alias{hamiltonian}
\title{
  Finds a Hamiltonian path or cycle
}
\description{
  A Hamiltionian path or cycle (a.k.a. Hamiltonian circuit) is a path through
  a graph that visits each vertex exactly once, resp. a closed path through
  the graph.
}
\usage{
hamiltonian(edges, start = 1, cycle = TRUE)
}
\arguments{
  \item{edges}{an edge list describing an undirected graph.}
  \item{start}{vertex number to start the path or cycle.}
  \item{cycle}{Boolean, should a path or a full cycle be found.}
}
\details{
  \code{hamiltonian()} applies a backtracking algorithm that is relatively
  efficient for graphs of up to 30--40 vertices. The edge list is first
  transformed to a list where the \code{i}-th component contains the list
  of all vertices connected to vertex \code{i}.

  The edge list must be of the form \code{c(v1, v2, v3, v2, ...)} meaning
  that there are edges \code{v1 --> v2, v3 --> v4}, etc., connecting these
  vertices. Therefore, an edge list has an even number of entries.

  If the function returns \code{NULL}, there is no Hamiltonian path or
  cycle. The function does not check if the graph is connected or not. And
  if \code{cycle = TRUE} is used, then there also exists an edge from the
  last to the first entry in the resulting path.

  Ifa Hamiltonian cycle exists in the graph it will be found whatever the
  starting vertex was. For a Hamiltonian path this is different and a 
  successful search may very well depend on the start.
}
\value{
  Returns a vector containing vertex number of a valid path or cycle,
  or \code{NULL} if no path or cycle has been found (i.e., does not exist);
  If a cycle was requested, there exists an edge from the last to the first
  vertex in this list of edges.
}
\references{
  Papadimitriou, Ch. H., and K. Steiglitz (1998). Optimization Problems:
  Algorithms and Complexity. Prentice-Hall/Dover Publications.
}
\author{
  Hans W. Borchers
}
\note{
  See the \code{igraph} package for more information about handling graphs and
  defining them through edge lists or other constructs.
}
\seealso{
  Package \code{igraph}
}
\examples{
## Dodekaeder graph
D20_edges <- c(
     1,  2,  1,  5,  1,  6,  2,  3,  2,  8,  3,  4,  3, 10,  4,  5,  4, 12,
     5, 14,  6,  7,  6, 15,  7,  8,  7, 16,  8,  9,  9, 10,  9, 17, 10, 11,
    11, 12, 11, 18, 12, 13, 13, 14, 13, 19, 14, 15, 15, 20, 16, 17, 16, 20, 
    17, 18, 18, 19, 19, 20)
hamiltonian(D20_edges, cycle = TRUE)
# [1]  1  2  3  4  5 14 13 12 11 10  9  8  7 16 17 18 19 20 15  6
hamiltonian(D20_edges, cycle = FALSE)
# [1]  1  2  3  4  5 14 13 12 11 10  9  8  7  6 15 20 16 17 18 19

## Herschel graph
# The Herschel graph the smallest non-Hamiltonian polyhedral graph.
H11_edges <- c(
     1,  2,  1,  8,  1,  9,  1, 10,  2,  3,  2, 11,  3,  4,  3,  9,  4,  5,
     4, 11,  5,  6,  5,  9,  5, 10,  6,  7,  6, 11,  7,  8,  7, 10,  8, 11)
hamiltonian(H11_edges, cycle = FALSE)
# NULL

\dontrun{
## Example: Graph constructed from squares
N <- 45  # 23, 32, 45
Q <- (2:trunc(sqrt(2*N-1)))^2
sq_edges <- c()
for (i in 1:(N-1)) {
    for (j in (i+1):N) {
        if ((i+j) %in% Q) {
            sq_edges <- c(sq_edges, i, j)
    }
}

require(igraph)
sq_graph <- make_graph(sq_edges, directed=FALSE)
plot(sq_graph)

if (N == 23) {
    # does not find a path with start=1 ...
    hamiltonian(sq_edges, start=18, cycle=FALSE)
    # hamiltonian(sq_edges)                     # NULL
} else if (N == 32) {
    # the first of these graphs that is Hamiltonian ...
    # hamiltonian(sq_edges, cycle=FALSE)
    hamiltonian(sq_edges)
} else if (N == 45) {
    # takes much too long ...
    # hamiltonian(sq_edges, cycle=FALSE)
    hamiltonian(sq_edges)
}}
}
\keyword{ graph-algorithms }
