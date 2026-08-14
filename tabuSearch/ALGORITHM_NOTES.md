# Algorithm notes

## Upstream algorithm

The translation follows the binary tabu-search implementation in
`R/tabuSearch.R` from tabuSearch 1.2.0.  A one-bit move defines a neighbor.
Each preliminary search records the starting state and then performs
`iters-1` moves.  Non-improving moves enter a FIFO tabu list.  A tabu move can
be accepted when it exceeds both the best non-tabu neighbor and the aspiration
(best objective seen so far).  Intensification restarts from the best state;
diversification tabus the most frequently changed coordinates and starts from
a new random binary configuration.

## Intentional corrections

### Sampled-neighborhood utility

Upstream initializes the length-`size` neighbor-utility vector to zero, but
when `neigh < size` evaluates only the sampled coordinates.  The unevaluated
entries can therefore be selected as if their objective were exactly zero,
particularly when all actual objective values are negative.

The Fortran port tracks an explicit `evaluated` mask.  Only evaluated moves are
eligible.  If the random sample happens to contain only tabu moves and none is
admissible by aspiration, one non-tabu move is sampled and evaluated so the
search can continue.

### History used by diversification

The R implementation preallocates history storage substantially beyond the
number of completed iterations and then forms move frequencies over that whole
matrix.  The Fortran implementation computes frequencies using only the
`n_records` completed rows.

### Negative objectives during intensification

Upstream initializes the intensification sentinel `tempo` to zero.  This means
an all-negative objective can skip intensification even when a restart should
be attempted.  The Fortran port starts from `-huge(1.0_dp)`, making the logic
independent of the objective's sign.

### Diversification tie handling

The upstream code invokes random-tie `rank()` separately when constructing the
tabu indicator and the FIFO order, so tied move frequencies can produce two
different selected sets.  The Fortran port makes one random tie-breaking
ordering and uses the same selected set for both structures.

## RNG

R's global RNG is not reproduced.  `tabu_rng` is a local Park-Miller minimal
standard generator with explicit state.  A supplied `tabu_control%seed` makes
runs reproducible across Fortran compilers without relying on integer overflow.
This changes stochastic trajectories relative to an R run with `set.seed`, but
not the search mathematics.
