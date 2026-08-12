# Algorithm notes

## Insertion heuristics

The port follows the upstream Rosenkrantz-style construction routines. For
ATSP input, nearest/farthest candidate selection considers both the distance to
and from already placed cities. Cheapest insertion uses the translated C
insertion-delta kernel.

`arbitrary_insertion` follows the optimized upstream strategy: randomize the
city order once, seed the partial tour with the first two cities, and insert
remaining cities at their cheapest positions.

## 2-opt

The algorithm supports asymmetric costs by including the direction changes of
all internal edges when a subtour is reversed. The upstream closing-edge bug
described in `README.md` is fixed.

Only improvements larger than `1e-7`, matching the upstream C constant, are
accepted.

## Simulated annealing

The default local move is subtour reversal. Swap and mixed modes are also
provided. The initial temperature defaults to initial-tour length divided by
number of cities, and the temperature uses a logarithmic schedule consistent
with R's SANN documentation/control semantics.

## Infinite distances

`solve_tsp` replaces infinities before heuristics exactly in the spirit of
upstream `.replaceInf`: positive infinity becomes `max(finite)+2*range`, and
negative infinity becomes `min(finite)-2*range`. The final reported length is
computed on the original matrix, so a solution using a prohibited `Inf` edge
still reports an infinite length.

## TSPLIB

Reading supports the upstream explicit formats:
`FULL_MATRIX`, `UPPER_ROW`, `LOWER_COL`, `UPPER_COL`, `LOWER_ROW`,
`UPPER_DIAG_ROW`, `LOWER_DIAG_COL`, `UPPER_DIAG_COL`, and `LOWER_DIAG_ROW`,
plus coordinate types `EUC_2D`, `EUC_3D`, `ATT`, and `GEO`.
