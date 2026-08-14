# Algorithm notes

## Unsigned 32-bit arithmetic

The original WELL implementation uses C `unsigned int`. Standard Fortran has
no portable unsigned 32-bit integer type, so this port stores each word in
`integer(int32)` and treats it as a bit pattern. Logical right shifts use
`shiftr`; left shifts use `shiftl`; XOR/AND/OR use the standard bit intrinsics.
This reproduces the recurrence independently of the sign bit.

Scalar seed expansion needs a multiplication modulo 2^32. That multiplication
is performed in `int64`, where the largest product fits, followed by explicit
reduction modulo 2^32. The implementation therefore does not depend on signed
integer overflow.

## Circular state

The upstream C sources contain specialized case functions to avoid `% R`
operations for non-power-of-two state lengths. This port expresses the same
recurrence directly with modular circular indexing. Cross-language validation
covers enough draws to traverse every state boundary multiple times.

## Tempering

`WELL19937c` is the tempered form of the WELL19937a recurrence.
`WELL44497b` is the tempered form of WELL44497a. The two post-recurrence
bit masks are reproduced exactly.

The R `WELL2test()` wrapper accepts `temper=TRUE` for several orders where its
C dispatcher actually ignores the flag. `well_from_options()` preserves those
observable choices: 19937 selects `c`, 44497 selects `b`, while 800/21701/23209
retain their selected ordinary recurrence.

## Modernization

Unlike the upstream implementation, generator state is held in a derived type
rather than file-scope globals. This makes independent simultaneous streams
natural and avoids process-global mutation while preserving the exact sequence
for a given variant, seed, and state.
