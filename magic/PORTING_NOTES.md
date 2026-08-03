# Porting notes

## Array representation

R permits arrays of runtime rank. Standard Fortran allocatable arrays have a
compile-time rank, so this port uses `integer_tensor`, a flat column-major
array with a runtime `shape(:)`. Square-specific routines retain ordinary
rank-two arrays for convenient `matmul`, slicing, and printing.

## Indexing

The public API is one-based, matching both R and Fortran. `tensor_offset` and
`unravel_index` preserve Fortran/R column-major ordering. Circular operations
use mathematical modulo so negative shifts behave like R's `%%`-based code.

## Integer range

The upstream package mostly uses R integers and doubles. This port uses signed
64-bit integers throughout. Very large orders can still overflow products,
magic constants, or multiplicative tests; callers should choose orders that
fit in `integer(ik)` and available memory.

## Random Latin squares

The incidence-array Markov move is translated directly. Reproducibility is
controlled with the standard Fortran `random_seed` intrinsic. The exact stream
will not match R's RNG.

## Error handling

Invalid-input procedures can receive an optional `magic_error`. Constructors
return an unallocated or zero-sized result when an error prevents creation.
The examples use valid inputs and therefore omit the status object.
