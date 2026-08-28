# Algorithm notes

## Sobol and Owen scrambling

The ordinary Sobol constants are translated from the package's bundled
Joe-Kuo/PBRT tables. The extended Owen-scrambled implementation retains all
21201 x 32 direction words from `sobol_directions.h`. Fast Owen scrambling
uses the same 32-bit bit reversal, hash, xor, multiplication, and shift
operations as upstream.

The large tables are represented by free-form Fortran module arrays initialized
with per-dimension `DATA` statements. This avoids gfortran's constructor and
continuation limits without changing any table bits.

## PCG32

The translated generator implements the upstream default-stream PCG32
set-sequence XSH-RR generator. Exact modulo-2^64 state arithmetic is performed
using a processor-supported integer kind with at least 38 decimal digits and
then reduced explicitly to 64 bits, avoiding reliance on signed integer
overflow.

## Halton

Faure digit permutations are generated recursively. The generated upstream
Halton kernels use a fixed number of reversed base digits and a single-
precision scale; the Fortran implementation reproduces those details, including
all per-prime digit counts through dimension 256.

Random Halton differs only in permutation generation: upstream uses
`std::shuffle`, whose seeded result is standard-library dependent. The Fortran
port specifies Fisher-Yates/PCG32 shuffling for portable reproducibility.

The upstream public C++ wrappers for Halton single values reverse the sampler
arguments. The port fixes this and exposes the intended `(dimension,index)`
low-level semantics.

## PJ / PMJ / PMJ(0,2)

PJ, PMJ, PMJ blue-noise, PMJ(0,2), and PMJ(0,2) blue-noise are direct
translations of the bundled progressive-sampling implementations. Blue-noise
variants use 100 best candidates and toroidal nearest-neighbor distance.

A subtle upstream C++ behavior is intentionally preserved: the sample-set
objects copy the PCG generator, while subquadrant selection continues to use
the caller's generator. Therefore jitter/candidate generation and subquadrant
selection advance two independent PCG streams initialized from the same seed.
This is required for exact sequence compatibility.
