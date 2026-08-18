# tensorA-fortran

A free-format modern Fortran translation of the computational core of the R
package `tensorA` 0.36.2.1 (K. Gerald van den Boogaart).

`tensorA` provides advanced tensor arithmetic with named indices, including
Einstein and Riemann summation, covariance/contravariance operations,
batched linear algebra, and tensor statistics.

## Design

The Fortran representation is `type(tensor_t)`. A tensor stores:

- complex double-precision data in Fortran column-major order;
- an integer shape;
- a name for every tensor axis.

Real tensors are represented as complex tensors with zero imaginary part so
one implementation supports both the real and complex paths from the original
R/C code.

Axis names are computational: named contraction, broadcasting, reordering,
Einstein summation, and Riemann summation all use them. Duplicate axis names
are allowed when an operation naturally creates them, but a later name lookup
fails explicitly when the requested name is ambiguous.

## Main translated functionality

- tensor construction and conversion helpers;
- axis positions, renaming, marking, reordering, repetition, slicing, binding,
  dropping/undropping and axis collapsing;
- named and positional tensor contraction (`mul.tensor` equivalent);
- trace, margins, diagonal multiplication, delta, diagonal, triple-delta and
  one tensors;
- named broadcasting for addition, subtraction, multiplication and division;
- Einstein and Riemann pair contraction;
- covariant/contravariant name conversion and metric index dragging;
- inverse and Moore-Penrose pseudoinverse;
- linear solve with batch (`by`) dimensions;
- complex-capable SVD, Cholesky, tensor powers and operator norms;
- matrix reshaping;
- tensor norm, means, covariance and variance.

The contraction kernel follows the performance intent of the original C
`tensoraCmulhelper`/`tensoramulhelper`: tensors are reordered into batched
matrix blocks and Fortran `matmul` is used for each `by` slice. No C code is
called.

## Build

```sh
fpm build
fpm test
fpm run --example tensor_example
```

The library requires only a Fortran 2018 compiler. It has no BLAS, LAPACK, C,
or C++ dependency. The SVD, inverse/pseudoinverse, and Cholesky support are
implemented in native Fortran.

## Small example

```fortran
program demo
    use tensora
    implicit none
    type(tensor_t) :: a, b, c

    a = tensor([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp], &
               [2,3], ['i','j'])
    b = tensor([1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp], &
               [3,2], ['j','k'])
    c = einstein_pair(a,b)

    print *, real(c%data,dp)
end program demo
```

See `TRANSLATION_NOTES.md` for the R-to-Fortran mapping and intentional
interface differences.
