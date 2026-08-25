# Why Fortran rather than C++?

This repository translates computational code from R packages into modern
Fortran. C++ would also be a reasonable implementation language, particularly
because Rcpp, pybind11, and MATLAB's C++ interfaces make language integration
convenient. The choice of Fortran is an engineering tradeoff rather than a
claim that Fortran is universally preferable.

## Comparison

| Area | Modern Fortran | C++ |
| --- | --- | --- |
| Numerical arrays | Native multidimensional arrays and array operations | Usually uses loops or an additional array library |
| Similarity to vectorized R | Mathematical expressions often translate directly | Commonly requires more structural rewriting |
| Performance | Excellent for dense numerical kernels | Equally capable when algorithms and memory access are comparable |
| R integration | Requires an interoperability layer such as `bind(C)`, `.Call()`, or RFI | Rcpp provides a mature high-level interface |
| Python integration | Uses a C ABI, f2py, Cython, `ctypes`, or another extension layer | pybind11 and nanobind provide mature interfaces |
| MATLAB and Octave integration | Uses C-compatible Fortran wrappers and gateway code | C and C++ MEX interfaces are directly supported |
| Build and packaging | fpm is simple but relatively young | CMake and the wider C++ packaging ecosystem are more established |
| Numerical libraries | Direct access to BLAS, LAPACK, and many established numerical codes | Broad access to Eigen, Boost, and other modern libraries |
| Complex data structures | Best suited to regular numeric arrays and procedural kernels | Strong support for containers, graphs, trees, strings, and object models |
| Generic programming | Available but comparatively limited | Templates and concepts are powerful, though more complex |
| Portability across languages | A stable C ABI can be provided with `bind(C)` | A C ABI is also preferable because the native C++ ABI is not universally stable |
| Tool availability | Fortran toolchains can present an additional installation hurdle | C++ compilers and development tools are more widely available |

Neither language has an inherent speed advantage for every program. Algorithm
selection, data layout, allocation, compiler optimization, and library quality
usually matter more than the implementation language.

## Why Fortran fits these translations

Most packages in this repository are dominated by statistics, probability
distributions, optimization, time-series analysis, and matrix computation.
These workloads map naturally to Fortran because it offers:

- Native multidimensional arrays and concise whole-array expressions.
- Column-major storage compatible with R, MATLAB, and Octave.
- Strong optimization of numerical loops and array operations.
- Straightforward use of BLAS and LAPACK.
- Compact procedural implementations that remain close to the mathematical
  structure of the original R code.

Shared R-like numerical helpers can also be placed in `r_mod.f90`, reducing
the need for each translation to implement basic statistical operations
independently.

## Where C++ has advantages

C++ is often a better choice when the central challenge is integration or
complex data organization rather than numerical computation. Examples include
code dominated by graphs, trees, hash tables, strings, callbacks, sparse-data
frameworks, or rich model objects. Its established R, Python, and MATLAB
binding ecosystems can also simplify distribution.

Those advantages do not require the numerical translations to be rewritten.
C or C++ can be used in small language-specific adapters around a common
Fortran interface:

```text
R / Python / MATLAB / Octave
             |
      language adapters
             |
      stable bind(C) API
             |
 existing Fortran translation
```

The public interoperability layer should use C-compatible scalar types,
explicit dimensions and shapes, caller-provided output storage, documented
column-major layout, and status codes rather than allowing Fortran termination
or C++ exceptions to cross the boundary.

## Project direction

Fortran remains the primary language for translated computational kernels.
Thin C or C++ adapters may be added where they improve integration with R,
Python, MATLAB, or Octave. A new translation may reasonably use C++ instead
when its main requirements are complex data structures, callbacks, or
language integration rather than array-oriented numerical computation.
