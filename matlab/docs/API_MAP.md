# API map

| R package routine | Fortran routine | Status / notes |
|---|---|---|
| `ceil` | `ceil` | Elemental real function. |
| `cell` | - | Omitted heterogeneous R container. |
| `colorbar` | - | Plotting omitted. |
| `eye` | `eye` | Rectangular or square real matrix. |
| `factors` | `factors` | Uses 64-bit integers. |
| `fileparts` | `fileparts` | Returns `type(fileparts_result)`. |
| `filesep` | `filesep` | Runtime platform separator. |
| `find` | `find` | Real, integer, and logical vectors/matrices; column-major linear indices. |
| `fix` | `fix` | Elemental truncation toward zero. |
| `fullfile` | `fullfile` | Two or three path components. |
| `hilb` | `hilb` | Real Hilbert matrix. |
| `imagesc` | - | Plotting omitted. |
| `isempty` | `isempty` | Assumed-rank, unlimited-polymorphic inquiry. |
| `isprime` | `isprime` | Elemental 64-bit integer function. |
| `jet.colors` | - | Plotting/color helper omitted. |
| `linspace` | `linspace` | Preserves `n < 2` behavior. |
| `logspace` | `logspace` | Preserves special `b == pi` behavior. |
| `magic` | `magic` | Odd, doubly-even, and singly-even orders. |
| `meshgrid` | `meshgrid`, `meshgrid3` | Typed 2-D and 3-D results. |
| `mod` | `matlab_mod` | Renamed to avoid collision with Fortran intrinsic `mod`. |
| `multiline.plot.colors` | - | Plotting/color helper omitted. |
| `ndims` | `ndims` | MATLAB-style vectors have shape `[1,n]`. |
| `nextpow2` | `nextpow2` | Elemental real-to-integer function. |
| `numel` | `numel` | Assumed-rank inquiry. |
| `ones` | `ones` | Rank-2 constructor. |
| `padarray` | `padarray` | Rank-2 real arrays; constant, circular, replicate, symmetric. |
| `pascal` | `pascal` | Supports `k=0,1,2`. |
| `pathsep` | `pathsep` | Runtime platform separator. |
| `pow2` | `pow2`, `pow2_scaled` | One- and two-argument forms. |
| `primes` | `primes` | Sieve of Eratosthenes. |
| `rem` | `rem` | MATLAB sign convention; division by zero returns NaN. |
| `repmat` | `repmat` | Real vectors/matrices; source-compatible branch by default. |
| `reshape` | `reshape2d` | Vector/matrix input to rank-2 output. |
| `rosser` | `rosser` | Exact integer 8x8 matrix. |
| `rot90` | `rot90` | Any integer number of quarter-turns. |
| `size` | `shape_of`, `size_dim` | Renamed to avoid shadowing Fortran intrinsic `size`. |
| `std` | `std`, `std_cols` | Sample SD by default; population flag is an extension. |
| `strcmp` | `strcmp` | Scalar character comparison. |
| `sum` | Fortran `sum`, `sum_cols` | Intrinsic `sum` handles vectors; `sum_cols` matches matrix behavior. |
| `tic`, `toc` | `tic`, `toc` | Uses `system_clock`. |
| `vander` | `vander`, `vander_complex` | Real and complex forms. |
| `zeros` | `zeros` | Rank-2 constructor. |
| `fliplr`, `flipud` | `fliplr`, `flipud` | Real vectors and matrices. |
