# API inventory

This file lists the public procedures by module.  For exact argument ranks, optional arguments, and result types, use the explicit interfaces in `src/`.  All real-valued numerical work uses `dp = kind(1.0d0)` unless a procedure is explicitly integer or complex.

## Aggregate use

```fortran
use pracma
```

The aggregate module re-exports the public entities of all component modules.

## Result types

`root_result`, `optimization_result`, `quadrature_result`, `ode_result`, `linear_solve_result`, `eigen_result`, `symmetric_eigen_result`, `polynomial_division_result`, `peak_result`, `circle_result`, `pchip_result`, `qp_result`, `polynomial_fit_result`, `regression_result`, `qpspecial_result`, `linprog_result`.

## Status codes

`pracma_ok`, `pracma_invalid_argument`, `pracma_dimension_mismatch`, `pracma_singular`, `pracma_not_converged`, `pracma_nonfinite`, and related procedure-specific statuses are returned through result fields or optional `status` arguments.

## `pracma_callbacks` - Callback interfaces

- `scalar_function`, `scalar_derivative`, `objective_function`, `vector_function`, `vector_field`, `complex_scalar_function`, `complex_objective_function`, `complex_vector_function`
- `bivariate_function`, `trivariate_function`, `vector_curve`, `regression_model`

## `pracma_basic` - Array, elementary, and statistical utilities

- `linspace`, `logspace`, `logseq`, `eye`, `ones`, `zeros`, `Diag`, `blkdiag`
- `hilb`, `pascal`, `vander`, `hankel`, `Toeplitz`, `flipud`, `fliplr`, `rot90`
- `circshift`, `repmat`, `meshgrid`, `size2`, `numel`, `nnz`, `isempty`, `dot`
- `cross`, `crossn`, `Norm`, `fnorm`, `Trace`, `harmmean`, `geomean`, `trimmean`
- `Mode`, `std`, `std_err`, `ceil`, `Fix`, `mod`, `rem`, `idivide`
- `gcd`, `Lcm`, `nextpow2`, `pow2`, `distmat`, `pdist`, `pdist2`, `hausdorff_dist`
- `accumarray`, `bsxfun`, `uniq`, `find`, `finds`, `findintervals`, `deg2rad`, `rad2deg`
- `sind`, `cosd`, `tand`, `cotd`, `asind`, `acosd`, `atand`, `acotd`
- `atan2d`, `secd`, `cscd`, `asecd`, `acscd`, `cot`, `csc`, `sec`
- `acot`, `acsc`, `asec`, `coth`, `csch`, `sech`, `acoth`, `acsch`
- `asech`, `sigmoid`, `logit`, `hypot_pracma`, `eps`, `real_part`, `imag_part`, `angle`
- `sort_real`

## `pracma_linalg` - Linear algebra and matrix functions

- `outer_product`, `identity_matrix`, `symmetrize`, `solve_linear`, `inverse_matrix`, `inv`, `mldivide`, `mrdivide`
- `determinant`, `logdet_spd`, `cholesky`, `isposdef`, `nearest_spd`, `symmetric_eigen`, `eigjacobi`, `eig`
- `pinv`, `Rank`, `nullspace`, `cond`, `rref`, `gramSchmidt`, `qrSolve`, `givens`
- `householder`, `lu`, `lu_crout`, `lufact`, `lusys`, `expm`, `logm`, `sqrtm`
- `signm`, `rootm`, `hessenberg`, `arnoldi`, `gmres`, `orth`, `subspace`, `trace_matrix`
- `frobenius_norm`, `charpoly`, `linearproj`, `affineproj`, `procrustes`, `kabsch`, `nearest_symmetric`

## `pracma_polynomial` - Polynomial and approximation routines

- `polyval`, `polyvalm`, `polyadd`, `polymul`, `polypow`, `polydiv`, `polyder`, `polyint`
- `polytrans`, `poly2str`, `roots`, `rootsmult`, `polyroots`, `compan`, `horner`, `hornerdefl`
- `polyfit`, `polyfix`, `chebPoly`, `chebCoeff`, `chebApprox`, `legendre`, `laguerre`, `bernstein`
- `bernsteinb`, `pade`, `Poly`, `polygcf`, `rationalfit`, `trigPoly`, `trigApprox`, `polyval_vector`

## `pracma_special` - Special mathematical functions

- `sinc`, `psinc`, `agmean`, `gammaz`, `gammainc`, `incgam`, `psi`, `erf_pracma`
- `erfc_pracma`, `erfinv`, `erfcinv`, `erfcx`, `erfz`, `erfi`, `expint`, `expint_Ei`
- `li`, `Si`, `Ci`, `fresnel`, `lambertWp`, `lambertWn`, `zeta`, `eta`
- `polylog`, `ellipke`, `ellipj`, `nthroot`, `einsteinF`, `bernoulli`, `factorial2`, `golden_ratio`
- `humps`

## `pracma_differentiation` - Numerical differentiation

- `gradient`, `grad`, `jacobian`, `hessian`, `hessvec`, `hessdiag`, `laplacian`, `numderiv`
- `numdiff`, `complexstep`, `grad_csd`, `jacobian_csd`, `hessian_csd`, `laplacian_csd`, `fornberg`, `fderiv`
- `taylor`

## `pracma_integration` - Numerical integration

- `trapz`, `cumtrapz`, `trapzfun`, `midpoint`, `simpson`, `romberg`, `quad`, `quadl`
- `quadgk`, `quadcc`, `quadgr`, `gauss_kronrod`, `integral`, `simpadpt`, `quadv`, `clenshaw_curtis`
- `cotes`, `gaussLegendre`, `gaussHermite`, `gaussLaguerre`, `rectint`, `simpson2d`, `triquad`, `quad2d`
- `dblquad`, `triplequad`, `integral2`, `integral3`, `line_integral`, `quadinf`

## `pracma_roots` - Root finding and nonlinear systems

- `bisect`, `secant`, `regulaFalsi`, `brentDekker`, `fzero`, `newtonRaphson`, `halley`, `newtonHorner`
- `muller`, `ridders`, `findzeros`, `findmins`, `fsolve`, `fzsolve`, `newtonsys`, `broyden`
- `itersolve`, `aitken`

## `pracma_optimization` - Optimization and least squares

- `fminbnd`, `fibsearch`, `nelder_mead`, `fminsearch`, `anms`, `hooke_jeeves`, `steep_descent`, `fminunc`
- `fletcher_powell`, `gaussNewton`, `curvefit`, `lsqnonlin`, `lsqcurvefit`, `lsqnonneg`, `L1linreg`, `geo_median`
- `qpspecial`, `qpsolve`, `quadprog`, `linprog`, `fmincon`, `qpsolve_projection`

## `quadprog` - Quadratic programming

- `solve_qp`, `solve_qp_compact`

## `pracma_interpolation` - Interpolation and splines

- `interp1`, `interp1_linear`, `interp1_nearest`, `interp1_spline`, `interp1_pchip`, `interp2`, `neville`, `newton_interp`
- `lagrange_interp`, `barycentric`, `barycentric_weights`, `barylag`, `barylag2d`, `cubicspline`, `pchip`, `akima`
- `ppval`, `mkpp`, `deval`, `ratinterp`, `lebesgue_function`, `spinterp`, `cutpoints`, `piecewise_linear`

## `pracma_ode` - ODE and time-stepping methods

- `rk4`, `rk4sys`, `euler_heun`, `abm3`, `rkf54`, `ode23`, `ode45`, `ode78`
- `ode23s`, `bulirsch_stoer`, `newmark`, `cranknic`, `shooting_ivp`

## `pracma_signal_stats` - Signal processing and statistics

- `conv`, `deconv`, `fft`, `ifft`, `fftshift`, `ifftshift`, `detrend`, `movavg`
- `savgol`, `savgol_filter`, `findpeaks`, `hampel`, `entropy`, `hurst`, `histc`, `rmserr`
- `runge`, `humps_function`, `mexpfit`, `odregress`, `trigregress`, `whittaker`, `andrews_curve`, `normest`
- `autocorrelation`, `crosscorrelation`, `moving_median`, `moving_mean`, `periodogram`

## `pracma_geometry` - Geometry and spatial routines

- `cart2sph`, `sph2cart`, `cart2pol`, `pol2cart`, `haversine`, `inpolygon`, `polyarea`, `polycenter`
- `poly_length`, `polygon_crossings`, `circlefit`, `segment_intersection`, `segment_distance`, `point_segment_distance`, `arclength`, `stereographic_project`
- `stereographic_inverse`, `fractalcurve`, `triangle_area`, `triarea`, `line_intersection`, `nearest_point_polyline`, `poisson2disk`, `kriging`
- `project_coordinates`, `plane_projection`

## `pracma_combinatorics` - Combinatorics, number theory, and matrix generators

- `isprime`, `primes`, `factors`, `fact`, `factorial`, `nchoosek`, `perms`, `combs`
- `randperm`, `randcomb`, `dec2bin`, `bin2dec`, `bitget`, `bitset`, `bitand_vec`, `bitor_vec`
- `bitxor_vec`, `hadamard`, `magic`, `moler`, `rosser`, `wilkinson`, `randortho`, `set_random_seed`
- `rand_uniform`, `randn`, `randi`, `sample_without_replacement`, `sumalt`, `andor`, `modular_power`, `chinese_remainder`

## `pracma_compat` - Compatibility names and supplemental wrappers

- `invlap`, `rand`, `rand_matrix`, `randn_matrix`, `rands`, `randp`, `randsample`, `bits`
- `akimaInterp`, `ppfit`, `lsqlin`, `lsqlincon`, `kron`, `trisolve`, `segm_intersect`, `segm_distance`
- `poly_crossings`, `tril`, `triu`, `squareform`, `flipdim`, `softline`, `piecewise`, `divisors`
- `rssimple`, `hurstexp`, `approx_entropy`, `sample_entropy`, `bubbleSort`, `insertionSort`, `selectionSort`, `shellSort`
- `heapSort`, `mergeSort`, `quickSort`, `quickSortx`, `mergeOrdered`, `is_sorted`, `sortrows`, `normF`
- `abm3pc`, `shubert`, `bvp`, `stereographic`, `stereographic_inv`, `lebesgue`, `half`

## Conventions

- R periods in procedure names generally become underscores.
- Matrices use ordinary Fortran column-major storage.
- Procedures do not perform R-style vector recycling.
- Callback procedures must match the abstract interfaces in `pracma_callbacks`.
- Failure is reported through status codes or typed result fields rather than R conditions.
