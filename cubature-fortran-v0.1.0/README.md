# cubature-fortran

Modern free-format Fortran translation of the computational interface of the
R package `cubature` 2.1.4-1.

The package provides adaptive multivariate integration over hyperrectangles
with scalar or vector-valued integrands:

- `hcubature` / `adaptIntegrate`: globally h-adaptive cubature using the
  15-point Gauss-Kronrod rule in one dimension and the embedded Genz-Malik
  degree-7/degree-5 rule in multiple dimensions.
- `pcubature`: p-adaptive tensor-product Clenshaw-Curtis integration.
- `cuhre`: deterministic Cuba-style adaptive cubature.
- `divonne`: adaptive stratified low-discrepancy integration.
- `suave`: adaptive randomized stratified integration.
- `vegas`: adaptive separable importance-sampling integration.
- `cubintegrate`: one uniform method-selection interface, including the
  tangent transform used by the R package for infinite limits.

All Fortran source is free-format `.f90`; implicit typing and implicit external
interfaces are disabled in `fpm.toml`.

## Example

```fortran
program demo
    use cubature, only : dp, cubature_result, hcubature
    implicit none
    type(cubature_result) :: ans

    call hcubature(f, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], &
        1, ans, rel_tol=1.0e-8_dp, abs_tol=1.0e-12_dp)
    print *, ans%integral(1), ans%error(1), ans%evaluations
contains
    subroutine f(x, value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value(:)
        value(1) = product(cos(x))
    end subroutine f
end program demo
```

The exact integral in this example is `sin(1)^2 = 0.7080734182735712...`.

## Build

```text
fpm build
fpm test
fpm run --example integrate_example
```

The release was also validated directly with GNU Fortran 14.2 using
Fortran 2018, optimization, warnings-as-errors, and full runtime checking.

## Result object

`type(cubature_result)` contains:

- `integral(:)` - integral estimate for each component;
- `error(:)` - estimated absolute error;
- `prob(:)` - reliability field corresponding to the Cuba result surface;
- `evaluations` - number of integrand evaluations;
- `return_code` - zero on convergence, nonzero otherwise;
- `nregions` - number of adaptive regions when applicable.

See `TRANSLATION_NOTES.md` for the precise correspondence with the two
upstream C libraries and for deliberate implementation differences.
