# API mapping

| R caRamel 1.5 | Fortran | Notes |
|---|---|---|
| `caRamel()` | `caramel_optimize()` | Sequential native callback API; no R cluster/UI/file-output plumbing |
| `pareto()` | `pareto()` | Integer 0/1 front flags; maximization convention |
| C `c_pareto_2d` / Fortran `pareto_2d` | `pareto_2d()` | Public compatibility entry point |
| C `c_pareto_3d` / Fortran `pareto_3d` | `pareto_3d()` | Public compatibility entry point |
| `dominate()` | `dominate()` | Onion-peel Pareto rank |
| `dominated()` | `dominated()` | Logical mask |
| `val2rank()` | `val2rank()` | `opt=1/2/3` semantics retained |
| `boxes()` | `boxes()` | Returns 64-bit integer box identifiers |
| `downsize()` | `downsize()` | One retained representative per precision box |
| `decrease_pop()` | `decrease_pop()` | Archive/population index selection |
| `vol_splx()` | `vol_splx()` | Native determinant calculation |
| `Dimprove()` | `dimprove()` | Oriented improvement edges and lengths |
| `rselect()` | `rselect()` | Weighted sampling with replacement |
| `matvcov()` | `matvcov()` | Upstream standardized covariance convention |
| `Cusecovar()` | `cusecovar()` | Native normal RNG + Cholesky |
| `Cinterp()` | `cinterp()` | Same normalized-uniform barycentric rule as R |
| `Cextrap()` | `cextrap()` | Same exponential extrapolation length rule |
| `Crecombination()` | `crecombination()` | `index_block` represents functional blocks |
| `newXval()` | `new_xval()` | All four generation rules + fireworks rule |
| `geometry::delaunayn()` | `delaunay_nd()` | Native n-D Bowyer-Watson-style incremental triangulation |
| `plot_caramel()` | omitted | Plotting code |
| `plot_pareto()` | omitted | Plotting code |
| `plot_population()` | omitted | Plotting code |
| `.onAttach()` | omitted | R package UI only |

## Main callback difference

The R package's sequential objective function receives an integer row index and reads a global matrix named `x`. The Fortran callback instead receives the parameter vector directly:

```fortran
subroutine objective(x, values)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: values(:)
end subroutine objective
```

This removes global mutable state and makes the computational interface usable from ordinary Fortran applications.
