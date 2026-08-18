# Translation notes

## Scope

`MultiRNG` 1.2.4 contains 13 R source files and no compiled C/C++/Fortran backend.
The full computational surface has been translated to native free-format Fortran.
There was no plotting or graphical code to omit.

## Mapping

| R routine | Fortran routine |
|---|---|
| `draw.d.variate.normal` | `draw_d_variate_normal` |
| `draw.d.variate.t` | `draw_d_variate_t` |
| `draw.d.variate.uniform` | `draw_d_variate_uniform` |
| `draw.correlated.binary` | `draw_correlated_binary` |
| `draw.dirichlet` | `draw_dirichlet` |
| `draw.multinomial` | `draw_multinomial` |
| `draw.dirichlet.multinomial` | `draw_dirichlet_multinomial` |
| `draw.multivariate.hypergeometric` | `draw_multivariate_hypergeometric` |
| `draw.multivariate.laplace` | `draw_multivariate_laplace` |
| `draw.wishart` | `draw_wishart`, `draw_wishart_flat` |
| `draw.inv.wishart` | `draw_inv_wishart`, `draw_inv_wishart_legacy` |
| `generate.point.in.sphere` | `generate_point_in_sphere` |
| `loc.min` | `loc_min` |

## Numerical implementation

The Fortran package is self-contained. It supplies native implementations of normal,
gamma/chi-square, Poisson, binomial, and hypergeometric random generation, Cholesky
factorization, SPD inversion, and the normal CDF needed by the Gaussian-copula uniform
generator.

The correlated-binary generator retains the Park-Park-Shin common-Poisson construction.
The R matrix-indexing state machine is expressed as an explicit clique/incidence
decomposition of the same `alpha` overlap matrix, avoiding R's chained subassignment
semantics while retaining the target means and correlations.

## Upstream behaviors deliberately preserved

1. **Multivariate t**: the upstream function draws one chi-square random variable for a
   complete call and applies that common scale to every generated row. This makes rows
   dependent within a call. The Fortran translation preserves that behavior.

2. **Dirichlet-multinomial**: the upstream function generates `no.row` Dirichlet draws,
   averages them columnwise, then uses the resulting single probability vector for all
   multinomial rows. It therefore is not the usual independent Dirichlet mixture per
   row. The Fortran routine preserves the upstream algorithm.

3. **Wishart output**: upstream flattens each matrix into a row. The primary Fortran
   routine returns `(replicate,d,d)`; `draw_wishart_flat` and
   `draw_inv_wishart_flat` provide row-flattened compatibility output.

## Corrected inverse-Wishart routine

The upstream `draw.inv.wishart()` contains a computational defect: after producing a
Wishart draw it never takes the matrix inverse required by the inverse-Wishart law.
`draw_inv_wishart()` implements the density documented in `draw.inv.wishart.Rd`:
if `inv_sigma = Sigma^{-1}`, it draws

    W ~ Wishart(nu, inv_sigma)
    X = W^{-1}

so `X ~ InvWishart(nu, Sigma)`.

`draw_inv_wishart_legacy()` reproduces the original R implementation exactly for users
who need behavioral compatibility with MultiRNG 1.2.4.

## Validation

The tests check:

- multivariate-normal means and covariance;
- Gaussian-copula uniform margins;
- Dirichlet, multinomial, Dirichlet-multinomial, and multivariate-hypergeometric means;
- Wishart and inverse-Wishart theoretical means;
- unit-sphere geometry and the `gamma=2` multivariate-Laplace covariance;
- multivariate-t means and `loc_min`;
- the 3-D correlated-Bernoulli example's requested means and correlations.
