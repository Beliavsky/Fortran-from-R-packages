# API mapping

| Upstream R API | Fortran API | Notes |
| --- | --- | --- |
| `dtruncnorm(x,a,b,mean,sd)` | `dtruncnorm` | Scalar or vector `x`; `dtruncnorm_recycle` for fully recycled arrays. |
| `ptruncnorm(q,a,b,mean,sd)` | `ptruncnorm` | Scalar or vector `q`; recycled-array form supplied. |
| `qtruncnorm(p,a,b,mean,sd)` | `qtruncnorm` | Uses the translated R/NETLIB `zeroin` logic; vector and recycled forms supplied. |
| `rtruncnorm(n,a,b,mean,sd)` | `rtruncnorm(n,...)` | Returns an allocatable vector. `rtruncnorm_recycle` reproduces recycling of parameters. |
| `etruncnorm(a,b,mean,sd)` | `etruncnorm` | Scalar; `etruncnorm_recycle` for array inputs. |
| `vtruncnorm(a,b,mean,sd)` | `vtruncnorm` | Scalar; `vtruncnorm_recycle` for array inputs. |

A scalar one-draw overload of `rtruncnorm(a,b,mean,sd)` is also supplied for
Fortran convenience.

## Internal native code

The following upstream computational support is incorporated into
`truncnorm_core` rather than exposed as R-style native entry points:

- exponential rejection sampling on `(a,Inf)` and `(a,b)`;
- ordinary-normal, half-normal and uniform rejection samplers;
- left/right/two-sided sampler dispatch using `t1=0.15`, `t2=2.18`,
  `t3=0.725`, `t4=0.45`;
- R/NETLIB `truncnorm_zeroin` Brent root finding;
- left/right/two-sided expectation and variance helpers.

`exports.c` and `sexp_macros.h` are R interface glue and are retained only in
`upstream/`.
