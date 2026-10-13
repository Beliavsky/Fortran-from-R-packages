# Numerical API

Import `ineq_api`, or `ineq_indices` and `ineq_lorenz`. `dp` is re-exported.
All sample inputs are rank-one `real(dp)` arrays. All index results are scalar
`real(dp)`. All computational procedures are pure. Optional `na_rm` defaults
true. Missing observations are IEEE NaNs; NA and NaN are not distinguished.
Theoretical scalar procedures `lc_exp`, `lc_lognorm`, and `lc_pareto` are elemental.

## Scalar indices

| Procedure | Optional arguments | Default/result details |
|---|---|---|
| `ineq(x,parameter,type,na_rm)` | `parameter`, character `type`, `na_rm` | Type Gini by default; exact R names only |
| `gini(x,corr,na_rm)` | logical `corr`, `na_rm` | Correction false; when true multiply by n/(n-1) |
| `rs(x,na_rm)` | `na_rm` | Mean absolute deviation divided by twice the mean |
| `atkinson(x,parameter,na_rm)` | `parameter=0.5`, `na_rm` | Parameter one selects geometric mean |
| `theil(x,parameter,na_rm)` | `parameter=0`, `na_rm` | Zero is Theil T; any other finite parameter is mean log deviation |
| `kolm(x,parameter,na_rm)` | `parameter=1`, `na_rm` | Same units as income; zero parameter gives NaN |
| `var_coeff(x,square,na_rm)` | logical `square=.false.`, `na_rm` | Population SD/mean; singleton gives NaN as R's var is undefined |
| `entropy(x,parameter,na_rm)` | `parameter=0.5`, `na_rm` | At 0/1 dispatches to upstream Theil L/T |
| `conc(x,parameter,type,na_rm)` | `parameter`, `type='Herfindahl'`, `na_rm` | Herfindahl or Rosenbluth |
| `herfindahl(x,parameter,na_rm)` | `parameter=1`, `na_rm` | Generalized concentration, order parameter+1 |
| `rosenbluth(x,na_rm)` | `na_rm` | Rank concentration |
| `pov(x,k,parameter,type,na_rm)` | `parameter`, `type='Watts'`, `na_rm` | Watts, Sen, SST, or Foster |
| `watts(x,k,na_rm)` | `na_rm` | Average log poverty gap |
| `sen(x,k,na_rm)` | `na_rm` | Sen poverty |
| `sst(x,k,na_rm)` | `na_rm` | Sen-Shorrocks-Thon poverty |
| `foster(x,k,parameter,na_rm)` | `parameter=1`, `na_rm` | Power is parameter-1; parameters 1/2/3 give headcount/gap/squared gap |

The dispatchers accept the exact strings in the upstream namespace/examples.
`ineq` additionally accepts `var` and `square.var`. Unknown names yield NaN.
Poverty uses the strict comparison `x < k`; `k` must be finite. An empty poor
subsample returns zero, including empty retained inputs. Other empty sample
indices return NaN except Herfindahl (zero at nonnegative parameter, positive
infinity at negative parameter) and Rosenbluth (-1), matching R. Parameters are
finite unless explicitly noted for lognormal spread below. An exact zero
Herfindahl parameter returns one for a nonzero-total sample, retaining R's
`1^Inf` convention; its nonzero-parameter limit generally differs.

No general nonnegative-income restriction is imposed. Valid finite signed data
follow the numerical formulas, including signed coefficients or NaNs where real
logarithms/fractional powers are undefined. Zero-total Gini/RS/CV can be infinite;
Rosenbluth and Herfindahl zero totals are NaN. Infinities in empirical inputs are
rejected, including when `na_rm` is true. Missing-data deletion applies only to
NaNs. Nonfinite parameters produce NaN.

## Empirical curves and asymmetry

`curve = lc(x, n)` returns `type(lorenz_curve)`. Optional `n` is a real frequency
vector of length `size(x)`, defaulting to ones. A missing `x(i)` removes both
`x(i)` and `n(i)`. Stable ascending sorting retains corresponding frequencies.
`curve%p`, `%l`, and `%l_general` are allocatable vectors, beginning at the
origin. `p` is cumulative frequency/total frequency; `l` is cumulative
income/total income; `l_general` is cumulative income/total frequency. Frequencies
may be real (not just integer). Signed frequencies retain raw algebra but need
not describe a monotone Lorenz curve; application-domain validation is the
caller's responsibility. Zero frequencies are allowed in `lc`.

`curve = lc_mehran(x,n)` returns the same type, with `%p` and `%l` only;
`%l_general` has length zero. It applies the source's tangent slopes,
intersections, clipping and finite-abscissa filtering. At least two retained
classes are required. Zero-frequency classes or undefined slopes can make the
bound invalid. As in R, successful output can contain repeated abscissae; do not
silently deduplicate them. This is a bound computed from grouped means, not an
interpolated estimate of within-group data.

`curve%status` is `ineq_ok=0` or `ineq_invalid=1`. Invalid results have empty
coordinate arrays. This deliberately replaces R errors or partially NaN curves
for empty input, zero total frequency/income, invalid frequency length, or
nonfinite retained values. Coordinate arithmetic itself can overflow when the
mathematical result is not representable.

`values = lasym(x,n,interval,na_rm)` returns an allocatable vector. Optional
`n` defaults to ones; `interval` defaults false. When observations equal the
weighted mean, `interval=.true.` returns two bounds; otherwise even a requested
interval contains one value. With `na_rm=.false.` and a missing observation,
it returns one or two NaNs according to `interval`.

Two unusual upstream Lasym behaviors are preserved deliberately:

1. After removing missing x, the source indexes the original n with the new
   shorter permutation. Thus the first retained-count frequencies are used,
   rather than deleting their paired missing entries. `lc` behaves differently.
2. Interpolation between values surrounding the mean adds `delta` as one
   frequency unit, not `delta` times the adjacent class frequency.

The exact-mean branch uses original-scale arithmetic where representable, so
ordinary exact ties are not lost to normalization. Undefined inputs yield NaN
vectors. Frequency shape mismatch or infinite/NaN retained weights is invalid.

## Theoretical curves

All theoretical curves evaluate a population share `p` in [0,1]; outside this
range they return NaN. Call array-parameter functions in a loop for many p values.

| Procedure | Default parameter | Formula |
|---|---|---|
| `lc_exp(p)` | None | p + (1-p) log(1-p), with endpoint one filled in |
| `lc_lognorm(p,parameter)` | 1 | Phi(Phi-inverse(p) - parameter) |
| `lc_pareto(p,parameter)` | 2 | 1 - (1-p)^((parameter-1)/parameter) |
| `lc_singh(p,parameter)` | [2,2] | I(1-(1-p)^b; 1+d,b-d), b=parameter(1)-1, d=1/(parameter(2)+1) |
| `lc_dagum(p,parameter)` | [2,2] | I(p^b; a+1/b,1-a), a=1/parameter(1), b=parameter(2) |
| `theor_lc(p,type,parameter)` | Type Singh-Maddala | Direct evaluation in place of an R closure |

For `theor_lc`, accepted exact type names are `Singh-Maddala`, `Dagum`, `lognorm`,
`Pareto`, and `exponential`. The parameter is a vector even for one-parameter
families; only the first required components are used. An omitted or too-short
two-parameter vector falls back to [2,2]. Nonpositive lognormal spread uses one;
Pareto shape at most one uses two. Infinite lognormal spread is supported as
in R (zero in the interior and unchanged endpoints); its NaN parameter is NaN.
Infinite Pareto parameters yield NaN. The original Singh/Dagum parameter guards
are preserved for finite inputs; they are weaker than validity of the resulting
beta shapes. Negative beta shapes return NaN. R zero-shape limits, including
some economically inadmissible but algebraically defined parameters, are handled.
Missing/nonfinite components of two-parameter vectors fall back to [2,2], unlike
R cases that error or apply special infinite-shape limits. No warning stream is
produced for a fallback.

`lc_exp(1)=1` is the only deliberate change in the ordinary valid-probability,
valid-parameter formulas: upstream returns NaN at that endpoint. Stable
`log1p`/`expm1` transformations and the exponential small-p series can also give
a nonzero tail where upstream
subtraction rounds to zero.

## Majorization

```fortran
call major(x, y, dominates, status)
```

`dominates` is logical output. Optional integer `status` is zero for comparable
inputs, `ineq_incomparable=2` for unequal retained lengths or exact totals, and
`ineq_invalid=1` for nonfinite retained data. NaNs are removed, matching R's
`sort` default. The result is true if every ascending cumulative sum of x is
at most that of y; the direction matches the R function. Equal empty vectors
are comparable and return true. No equality tolerance is introduced. Original
arithmetic is used where representable, and common scaling is a fallback for
extreme values. R uses different accumulation precision, so borderline floating
comparisons are not guaranteed identical.
