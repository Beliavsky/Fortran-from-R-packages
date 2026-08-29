# Porting notes

## Numerical fidelity

The forward pass, loss functions, back-propagation gradient, and analytic
Hessian are translations of `src/nnet.c`.  Connection ordering follows the R
`norm.net` / `add.net` construction exactly for the supported exported
architecture:

- hidden unit: bias, then all inputs;
- output unit with hidden layer: bias, hidden units, then optional skip inputs;
- zero-hidden skip model: bias, then all inputs.

This means an upstream `Wts` vector for the same architecture can be supplied as
`initial_wts` without reordering.

`softmax=.true.` forces linear output units and disables the binary entropy
loss, as upstream.  `censored=.true.` forces softmax and uses the probability
sum over allowed categories in each row.

## Optimization

Upstream C calls R's `vmmin` BFGS implementation.  The Fortran port calls the
MIT-licensed `r_mod::optim_bfgs` with the exact translated objective and
analytic gradient.  Fixed parameters are removed from the optimization vector
and expanded back through the mask for each objective/gradient evaluation.
Consequently fitted objective values should agree, but iteration counts and the
particular local minimum of a non-convex hidden-layer network need not be
bit-for-bit identical to R `vmmin`.

## Multinomial models

The R formula method builds a design matrix containing any desired intercept,
then fixes the network bias for the baseline-category parameterization.  The
Fortran matrix API therefore expects callers to include an intercept column in
`x` if they want one.

For response count matrices, rows are normalized and case weights multiplied by
row totals exactly as in upstream `multinom`.  With `censored=.true.`, response
rows are treated as category-allowance indicators and are not normalized.

Offsets are retained.  Binary models append one offset column with a fixed
coefficient of one.  Multiclass models append one offset column per class and
fix the corresponding class-specific coefficient to one, matching the R
`Wts`/`mask` construction.

`vcov.multinom` upstream uses a simplified `MASS::ginv`.  This port uses LAPACK
`DGESVD` and the same `sqrt(epsilon)` singular-value threshold principle.

## Helper reuse

`src/r_mod.F90` is a formatting-only copy of the previously supplied
MIT-licensed `r_mod.f90`.  It is used for BFGS optimization, uniform RNG, and QR
rank calculation.  The original byte content is retained as
`upstream/r_mod-original.f90`.

The only new general numerical helper required by this port is the small SVD
pseudoinverse wrapper around LAPACK `DGESVD`, because the supplied helper does
not expose an R-compatible `ginv` routine.

## Omitted interface machinery

Formula parsing, model frames, contrasts, S3 classes/method dispatch, printing,
interactive trace output, and formula-driven add/drop/ANOVA workflows are not
part of the Fortran API.  `summarize_rows` preserves the native computational
kernel used by `multinom(summ=2/3)`, but the R `summ` convenience modes are not
recreated as separate wrappers.
