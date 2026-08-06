# Validation

The validation suite covers:

- synthetic recovery of a two-component lognormal-mixture option surface;
- forward-mean matching and option-price residuals;
- PDF normalization and ordered quantiles;
- the `100 - futures price` STIR-rate transformation;
- coupon-bond dirty pricing and `tvm::xirr` yield recovery;
- bond-yield density Jacobian transformation;
- CTD probability non-negativity and unit total mass;
- Gaussian-copula spread simulation and KDE normalization; and
- end-to-end calls to every translated R computational export.

Both a runtime-checked debug build and an optimized GNU Fortran build are run by
`scripts/validate.sh`.
