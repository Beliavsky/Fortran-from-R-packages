# Porting notes

## Representation

R6 classes, nested lists, formulas, hooks, and S3 dispatch are replaced by Fortran derived types and explicit procedures. The result is deterministic and suitable for batch actuarial calculations.

## Mortality tables

The upstream package obtains mortality tables from the R `MortalityTables` package. The Fortran interface takes annual `qx` and optional `ix` arrays directly. Select/period/cohort table generation is outside this project.

## Contract engine

Cash-flow conventions follow the upstream tests for annuities, deferred annuities, term-fix contracts, dread disease, endowments, pure endowments, whole life, and combined endowment/dread-disease contracts. Present values implement the upstream Thiele-style backward recursions and payment-frequency corrections.

The upstream premium and reserve machinery supports arbitrary company hooks and high-dimensional cost tensors. The Fortran engine implements the standard linear premium equation with security, acquisition, collection, administration, fixed, tax, rebate, premium-refund, and frequency loadings. Company-specific hooks should be implemented by modifying `insurance_tariff` inputs or wrapping the public procedures.

## Profit participation

The reusable calculation bases, rate-on-base formulas, accumulated profit account, terminal bonus, and terminal-bonus-fund allocation are implemented. Arbitrary R callbacks and data-frame column dispatch are not reproduced.

## Omitted infrastructure

Excel/openxlsx export, R Markdown/HTML output, RStudio project creation, plotting, data-frame formatting, contract-history snapshots, and dynamic R object inheritance are retained in `original/` but are not compiled.
