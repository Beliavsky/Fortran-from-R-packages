# Translation coverage

## Direct computational coverage

- `nodeSpacer.R`: complete
- `payoff.R`: complete
- `matrixFDM.R`: complete numerical equivalent using native CSR assembly
- `multiAssetOption.R`: complete computational equivalent

## Excluded

- `plotOptionValues.R`: plotting and animation
- R argument-count checks and R object/class behavior
- R `Matrix` sparse object construction

## Added for Fortran usability

- Typed configuration and result objects
- Status/error object
- Native sparse matrix type and matrix-vector product
- Sparse iterative solver
- Multilinear interpolation
- Solver and penalty iteration diagnostics
- Automated regression tests and examples
