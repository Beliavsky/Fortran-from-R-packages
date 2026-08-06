# Porting notes

## Array model

R can dynamically create arrays of arbitrary rank and element type. Standard
Fortran function results have a fixed declared rank and type. The port therefore
uses:

- unlimited-polymorphic assumed-rank arguments for inquiry functions;
- real rank-1/rank-2 overloads for most allocation-producing operations;
- typed derived results for two- and three-dimensional mesh grids.

This keeps the public API type-safe and usable with ordinary FPM compilers.

## Names that collide with Fortran intrinsics

The R names `mod`, `size`, `reshape`, and matrix `sum` overlap Fortran intrinsic
procedures. The port exposes `matlab_mod`, `shape_of`/`size_dim`, `reshape2d`,
and `sum_cols`, respectively. Fortran's intrinsic `sum` remains available for
ordinary vector or whole-array summation.

## `repmat` source compatibility

The upstream R implementation has a special branch when the final replication
dimension equals one. That branch transposes a Kronecker product and can produce
a shape different from standard MATLAB tiling. The Fortran routine preserves
this behavior by default. Pass `source_compatible=.false.` to obtain conventional
row/column tiling.

## Standard deviation

The upstream R routine accepts only `flag=0` and errors otherwise. The Fortran
routine preserves `flag=0` sample standard deviation and additionally supports a
nonzero flag for population standard deviation.

## Padding

`padarray` supports constant, circular, replicate, and symmetric modes and the
`pre`, `post`, and `both` directions. Symmetric padding duplicates edge values,
matching the index sequence used by the R source.

## Error handling

Fortran does not provide R-style dynamic exceptions. Functions that cannot
produce a conforming result, such as `reshape2d` with a changed element count or
`pascal` with an invalid mode, return an allocated `0 x 0` matrix. Dimension
validation should be performed by callers when inputs are not trusted.

## Omitted code

The following code is intentionally not translated:

- `imagesc`, `colorbar`, `jet.colors`, and `multiline.plot.colors`;
- R S4 classes/method dispatch and namespace hooks;
- heterogeneous list-backed `cell` arrays;
- R-specific expression capture and data-class coercion.
