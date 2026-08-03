# R-specific material not translated

The numerical algorithms are translated. The following R integration layers
are not separate Fortran procedures:

- S3 generics and formula methods;
- `ics` object methods and calls to `ICS::ics.components`;
- `htest` class construction and R printing conventions;
- `na.action` dispatch and data-frame coercion;
- package datasets, `.rda` serialization, and R documentation infrastructure.

The original package contains no plotting implementation to port.
