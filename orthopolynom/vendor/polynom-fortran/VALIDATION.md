# Validation

The deterministic suite contains five programs:

1. **Arithmetic and interpolation**
   - Published arithmetic fixture from the upstream manual
   - Quotient/remainder reconstruction
   - Lagrange interpolation
   - Real and complex root recovery
2. **Calculus and origin shift**
   - Published degree-five root polynomial
   - Derivative/integral inversion
   - Definite integration
   - Published `change.origin(..., 3)` fixture
3. **Orthogonal polynomials**
   - Orthonormal cross-product identity
   - Coefficients from the upstream `x = (0,1,2,4)` example
   - Monic unnormalized recurrence
4. **GCD, LCM, and summaries**
   - Pair and polylist reductions
   - Zeros and stationary points
5. **API and edge cases**
   - Coefficient transforms
   - Zero-polynomial monic diagnostic
   - Matrix interpolation
   - Polylist reductions and calculus
   - String formatting

The example additionally reports root recovery and an orthonormality error near
machine precision.
