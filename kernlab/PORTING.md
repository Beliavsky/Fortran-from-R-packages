# Porting notes

## Interface mapping

R S4 objects and formula/list interfaces are replaced by Fortran derived types
and explicit matrices. Dots and mixed case in R names become underscores and
lower case where needed by Fortran conventions; for example `kernelMatrix`
becomes `kernel_matrix`.

## Deliberate numerical adaptations

- Numeric kernel formulas follow the upstream R implementation.
- The string implementation currently provides the exported normalized or
  unnormalized fixed-length spectrum kernel. The upstream sequence,
  full-string, exponential, constant, and bound-range suffix-array engines are
  retained as original source but are not compiled.
- `kha` uses the deterministic converged batch kernel-PCA subspace rather than
  random online Hebbian iterations.
- `kfa` uses greedy incomplete-Cholesky basis selection and a regularized basis
  inverse.
- `ksvm` classification uses deterministic one-versus-rest SMO. The upstream R
  package uses one-versus-one native solvers and supports more SVM formulations.
- Real-target `ksvm` is a regularized kernel least-squares adaptation rather
  than epsilon-insensitive or nu-SVR.
- Integer-target `gausspr` uses the regularized multiclass kernel system also
  used by LS-SVM; Gaussian-process posterior variance is available for
  regression.
- `ipop` preserves the same objective, box bounds, and two-sided linear
  constraints but uses projected augmented-penalty optimization rather than the
  upstream LOQO interior-point implementation.
- `csi` retains incomplete-Cholesky reduction, pivoting, QR factors, and
  prediction-gain diagnostics but does not reproduce every look-ahead cache in
  the original R implementation.
- `rvm`, `kqr`, `onlearn`, and bootstrap MMD use deterministic portable
  iterations rather than R random-state or native backends.

These differences are exposed in names, result status, and documentation rather
than hidden behind claims of bit-for-bit equivalence.
