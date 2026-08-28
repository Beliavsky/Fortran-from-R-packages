# Testing

The test suite covers:

- Richardson, simple, one-sided, elementwise, and complex-step gradients
- Richardson and complex-step Jacobians
- Richardson and hybrid complex-step Hessians
- `r=1` operation
- Bates-Watts `D` matrix layout and polynomial derivatives
- The upstream Puromycin/Bates-Watts regression example at its original
  maximum relative-error threshold of `1e-6`

Run with FPM:

```text
fpm test
```

Or use the strict GNU Fortran scripts in `scripts/`. The strict configuration
uses Fortran 2018 conformance, all common warnings as errors, bounds/runtime
checking, backtraces, and traps for invalid operations, division by zero, and
overflow.
