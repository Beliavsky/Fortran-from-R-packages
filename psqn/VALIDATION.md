# Validation

The translation environment did not contain the `fpm` executable, so validation
was performed by compiling the same FPM source tree directly with GNU Fortran
14.2.0.

Compiler checks used:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -fcheck=all -O0
```

All eight included tests passed:

```text
test_auglag: PASS
test_bfgs: PASS
test_generic: PASS
test_hessian: PASS
test_interpolation: PASS
test_private: PASS
test_richardson: PASS
test_structured: PASS
```

Both examples also compiled and ran, producing the expected exact quadratic
minimizers and status code 0.

## 0.1.1 portability fix

`optimize_core` no longer declares `constraint_evaluator` as an optional
procedure dummy. Some GNU Fortran versions on Windows lose the explicit
interface of an optional procedure dummy when it is referenced by an internal
procedure through host association, producing `-Werror=implicit-interface` at
the augmented-Lagrangian constraint callback. The callback is now always a
non-optional `procedure(psqn_element_eval)`; unconstrained calls pass a private
no-op evaluator. This preserves behavior and makes the nested callback interface
explicit on those compilers.
