# The trustOptim package

[trustOptim](https://braunm.github.io/trustOptim/) is an R package for unconstrained optimization of nonlinear
objective functions with sparse Hessians. You should consider using
this package if:

-  you are optimizing over a large number of parameters;
-  you are able to provide a function that computes the gradient of the function analytically;
-  the Hessian of the objective function is sparse (i.e., there are
   relatively few nonzero cross-partial derivatives), and you can
   provide a function that computes the Hessian.
-  you desire the stability of a trust region optimization algorithm;
and
-  you want an optimization algorithm that uses the norm of the
   gradient being zero (and not heuristics like "is the objective function making
   progress) as a stopping rule.

A common use case for a sparse Hessian is a hierarchical model that assumes
conditional independences across heterogeneous units.  Even if the
Hessian is not sparse, and/or it is hard to compute the Hessian
analytically, this package does support BFGS and SR1 updates instead.
But the real benefit of this package comes from exploiting the
sparsity of the Hessian of the objective function.


# Getting the package

The latest release version is available on CRAN. You should be able to
install it with a simple

```
install.packages("trustOptim")`.
```

The source code is available at https://github.com/braunm/trustOptim .

# Using the package

The `trust.optim` function calls the optimizer.

The objective function is defined by three R functions, all of which
take a numeric vector as the first argument.

-  **f**:  provides the objective function as a numeric scalar value
-  **grad**:  provides the gradient as a numeric vector
-  **hess**:  provides the Hessian as a *dsCMatrix* object (from the
*Matrix* package).

Starting from parameter vector *x*, to minimize a function **f**, call

```
opt <- trust.optim(x, f, grad, hess, method="Sparse")
```

See the vignettes for examples, and the package manual for details, options, etc.
