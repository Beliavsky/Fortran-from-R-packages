# Validation

The automated test program covers:

- Rosenbrock minimization with analytic-gradient BFGS
- Rosenbrock minimization with finite-difference BFGS
- numerical Hessian recovery at the Rosenbrock optimum
- Nelder-Mead on Rosenbrock
- all three nonlinear-CG update types on an anisotropic quadratic
- active lower and upper bounds in L-BFGS-B
- deterministic simulated annealing on a multimodal function
- a custom discrete SANN swap proposal
- objective maximization using negative `fnscale`
- non-unit parameter scaling
- polymorphic user data
- invalid bounds
- standalone numerical derivatives

Build profiles used for release validation:

```console
gfortran -std=f2018 -O0 -g -fcheck=all \
  -ffpe-trap=invalid,zero,overflow ...
gfortran -std=f2018 -O3 -march=native ...
```
