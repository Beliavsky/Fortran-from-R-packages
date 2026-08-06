# Translation coverage

## Upstream files

- `R/NlcOptim.R`: translated computationally in full.
- `R/solqp.R`: its role as a fallback QP solver is covered by the supplied
  Goldfarb-Idnani `quadprog` implementation plus the elastic feasibility QP.
- `man/solnl.Rd`: examples and argument semantics are reflected in the API,
  demonstration, tests, and documentation.

## Export coverage

The upstream NAMESPACE exports only `solnl`; it is available with the same
name in module `nlcoptim`.

## Omitted R infrastructure

- R matrix/list coercion and shape restoration
- `tryCatch` exception mechanics
- printed warnings and formula syntax
- MASS `ginv` calls used only by the private R QP fallback
- R random perturbations used in a singular private-fallback branch

No plotting code exists in the upstream package.
