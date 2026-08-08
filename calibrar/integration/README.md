# External optimizer integration

The R package `calibrar` delegates many optimization methods to other packages. The Fortran core keeps these dependencies conceptually separate.

For exact external-method behavior, connect a Fortran implementation through the `scalar_objective` / `gradient_callback` interfaces and map its result to `optim_result`.

Relevant R dependencies include:

- BB (`spg`)
- cmaes
- DEoptim
- dfoptim
- GenSA
- minqa (`bobyqa`)
- optimx (`Rcgmin`, `hjn`)
- lbfgsb3c
- pso
- rgenoud
- soma

This design keeps `calibrar-fortran` self-contained while allowing separately licensed optimizer translations to be added without duplicating them inside the core library.
