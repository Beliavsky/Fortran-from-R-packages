# Translation coverage

## Included

The complete computational content of the exported OOR API is translated:

- Parallel Optimistic Optimization (`POO`).
- Its binary HOO-style tree construction, exploration, sampling and path
  updates.
- Stochastic StoSOO.
- Deterministic SOO (`control$type = "det"`).
- Arbitrary-dimensional box domains for StoSOO/SOO.
- Minimization and maximization sign handling.
- Stochastic repeated sampling and UCB updates.
- Longest-side ternary node splitting.
- Search-tree and evaluation-history output.
- All five exported benchmark/test functions.

## Intentionally omitted

- `plotStoSOO` and all graphics calls.
- S4 `Leaf` machinery and R environments.  Equivalent state is represented by
  Fortran derived types.
- R `...` argument forwarding.  Users capture auxiliary data in module state or
  their objective procedure instead.

## Source behaviors deliberately preserved

Several details may look unusual but are present in OOR 0.1.4 and are retained
for trajectory fidelity:

1. POO's metric expression is translated as `rho ** (depth ** k)`, matching the
   source expression `rhos[k]^depth^k`, rather than replacing it with the more
   usual `rho**depth` HOO term.
2. POO loops while its new-evaluation count is `<= horizon` and then performs a
   final sample.  Consequently `result%evaluations` can exceed the nominal
   horizon.
3. StoSOO samples both side children immediately when a node is split.  The R
   budget check can therefore overshoot `nb_iter` by one evaluation; the
   Fortran translation preserves that behavior.
4. In stochastic mode the R code selects its final point by scanning for the
   deepest node with `leaf == 0`, ranking candidates by accumulated `sums`, and
   only then reporting `sums/k`.  The same rule is used here.
5. StoSOO's `sample_when_created` is hard-coded to true in the R source and is
   likewise fixed here.

## Small language-level difference

`difficult2(0.5)` leads to an undefined `log2(0) %% 1` branch in the R source.
The Fortran function returns a quiet NaN at that exact point.  Likewise,
`difficult(0)` returns a quiet NaN for its singular expression.
