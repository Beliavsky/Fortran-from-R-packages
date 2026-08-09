# Translation coverage

Upstream: smoof 1.7.0.

## Translated

### Single-objective benchmarks

All ordinary numeric formulas from the 77 `sof.*` generators are translated,
except the MPM2 randomized peak-model generator. This includes Ackley, Adjiman,
Alpine01/02, Aluffi-Pentini, Bartels-Conn, Beale, Bent Cigar, Bird,
Bohachevsky N1, Booth, Branin, Brent, Brown, Bukin N2/N4/N6, Carrom Table,
Chichinadze, Chung-Reynolds, Complex, Cosine Mixture, Cross-in-Tray, Cube,
Deckkers-Aarts, Deflected Corrugated Spring, Dixon-Price, Double Sum,
Generalized Drop-Wave, Easom, Eggcrate, Eggholder, El-Attar-Vidyasagar-Dutta,
Engvall, Exponential, Freudenstein-Roth, Giunta, Goldstein-Price, Griewank,
Hansen, Hartmann 3/4/6, Himmelblau, Holder Table N1/N2, Hosaki,
Hyper-Ellipsoid, Inverted Vincent, Jennrich-Sampson, Judge, Keane, Kearfott,
Leon, Matyas, McCormick, Michalewicz, Modified Rastrigin, Periodic, Powell Sum,
Price N1/N2/N4, Rastrigin, Rosenbrock, Schaffer N2/N4, Schwefel, Shekel,
Shubert, Six-Hump Camel, Sphere, Styblinski-Tang, Sum of Different Powers,
Swiler2014, Three-Hump Camel, Trecanni, and Zettl.

### Multi-objective

- DTLZ1-DTLZ7
- ZDT1, ZDT2, ZDT3, ZDT4, ZDT6
- WFG1-WFG9, including the upstream transformation and shape functions
- MOP1-MOP7
- BK1, Viennet, Kursawe, Dent, bi-sphere
- CEC2009 UF1-UF10
- CEC2019 SYMPART simple/rotated, OMNI, MMF1-MMF15a
- ED1/ED2

### NK

The native NK evaluator is translated. Callers may generate/store their own
links and contribution tables and evaluate bit strings without R.

## Not translated in v0.1.0

- BBOB F1-F24 instance generation and the 92 BiObjBBOB combinations. These
  depend on a sizeable stateful C instance/rotation/seed subsystem, not just
  benchmark formulas.
- MPM2 randomized peak generation and its Python interoperability layer.
- NK/rMNK random landscape generation, correlation construction, and R import/
  export serialization. The evaluation kernel is present.
- GOMOP's wrapper that composes arbitrary other smoof function objects. In
  Fortran callers can simply call their selected functions directly.
- R objective-function classes, ParamHelpers parameter sets, tags, logging and
  counting wrappers, plotting/autoplot, Pareto-front visualization, and other
  R object/graphics infrastructure.

These omissions are not replaced by different algorithms under the same names.
The complete upstream source remains under `original/smoof-master/`.
