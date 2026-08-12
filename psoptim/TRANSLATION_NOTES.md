# Translation notes

## Source

- R package: `psoptim`
- version: 1.0
- date: 2016-01-30
- author: Krzysztof Ciupke
- package license: GPL (>= 2.0)

The package contains one computational R routine, `R/psoptim.R`, and no
compiled native code.

## Translated computational path

The Fortran routine preserves the package's global-best PSO equations:

`v = w*v + c1*r1*(pbest-x) + c2*r2*(gbest-x)`

followed by position update, personal-best update, and global-best update.
Position violations revert the affected coordinate to its previous value,
matching the R source.

## Deliberately omitted

- 2-D `image`, `contour`, `points`, convergence plot, and animation.
- `readline()` interaction.
- R list/column-name construction.
- R's global RNG implementation.

A self-contained deterministic RNG is used instead, so the same integer seed
does not imply the same random stream as R.

## Literal source quirks

Three quirks are preserved behind default-on controls and can be disabled;
see the README. This makes the Fortran translation useful both for source
compatibility and for conventional PSO use.

## Objective evaluation

The R implementation recomputes known personal/global-best fitness values
several times per generation. For a deterministic pointwise objective those
calls are redundant. The Fortran implementation caches those values and
reports the actual number of pointwise callback evaluations in `result%nfe`.
