# Translation notes

## Upstream

- R package: `BivGeo`
- Version: 2.1.1
- Authors: Ricardo Puziol de Oliveira and Jorge Alberto Achcar
- License declaration: GPL (>= 2)

The original package is pure R. Original source is retained in `upstream/` for provenance.

## API mapping

| R routine | Fortran routine |
|---|---|
| `dbivgeo1` | `dbivgeo1` |
| `dbivgeo2` | `dbivgeo2` |
| `pbivgeo(..., lower.tail=TRUE)` | `pbivgeo` |
| `pbivgeo(..., lower.tail=FALSE)` | `sbivgeo` or `pbivgeo(...,.false.)` |
| `covbivgeo` | `covbivgeo` |
| `corbivgeo` | `corbivgeo` |
| `cfbivgeo` | `cfbivgeo` |
| `mombivgeo` | `mombivgeo` |
| `rbivgeo1` | `rbivgeo1` |
| `rbivgeo2` | `rbivgeo2` |

## Numerical/interface differences

1. The R `rbivgeo1` implementation scans conditional-support values only up to `n`, so a draw can remain zero if its conditional quantile lies beyond that arbitrary limit. The Fortran version implements the same conditional distribution analytically and samples its three regions exactly, with no `n`-dependent support truncation.
2. At `theta3 = 1`, the shock component in `rbivgeo2` has infinite lifetime. The Fortran version handles this boundary directly by drawing independent marginal geometric variables rather than attempting a geometric draw with success probability zero.
3. Invalid parameters produce IEEE quiet NaNs in scalar numerical functions; sampling/estimation procedures provide an optional `ok` flag.
4. The Fortran API uses a `bivgeo_params` derived type instead of R numeric vectors.

## Source form

All compiled Fortran is free-format `.f90`, uses `implicit none`, and requires no C interoperability.
