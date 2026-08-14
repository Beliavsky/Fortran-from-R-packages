# Notices and provenance

This project translates the computational functionality of the R package
`mco` version 1.17 into modern Fortran.

Upstream package:

- Authors: Olaf Mersmann, Heike Trautmann, Detlef Steuer, Bernd Bischl
- NSGA-II copyright attribution: Kalyanmoy Deb
- License declared by upstream: GPL-2
- Source snapshot: `upstream/mco-master/`

The upstream hypervolume implementation credits Carlos M. Fonseca,
Luis Paquete, and Manuel Lopez-Ibanez. The Fortran port uses a new recursive
orthogonal-slicing implementation rather than copying the bundled C code.

All translated source is distributed under GPL-2.0-only. The complete license
text is in `LICENSE`.
