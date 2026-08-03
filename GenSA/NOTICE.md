# Notices and provenance

This project is a modern Fortran translation of the computational code in
**GenSA 1.1.15**, authored by Sylvain Gubian, Yang Xiang, Brian Suomela,
Julia Hoeng, and PMP SA.

The upstream package declares the license `GPL-2`. This translation is
distributed under **GNU GPL version 2** (`GPL-2.0-only`) and retains a complete
snapshot of the supplied upstream source under `upstream/GenSA-master/`.

The generalized simulated-annealing equations, control conventions, random
number generator, and visiting distribution are derived from the upstream
implementation. The Fortran callback API, type design, tests, examples, and
local-search implementations are part of this translation.

The upstream package cites the following primary algorithm references:

- Tsallis and Stariolo (1996), *Generalized Simulated Annealing*.
- Xiang, Sun, Fan, and Gong (1997), *Generalized Simulated Annealing Algorithm
  and Its Application to the Thomson Model*.
- Xiang and Gong (2000), *Efficiency of Generalized Simulated Annealing*.
- Xiang, Gubian, Suomela, and Hoeng (2013), *Generalized Simulated Annealing
  for Efficient Global Optimization: the GenSA Package for R*.
