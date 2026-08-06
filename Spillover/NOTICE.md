# Notice

This is an independent modern Fortran translation of the computational code in the R package `Spillover` version 0.1.1 by Jilber Urbina.

The upstream archive is retained unchanged at `upstream/Spillover-master.zip` for provenance. The upstream package declares the license `GPL-2`; this translated work is therefore distributed under GPL-2.0-only.

The translation does not include code from the R dependencies `vars`, `fastSOM`, `zoo`, `dplyr`, `tidyr`, or `ggplot2`. Required VAR, FEVD, permutation, rolling-window, and linear-algebra functionality was implemented directly in Fortran from the published formulas and observable package behavior.
