# Notices

## Upstream pracma

This project is a modern Fortran translation of `pracma` 2.4.6 by Hans W.
Borchers.  The upstream package is licensed under GPL version 3 or later.  Its
original source distribution is retained under `original/`.

## Quadratic programming

The Goldfarb-Idnani quadratic-programming implementation in
`src/quadprog_*.f90` was adapted from the GPL-2.0-or-later `quadprog` package
and modernized for this project.  GPL-2.0-or-later code may be distributed
under GPL-3.0-or-later; the combined project is therefore licensed as
GPL-3.0-or-later.

## Translation scope

R plotting, graphics-device, console, workspace, environment, and regular-
expression facilities are not part of the compiled Fortran library.  Original
files are retained for reference but do not imply that every R helper has a
Fortran counterpart.
