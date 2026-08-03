# Notice

This project is a modern Fortran translation of the computational code in the
R package `moments` version 0.14.1 by Lukasz Komsta and Frederick Novomestky.

The original package declares `GPL (>= 2)`. The translated code is therefore
distributed under `GPL-2.0-or-later`. Original R sources, manuals, namespace,
and package metadata are retained under `original/moments/`.

The translation changes language bindings and error reporting, adds explicit
matrix overloads and status codes, and corrects the upstream cumulant
initialization defect by default. A compatibility option retains the original
cumulant behavior.
