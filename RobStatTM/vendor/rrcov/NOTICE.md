# Notice

This work is a modern Fortran translation of computational ideas and code from
`rrcov` 1.7-8, authored primarily by Valentin Todorov and distributed under GPL
version 3 or later.

The full upstream source tree supplied for the translation is retained under
`upstream/rrcov-master`. Source files in this translation carry
`SPDX-License-Identifier: GPL-3.0-or-later`.

The translation omits graphics and R runtime infrastructure. It includes new
self-contained numerical implementations where upstream rrcov calls external R
packages such as robustbase or pcaPP.
