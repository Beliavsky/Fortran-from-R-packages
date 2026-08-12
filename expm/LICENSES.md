# Licensing and provenance

This project is a modern Fortran translation of the computational code in the
R package **expm 1.0-0**.

The R package declares `GPL (>= 2)`.  In addition, `src/matexp_MH09.c` carries
an explicit GPL version 3 or later notice.  Because this translation includes
code and algorithms derived from that file as well as the rest of the package,
**expm-fortran is distributed under GPL-3.0-or-later**.

The full GPL-3 text is in `LICENSE`; copies of GPL-2 and GPL-3 are retained in
`licenses/` for provenance.  The original `DESCRIPTION`, `NAMESPACE`, R source,
and relevant C/Fortran source files are retained under `original/` with their
original notices.

Original expm authors and contributors include Martin Maechler, Christophe
Dutang, Vincent Goulet, Douglas Bates, David Firth, Marina Shapira, Michael
Stadelmann, Roger B. Sidje, Ravi Varadhan, Drew Schmidt, and the other authors
credited in the retained source files.
