# Source provenance

## Supplied R package archive

- File: `lbfgs-master.zip`
- SHA-256: `32490cee4da803625a5055bd12eb7cf710e17b9a3bcc2a478bb0f789d4c66404`
- Package version in `DESCRIPTION`: 1.2.1.2
- Package license: GPL version 2 or later

The extracted source is retained in `upstream/lbfgs-master`.

## libLBFGS lineage

The supplied package embeds libLBFGS code by Naoaki Okazaki, based on Jorge
Nocedal's L-BFGS implementation. The upstream project page is:

`https://www.chokkan.org/software/liblbfgs/`

The relevant source code is under the MIT license reproduced in
`licenses/libLBFGS-MIT.txt`.

No separately downloaded Fortran implementation was incorporated. The native
Fortran code was translated from the supplied C/C++ implementation so that the
package's More-Thuente, backtracking, and OWL-QN behavior could be preserved in
one self-contained FPM project.
