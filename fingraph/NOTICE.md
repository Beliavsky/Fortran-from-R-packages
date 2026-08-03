# Notice and attribution

This project is a modern Fortran translation of the computational code in the R
package `fingraph` by Ze Vinicius and Daniel P. Palomar.

The graph operators and dense numerical support were adapted from the completed
modern Fortran translation of `spectralGraphTopology`, whose original authors
include Ze Vinicius and Daniel P. Palomar.

The deterministic random-number and multivariate Student-t simulation support
was adapted from the completed modern Fortran translation of `fitHeavyTail`, by
Daniel P. Palomar, Rui Zhou, Xiwen Wang, Frederic Pascal, Esa Ollila, and other
contributors identified in that package.

## License note

`fingraph/DESCRIPTION` declares `License: GPL-3`, while the original repository's
`LICENSE` file contains the MIT License. Both facts are preserved:

- the original MIT text is in `original/LICENSE-MIT`
- the original GPL-3 declaration is in `original/DESCRIPTION`

Because this combined translation incorporates GPL-3.0-only code from the two
prior translations, the Fortran project as distributed is GPL-3.0-only. Every
translated source, test, demo, and example contains an SPDX header.
