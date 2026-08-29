# Notices and provenance

This project is a modern Fortran translation of the computational portions of
R package `urca` 1.3-4.

Upstream package authors/contributors include Bernhard Pfaff, Eric Zivot, and
Matthieu Stigler. The upstream package is distributed under GPL version 2 or,
at the user's option, any later version. The translation is distributed under
the same GPL-2.0-or-later terms.

James G. MacKinnon is the author/copyright holder of response-surface material
used by the upstream package. The original `urca` distribution includes
permission correspondence regarding GPL use of that code. The complete
upstream files are preserved unchanged under:

- `upstream/src/UnitRootMacKinnon.f`
- `upstream/R/MacKinnonPValues.R`
- `upstream/inst/Licenses/MacKinnonLicense.txt`

Users redistributing or modifying this project should preserve the upstream
license, authorship, notices, and appropriate academic citations.

The user-supplied `nlme-fortran` translation is retained under
`reference/nlme-fortran/`. It is not incorporated into the compiled `urca`
library and remains subject to its own retained license/notices.
