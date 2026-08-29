# Notices and attribution

This project is a modern Fortran translation of computational code from the R
package **truncnorm 1.0-9**.

Upstream authors listed in `DESCRIPTION`:

- Olaf Mersmann
- Heike Trautmann
- Detlef Steuer
- Bjoern Bornkamp

The upstream package is licensed **GPL (>= 2)**, represented here by
SPDX identifier `GPL-2.0-or-later`. Original sources and metadata are retained
under `upstream/truncnorm-master/`.

The accept/reject random generator follows the reference cited by the package:

J. Geweke (1991), *Efficient simulation from the multivariate normal and
student-t distributions subject to linear constraints*, Computing Science and
Statistics: Proceedings of the 23rd Symposium on the Interface, pp. 571-578.

The two-sided truncated-normal variance stabilization in `src/truncnorm.c`
cites:

J. L. Foulley (2000), *A completion simulator for the two-sided truncated
normal distribution*, Genetics Selection Evolution 32(6), 631-635.

The upstream `zeroin.c` states that it was taken from the main R distribution,
with R Core Team copyright notices, and derives from NETLIB `c/brent.shar`.
That source and notice are retained verbatim in `upstream/`.

The supplied `r_mod.f90` was provided for these translations under the MIT
License. The supplied file did not contain a copyright-holder line, so none is
invented here; see `LICENSES/MIT-r_mod.txt`.
