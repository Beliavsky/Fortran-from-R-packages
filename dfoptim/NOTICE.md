# Notices and attribution

This package is a modern Fortran translation of the computational code in the
R package **dfoptim**, version 2023.1.0.

Original package authors:

- Ravi Varadhan
- Hans W. Borchers
- Vincent Bechard

The original package is licensed under GPL version 2 or any later version.
The Fortran translation is distributed under the same terms. See `COPYING`.

The Hooke-Jeeves and modified Nelder-Mead implementations in the original
package are based on algorithms and MATLAB code described by C. T. Kelley in
*Iterative Methods for Optimization* (SIAM, 1999). The original package
documentation states that this use was made with permission of Prof. Kelley
and, for the Nelder-Mead implementation, SIAM.

The MADS implementation is based on the lower-triangular polling method
described by C. Audet and J. E. Dennis, Jr., SIAM Journal on Optimization,
17(1), 188-217 (2006).

The complete supplied R package is retained under `upstream/dfoptim-master`
for provenance and license traceability.
