# Notices

This project is a modern Fortran translation of the computational portions of
PortfolioAnalytics 2.1.2, whose upstream source is distributed under the GNU
General Public License version 3. The original source archive is retained in
`upstream/PortfolioAnalytics-master.zip` for provenance.

The translation preserves the upstream GPL-3.0-only license. It does not claim
to be an official release of the PortfolioAnalytics maintainers.

The R package delegates several optimization paths to external packages such as
ROI, CVXR, DEoptim, GenSA, pso, mco, OSQP, GLPK, Symphony, and quadprog. This
Fortran project provides native projected-gradient, differential-evolution,
random-search, and scalar convex routines instead. Those are adapted numerical
interfaces, not copies of the external solver implementations.
