# RobStatTM lite dependency

This directory contains the minimal subset of the completed RobStatTM modern
Fortran translation required by RPEIF:

- double-precision kind and pi constant
- rho, psi, psi derivative, weight, M-scale, and tuning functions
- minimal utility interfaces needed by the psi module

The original translated `robstattm_psi.f90` and `robstattm_kinds.f90` files are
retained with GPL-3.0-or-later SPDX identifiers. The utility module is a reduced
implementation of the same interfaces.
