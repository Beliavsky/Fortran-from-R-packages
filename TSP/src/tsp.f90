! Modern Fortran translation of computational code from TSP 1.2.7.
! Original Copyright (C) Michael Hahsler and Kurt Hornik.
! SPDX-License-Identifier: GPL-3.0-only
! See LICENSE, COPYING, and UPSTREAM.md for provenance and licensing.

module tsp
    use tsp_kinds, only : dp
    use tsp_types
    use tsp_core
    use tsp_heuristics
    use tsp_transform
    use tsp_tsplib
    implicit none
    public
end module tsp
