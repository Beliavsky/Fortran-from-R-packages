! SPDX-License-Identifier: LGPL-2.0-only
module lpsolve
    use lpsolve_types
    use lpsolve_core, only: solve_lp, solve_lp_sparse
    use lpsolve_special, only: lp_assign, lp_transport, make_q8, q8_triplets
    implicit none
    public
end module lpsolve
