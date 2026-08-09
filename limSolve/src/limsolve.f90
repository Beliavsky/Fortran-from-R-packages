! Upstream license declaration: GPL (version unspecified)
module limsolve
    use limsolve_kinds
    use limsolve_types
    use limsolve_linalg, only: matrix_rank, null_space, pseudoinverse
    use limsolve_inverse, only: solve_generalized, nnls, ldp, ldei, lsei, resolution
    use limsolve_linear, only: solve_tridiag, solve_banded, solve_block
    use limsolve_ranges, only: linp, xranges, varranges, varsample
    use limsolve_sampling, only: xsample, seed_rng
    implicit none
    public
end module limsolve
