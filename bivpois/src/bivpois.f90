! Umbrella module for bivpois-fortran.
! License: GPL-2.0-or-later.
module bivpois
    use bivpois_kinds, only : dp
    use bivpois_math, only : seed_rng
    use bivpois_distribution, only : bp_grid, bp_table, dbp_scalar, dbp, rbp, &
                                     bp_probability_grid, make_bp_table
    use bivpois_fit, only : bp_mle2_result, bp_mle_result, bp_profile_result, &
                           profile_loglik, bp_mle2, bp_mle, lambda3_profile
    use bivpois_gof, only : bp_gof_result, bp_dispersion_statistic, bp_gof, bp_gof2
    implicit none
    public
end module bivpois
