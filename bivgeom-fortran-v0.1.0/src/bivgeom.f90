! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
module bivgeom
    use bivgeom_kinds, only : dp
    use bivgeom_math, only : seed_rng
    use bivgeom_distribution, only : feasible_roy, dbivgeom_roy, log_dbivgeom_roy, &
        fbivgeom_roy, sbivgeom_roy, fyxbivgeom_roy, eyxbivgeom_roy, lambda1_roy, &
        lambda2_roy, corbivgeom_roy, relbivgeom_roy, rbivgeom_roy, empirical_survival_roy
    use bivgeom_estimation, only : bivgeom_fit, negative_loglik_roy, fit_bivgeom_ml, &
        estbivgeom_roy, estimate_ls_roy, estimate_mmp_roy, estimate_mm1_roy, &
        estimate_mm2_roy, estimate_mm3_roy, estimate_mm4_roy
    use bivgeom_compat, only : dbivgeomroy, fbivgeomroy, sbivgeomroy, fyxbivgeomroy, &
        eyxbivgeomroy, corbivgeomroy, relbivgeomroy, rbivgeomroy, estbivgeomroy, &
        loglikgeomroy, minuslogroy, lambda1roy, lambda2roy, s_n
    implicit none
    public

end module bivgeom
