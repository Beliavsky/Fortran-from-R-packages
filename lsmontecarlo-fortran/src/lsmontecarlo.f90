! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmontecarlo
    use lsmc_european, only : eu_call_bs, eu_put_bs
    use lsmc_kinds, only : dp
    use lsmc_pricing, only : amer_put_lsm, amer_put_lsm_av, amer_put_lsm_cv
    use lsmc_pricing, only : american_put_lsmc, american_put_lsmc_antithetic
    use lsmc_pricing, only : american_put_lsmc_control
    use lsmc_pricing, only : asian_amer_put_lsm, asian_american_put_lsmc
    use lsmc_pricing, only : quanto_amer_put_lsm, quanto_amer_put_lsm_av
    use lsmc_pricing, only : quanto_american_put_lsmc, quanto_american_put_lsmc_antithetic
    use lsmc_random, only : seed_random_number
    use lsmc_simulation, only : fast_gbm, first_value_row
    use lsmc_surface, only : amer_put_lsm_price_surface
    use lsmc_surface, only : american_put_lsmc_price_surface
    use lsmc_surface, only : asian_amer_put_lsm_price_surface
    use lsmc_surface, only : asian_american_put_lsmc_price_surface
    use lsmc_surface, only : quanto_amer_put_lsm_price_surface
    use lsmc_surface, only : quanto_american_put_lsmc_price_surface
    use lsmc_types, only : option_result, price, price_surface
    use lsmc_types, only : surface_maximum, surface_mean, surface_minimum
    implicit none
    public
end module lsmontecarlo
