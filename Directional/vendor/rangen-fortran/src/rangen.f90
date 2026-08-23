module rangen
    use rangen_kinds, only : dp, i8
    use rangen_pcg32, only : pcg32_state
    use rangen_distributions, dist_set_seed => set_seed
    use rangen_sampling, only : set_sampling_seed, sample_int, sample_real, col_sample, row_sample, nano_time
    implicit none
    private

    public :: dp, i8, pcg32_state
    public :: set_seed, set_sampling_seed, seed_all
    public :: runif, rbeta, rexp, rchisq, rgamma, rgeom, rcauchy, rt
    public :: rpareto, rfrechet, rlaplace, rgumbel, rgumble, rarcsine, rnorm
    public :: runif_mat, rbeta_mat, rexp_mat, rchisq_mat, rgamma_mat, rgeom_mat
    public :: rcauchy_mat, rt_mat, rpareto_mat, rfrechet_mat, rlaplace_mat
    public :: rgumbel_mat, rgumble_mat, rarcsine_mat, rnorm_mat
    public :: col_runif, col_rbeta, col_rexp, col_rchisq, col_rgamma, col_rgeom
    public :: col_rcauchy, col_rt, col_rpareto, col_rfrechet, col_rlaplace
    public :: col_rgumbel, col_rgumble, col_rarcsine, col_rnorm
    public :: sample_int, sample_real, col_sample, row_sample, nano_time
    public :: euler_gamma

contains

    subroutine set_seed(seed_value)
        integer(i8), intent(in) :: seed_value
        call dist_set_seed(seed_value)
        call set_sampling_seed(seed_value)
    end subroutine set_seed

    subroutine seed_all(seed_value)
        integer(i8), intent(in) :: seed_value
        call set_seed(seed_value)
    end subroutine seed_all

end module rangen
