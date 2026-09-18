module otrimle_mod
    use otrimle_kinds, only : dp
    use otrimle_types, only : otrimle_fit, otrimle_grid_point, otrimle_grid_result, &
        otrimle_sim_result, otrimle_sim_summary, kernel_density_result
    use otrimle_rng, only : set_otrimle_seed
    use otrimle_core, only : init_clust, rimle, otrimle_fit_grid
    use otrimle_density, only : kmeanfun, ksdfun, kerndensmeasure, kerndensp, kerndenscluster
    use otrimle_simulation, only : generator_otrimle, otrimleg, otrimlesimg, summarize_otrimlesimgdens
    implicit none
    private

    public :: dp
    public :: otrimle_fit
    public :: otrimle_grid_point
    public :: otrimle_grid_result
    public :: otrimle_sim_result
    public :: otrimle_sim_summary
    public :: kernel_density_result
    public :: set_otrimle_seed
    public :: init_clust
    public :: rimle
    public :: otrimle_fit_grid
    public :: kmeanfun
    public :: ksdfun
    public :: kerndensmeasure
    public :: kerndensp
    public :: kerndenscluster
    public :: generator_otrimle
    public :: otrimleg
    public :: otrimlesimg
    public :: summarize_otrimlesimgdens
end module otrimle_mod
