module tclust_mod
    use tclust_kinds, only : dp
    use tclust_types, only : tclust_result, tkmeans_result, rlg_result, simulation_result, &
                             rand_index_result, fm_index_result, discr_fact_result, ctlcurves_result, &
                             tclust_ic_result, tclust_ic_solution_result
    use tclust_core, only : tclust, tclust_refine, tclust_information_criteria
    use tclust_tkmeans, only : tkmeans, tkmeans_refine
    use tclust_rlg, only : rlg, rlg_refine
    use tclust_metrics, only : rand_index_labels, rand_index_table, &
                               fowlkes_mallows_labels, fowlkes_mallows_table
    use tclust_simulation, only : simulate_tclust, simulate_rlg
    use tclust_diagnostics, only : discr_fact, ctlcurves
    use tclust_grid, only : tclust_ic, tclust_ic_solutions
    implicit none
    private

    public :: dp
    public :: tclust_result
    public :: tkmeans_result
    public :: rlg_result
    public :: simulation_result
    public :: rand_index_result
    public :: fm_index_result
    public :: discr_fact_result
    public :: ctlcurves_result
    public :: tclust_ic_result
    public :: tclust_ic_solution_result
    public :: tclust
    public :: tclust_refine
    public :: tclust_information_criteria
    public :: tkmeans
    public :: tkmeans_refine
    public :: rlg
    public :: rlg_refine
    public :: rand_index_labels
    public :: rand_index_table
    public :: fowlkes_mallows_labels
    public :: fowlkes_mallows_table
    public :: simulate_tclust
    public :: simulate_rlg
    public :: discr_fact
    public :: ctlcurves
    public :: tclust_ic
    public :: tclust_ic_solutions

end module tclust_mod
