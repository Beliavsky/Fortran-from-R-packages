module otrimle_types
    use otrimle_kinds, only : dp
    implicit none
    private

    type, public :: otrimle_fit
        integer :: code = 0
        integer :: iter = 0
        integer :: g = 0
        logical :: flags(4) = .false.
        real(dp) :: logicd = 0.0_dp
        real(dp) :: iloglik = -huge(1.0_dp)
        real(dp) :: criterion = huge(1.0_dp)
        real(dp), allocatable :: pi(:)
        real(dp), allocatable :: mean(:, :)
        real(dp), allocatable :: cov(:, :, :)
        real(dp), allocatable :: tau(:, :)
        real(dp), allocatable :: smd(:, :)
        real(dp), allocatable :: exproportion(:)
        integer, allocatable :: cluster(:)
        integer, allocatable :: size(:)
    end type otrimle_fit

    type, public :: otrimle_grid_point
        real(dp) :: logicd = 0.0_dp
        real(dp) :: criterion = huge(1.0_dp)
        real(dp) :: iloglik = -huge(1.0_dp)
        real(dp) :: noise_proportion = 0.0_dp
        integer :: code = 0
        logical :: flags(4) = .false.
    end type otrimle_grid_point

    type, public :: otrimle_grid_result
        integer, allocatable :: g_values(:)
        type(otrimle_fit), allocatable :: solution(:)
        real(dp), allocatable :: npar(:)
        real(dp), allocatable :: ibic(:)
        real(dp), allocatable :: criterion(:)
        real(dp), allocatable :: iloglik(:)
        real(dp), allocatable :: noiseprob(:)
        real(dp), allocatable :: logicd(:)
        real(dp), allocatable :: denscrit(:)
        real(dp), allocatable :: ddpm(:, :)
    end type otrimle_grid_result

    type, public :: otrimle_sim_result
        type(otrimle_grid_result) :: result
        integer :: simruns = 0
        real(dp), allocatable :: sim_denscrit(:, :)
        logical, allocatable :: sim_ok(:, :)
    end type otrimle_sim_result

    type, public :: otrimle_sim_summary
        integer, allocatable :: g_values(:)
        integer :: best_g = 0
        real(dp) :: sdcutoff = 2.0_dp
        real(dp), allocatable :: npr(:)
        real(dp), allocatable :: nprdiff(:)
        real(dp), allocatable :: logicd(:)
        real(dp), allocatable :: denscrit(:)
        real(dp), allocatable :: mean_dens(:)
        real(dp), allocatable :: sd_dens(:)
        real(dp), allocatable :: standardized_dens(:)
        real(dp), allocatable :: penalized_g(:)
        integer, allocatable :: penorder(:)
        integer, allocatable :: cluster(:)
    end type otrimle_sim_summary

    type, public :: kernel_density_result
        real(dp) :: measure = 0.0_dp
        real(dp), allocatable :: cp(:)
        real(dp), allocatable :: cpx(:)
    end type kernel_density_result

contains
end module otrimle_types
