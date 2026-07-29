! SPDX-License-Identifier: GPL-2.0-or-later
module evir_types
    use evir_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: evir_ok = 0
    integer, parameter, public :: evir_invalid_input = 1
    integer, parameter, public :: evir_no_exceedances = 2
    integer, parameter, public :: evir_optimization_failed = 3
    integer, parameter, public :: evir_singular_hessian = 4
    integer, parameter, public :: evir_domain_error = 5

    type, public :: evir_rng
        integer(kind=8) :: state = 88172645463325252_8
    end type evir_rng

    type, public :: gev_fit_result
        integer :: status = evir_ok
        logical :: converged = .false.
        logical :: gumbel = .false.
        integer :: n = 0
        integer :: n_all = 0
        integer :: block_size = 0
        real(dp) :: xi = 0.0_dp
        real(dp) :: sigma = 1.0_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: nllh = huge(1.0_dp)
        real(dp) :: se(3) = 0.0_dp
        real(dp) :: varcov(3, 3) = 0.0_dp
        real(dp), allocatable :: data(:)
    end type gev_fit_result

    type, public :: gpd_fit_result
        integer :: status = evir_ok
        logical :: converged = .false.
        integer :: n = 0
        integer :: n_exceed = 0
        character(len=8) :: method = 'ml'
        character(len=8) :: information = 'observed'
        real(dp) :: threshold = 0.0_dp
        real(dp) :: p_less_threshold = 0.0_dp
        real(dp) :: xi = 0.0_dp
        real(dp) :: beta = 1.0_dp
        real(dp) :: nllh = huge(1.0_dp)
        real(dp) :: se(2) = 0.0_dp
        real(dp) :: varcov(2, 2) = 0.0_dp
        real(dp), allocatable :: exceedances(:)
    end type gpd_fit_result

    type, public :: pot_fit_result
        integer :: status = evir_ok
        logical :: converged = .false.
        integer :: n = 0
        integer :: n_exceed = 0
        real(dp) :: period(2) = 0.0_dp
        real(dp) :: span = 0.0_dp
        real(dp) :: threshold = 0.0_dp
        real(dp) :: p_less_threshold = 0.0_dp
        real(dp) :: intensity = 0.0_dp
        real(dp) :: run = -1.0_dp
        real(dp) :: xi = 0.0_dp
        real(dp) :: sigma = 1.0_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: beta = 1.0_dp
        real(dp) :: nllh = huge(1.0_dp)
        real(dp) :: se(3) = 0.0_dp
        real(dp) :: varcov(3, 3) = 0.0_dp
        real(dp), allocatable :: exceedances(:)
        real(dp), allocatable :: times(:)
    end type pot_fit_result

    type, public :: gpdbiv_fit_result
        integer :: status = evir_ok
        logical :: converged = .false.
        logical :: global_fit = .false.
        real(dp) :: u1 = 0.0_dp
        real(dp) :: u2 = 0.0_dp
        integer :: ne1 = 0
        integer :: ne2 = 0
        real(dp) :: lambda1 = 0.0_dp
        real(dp) :: lambda2 = 0.0_dp
        real(dp) :: alpha = 0.8_dp
        real(dp) :: alpha_se = 0.0_dp
        real(dp) :: par1(2) = 0.0_dp
        real(dp) :: par2(2) = 0.0_dp
        real(dp) :: se1(2) = 0.0_dp
        real(dp) :: se2(2) = 0.0_dp
        real(dp) :: nllh = huge(1.0_dp)
        real(dp), allocatable :: data1(:)
        real(dp), allocatable :: data2(:)
        logical, allocatable :: joint1(:)
        logical, allocatable :: joint2(:)
    end type gpdbiv_fit_result

    type, public :: xy_result
        integer :: status = evir_ok
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
    end type xy_result

    type, public :: band_result
        integer :: status = evir_ok
        real(dp), allocatable :: index(:)
        real(dp), allocatable :: estimate(:)
        real(dp), allocatable :: lower(:)
        real(dp), allocatable :: upper(:)
        real(dp), allocatable :: threshold(:)
    end type band_result

    type, public :: records_result
        integer :: status = evir_ok
        integer, allocatable :: number(:)
        integer, allocatable :: trial(:)
        real(dp), allocatable :: record(:)
        real(dp), allocatable :: expected(:)
        real(dp), allocatable :: se(:)
    end type records_result

    type, public :: decluster_result
        integer :: status = evir_ok
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: times(:)
    end type decluster_result


    type, public :: matrix_result
        integer :: status = evir_ok
        real(dp), allocatable :: values(:, :)
    end type matrix_result

    type, public :: tail_curve_result
        integer :: status = evir_ok
        real(dp) :: threshold = 0.0_dp
        real(dp) :: location = 0.0_dp
        real(dp) :: shape = 0.0_dp
        real(dp) :: scale = 1.0_dp
        real(dp), allocatable :: empirical_x(:)
        real(dp), allocatable :: empirical_y(:)
        real(dp), allocatable :: model_x(:)
        real(dp), allocatable :: model_y(:)
    end type tail_curve_result

    type, public :: profile_result
        integer :: status = evir_ok
        real(dp) :: lower = 0.0_dp
        real(dp) :: estimate = 0.0_dp
        real(dp) :: upper = 0.0_dp
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: loglik(:)
    end type profile_result

end module evir_types
