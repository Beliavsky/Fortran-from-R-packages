! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
module lsmc_simulation
    use lsmc_kinds, only : dp
    use lsmc_random, only : fill_correlated_normals, fill_normal, seed_random_number
    implicit none
    private

    public :: fast_gbm
    public :: simulate_gbm_paths
    public :: simulate_antithetic_gbm_paths
    public :: simulate_correlated_gbm_paths
    public :: simulate_antithetic_correlated_gbm_paths
    public :: first_value_row

contains

    function fast_gbm(spot, sigma, n, m, rate, dividend, maturity, seed) result(paths)
        real(dp), intent(in), optional :: spot
        real(dp), intent(in), optional :: sigma
        integer, intent(in), optional :: n
        integer, intent(in), optional :: m
        real(dp), intent(in), optional :: rate
        real(dp), intent(in), optional :: dividend
        real(dp), intent(in), optional :: maturity
        integer, intent(in), optional :: seed
        real(dp), allocatable :: paths(:, :)
        real(dp) :: d
        real(dp) :: q
        real(dp) :: s
        real(dp) :: t
        real(dp) :: v
        integer :: np
        integer :: nt

        s = 1.0_dp
        v = 0.2_dp
        np = 1000
        nt = 365
        d = 0.06_dp
        q = 0.0_dp
        t = 1.0_dp
        if (present(spot)) s = spot
        if (present(sigma)) v = sigma
        if (present(n)) np = n
        if (present(m)) nt = m
        if (present(rate)) d = rate
        if (present(dividend)) q = dividend
        if (present(maturity)) t = maturity
        if (present(seed)) call seed_random_number(seed)

        allocate(paths(np, nt))
        call simulate_gbm_paths(s, v, np, nt, d, q, t, paths)
    end function fast_gbm

    subroutine simulate_gbm_paths(spot, sigma, n_paths, n_steps, rate, dividend, maturity, paths)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: sigma
        integer, intent(in) :: n_paths
        integer, intent(in) :: n_steps
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: dividend
        real(dp), intent(in) :: maturity
        real(dp), intent(out) :: paths(n_paths, n_steps)
        real(dp), allocatable :: z(:)
        real(dp) :: drift
        real(dp) :: diffusion
        real(dp) :: dt
        integer :: i
        integer :: j

        call validate_simulation_inputs(spot, sigma, n_paths, n_steps, maturity)
        allocate(z(n_steps))
        dt = maturity / real(n_steps, dp)
        drift = (rate - dividend - 0.5_dp * sigma * sigma) * dt
        diffusion = sigma * sqrt(dt)
        do i = 1, n_paths
            call fill_normal(z)
            paths(i, 1) = spot * exp(drift + diffusion * z(1))
            do j = 2, n_steps
                paths(i, j) = paths(i, j - 1) * exp(drift + diffusion * z(j))
            end do
        end do
    end subroutine simulate_gbm_paths

    subroutine simulate_antithetic_gbm_paths(spot, sigma, n_paths, n_steps, rate, dividend, maturity, paths)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: sigma
        integer, intent(in) :: n_paths
        integer, intent(in) :: n_steps
        real(dp), intent(in) :: rate
        real(dp), intent(in) :: dividend
        real(dp), intent(in) :: maturity
        real(dp), intent(out) :: paths(2 * n_paths, n_steps)
        real(dp), allocatable :: z(:)
        real(dp) :: drift
        real(dp) :: diffusion
        real(dp) :: dt
        integer :: i
        integer :: j

        call validate_simulation_inputs(spot, sigma, n_paths, n_steps, maturity)
        allocate(z(n_steps))
        dt = maturity / real(n_steps, dp)
        drift = (rate - dividend - 0.5_dp * sigma * sigma) * dt
        diffusion = sigma * sqrt(dt)
        do i = 1, n_paths
            call fill_normal(z)
            paths(i, 1) = spot * exp(drift + diffusion * z(1))
            paths(n_paths + i, 1) = spot * exp(drift - diffusion * z(1))
            do j = 2, n_steps
                paths(i, j) = paths(i, j - 1) * exp(drift + diffusion * z(j))
                paths(n_paths + i, j) = paths(n_paths + i, j - 1) * exp(drift - diffusion * z(j))
            end do
        end do
    end subroutine simulate_antithetic_gbm_paths

    subroutine simulate_correlated_gbm_paths(spot1, sigma1, rate1, dividend1, spot2, sigma2, rate2, dividend2, &
            rho, n_paths, n_steps, maturity, paths1, paths2)
        real(dp), intent(in) :: spot1
        real(dp), intent(in) :: sigma1
        real(dp), intent(in) :: rate1
        real(dp), intent(in) :: dividend1
        real(dp), intent(in) :: spot2
        real(dp), intent(in) :: sigma2
        real(dp), intent(in) :: rate2
        real(dp), intent(in) :: dividend2
        real(dp), intent(in) :: rho
        integer, intent(in) :: n_paths
        integer, intent(in) :: n_steps
        real(dp), intent(in) :: maturity
        real(dp), intent(out) :: paths1(n_paths, n_steps)
        real(dp), intent(out) :: paths2(n_paths, n_steps)
        real(dp), allocatable :: z1(:)
        real(dp), allocatable :: z2(:)
        real(dp) :: drift1
        real(dp) :: drift2
        real(dp) :: diffusion1
        real(dp) :: diffusion2
        real(dp) :: dt
        integer :: i
        integer :: j

        call validate_simulation_inputs(spot1, sigma1, n_paths, n_steps, maturity)
        call validate_simulation_inputs(spot2, sigma2, n_paths, n_steps, maturity)
        if (abs(rho) > 1.0_dp) error stop "simulate_correlated_gbm_paths: invalid rho"
        allocate(z1(n_steps), z2(n_steps))
        dt = maturity / real(n_steps, dp)
        drift1 = (rate1 - dividend1 - 0.5_dp * sigma1 * sigma1) * dt
        drift2 = (rate2 - dividend2 - 0.5_dp * sigma2 * sigma2) * dt
        diffusion1 = sigma1 * sqrt(dt)
        diffusion2 = sigma2 * sqrt(dt)
        do i = 1, n_paths
            call fill_correlated_normals(z1, z2, rho)
            paths1(i, 1) = spot1 * exp(drift1 + diffusion1 * z1(1))
            paths2(i, 1) = spot2 * exp(drift2 + diffusion2 * z2(1))
            do j = 2, n_steps
                paths1(i, j) = paths1(i, j - 1) * exp(drift1 + diffusion1 * z1(j))
                paths2(i, j) = paths2(i, j - 1) * exp(drift2 + diffusion2 * z2(j))
            end do
        end do
    end subroutine simulate_correlated_gbm_paths

    subroutine simulate_antithetic_correlated_gbm_paths(spot1, sigma1, rate1, dividend1, spot2, sigma2, rate2, &
            dividend2, rho, n_paths, n_steps, maturity, paths1, paths2)
        real(dp), intent(in) :: spot1
        real(dp), intent(in) :: sigma1
        real(dp), intent(in) :: rate1
        real(dp), intent(in) :: dividend1
        real(dp), intent(in) :: spot2
        real(dp), intent(in) :: sigma2
        real(dp), intent(in) :: rate2
        real(dp), intent(in) :: dividend2
        real(dp), intent(in) :: rho
        integer, intent(in) :: n_paths
        integer, intent(in) :: n_steps
        real(dp), intent(in) :: maturity
        real(dp), intent(out) :: paths1(2 * n_paths, n_steps)
        real(dp), intent(out) :: paths2(2 * n_paths, n_steps)
        real(dp), allocatable :: z1(:)
        real(dp), allocatable :: z2(:)
        real(dp) :: drift1
        real(dp) :: drift2
        real(dp) :: diffusion1
        real(dp) :: diffusion2
        real(dp) :: dt
        integer :: i
        integer :: j

        call validate_simulation_inputs(spot1, sigma1, n_paths, n_steps, maturity)
        call validate_simulation_inputs(spot2, sigma2, n_paths, n_steps, maturity)
        if (abs(rho) > 1.0_dp) error stop "simulate_antithetic_correlated_gbm_paths: invalid rho"
        allocate(z1(n_steps), z2(n_steps))
        dt = maturity / real(n_steps, dp)
        drift1 = (rate1 - dividend1 - 0.5_dp * sigma1 * sigma1) * dt
        drift2 = (rate2 - dividend2 - 0.5_dp * sigma2 * sigma2) * dt
        diffusion1 = sigma1 * sqrt(dt)
        diffusion2 = sigma2 * sqrt(dt)
        do i = 1, n_paths
            call fill_correlated_normals(z1, z2, rho)
            paths1(i, 1) = spot1 * exp(drift1 + diffusion1 * z1(1))
            paths1(n_paths + i, 1) = spot1 * exp(drift1 - diffusion1 * z1(1))
            paths2(i, 1) = spot2 * exp(drift2 + diffusion2 * z2(1))
            paths2(n_paths + i, 1) = spot2 * exp(drift2 - diffusion2 * z2(1))
            do j = 2, n_steps
                paths1(i, j) = paths1(i, j - 1) * exp(drift1 + diffusion1 * z1(j))
                paths1(n_paths + i, j) = paths1(n_paths + i, j - 1) * exp(drift1 - diffusion1 * z1(j))
                paths2(i, j) = paths2(i, j - 1) * exp(drift2 + diffusion2 * z2(j))
                paths2(n_paths + i, j) = paths2(n_paths + i, j - 1) * exp(drift2 - diffusion2 * z2(j))
            end do
        end do
    end subroutine simulate_antithetic_correlated_gbm_paths

    pure function first_value_row(x) result(values)
        real(dp), intent(in) :: x(:, :)
        real(dp) :: values(size(x, 1), size(x, 2))
        logical :: found
        integer :: i
        integer :: j

        values = 0.0_dp
        do i = 1, size(x, 1)
            found = .false.
            do j = 1, size(x, 2)
                if (.not. found .and. x(i, j) > 0.0_dp) then
                    values(i, j) = x(i, j)
                    found = .true.
                end if
            end do
        end do
    end function first_value_row

    pure subroutine validate_simulation_inputs(spot, sigma, n_paths, n_steps, maturity)
        real(dp), intent(in) :: spot
        real(dp), intent(in) :: sigma
        integer, intent(in) :: n_paths
        integer, intent(in) :: n_steps
        real(dp), intent(in) :: maturity

        if (spot <= 0.0_dp) error stop "GBM simulation: spot must be positive"
        if (sigma < 0.0_dp) error stop "GBM simulation: sigma must be nonnegative"
        if (n_paths <= 0) error stop "GBM simulation: n_paths must be positive"
        if (n_steps <= 0) error stop "GBM simulation: n_steps must be positive"
        if (maturity <= 0.0_dp) error stop "GBM simulation: maturity must be positive"
    end subroutine validate_simulation_inputs

end module lsmc_simulation
