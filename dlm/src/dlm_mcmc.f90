! SPDX-License-Identifier: GPL-2.0-or-later
module dlm_mcmc
    use dlm_types, only : dp, dlm_success, dlm_invalid_argument
    implicit none
    private

    public :: erg_mean
    public :: mcmc_mean
    public :: mcmc_sd

    interface erg_mean
        module procedure erg_mean_vector
        module procedure erg_mean_matrix
    end interface erg_mean

    interface mcmc_sd
        module procedure mcmc_sd_vector
        module procedure mcmc_sd_matrix
    end interface mcmc_sd

contains

    pure subroutine sokal_sd(x, value, clamp_negative_tau)
        real(dp), intent(in) :: x(:) !! Scalar MCMC series analyzed using the Sokal autocorrelation-window rule.
        real(dp), intent(out) :: value !! Monte Carlo standard deviation of the sample mean.
        logical, intent(in) :: clamp_negative_tau !! Replace a negative integrated autocorrelation estimate by one when true.

        real(dp), allocatable :: centered(:), ac(:)
        real(dp) :: denom, sample_variance, tau
        integer :: k, lag, lmax, n

        n = size(x)
        if (n <= 1) then
            value = 0.0_dp
            return
        end if
        lmax = min(n - 1, int(floor(30.0_dp * log10(real(n, dp)))))
        allocate(centered(n), ac(max(0, lmax)))
        centered = x - sum(x) / real(n, dp)
        denom = dot_product(centered, centered)
        if (denom <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        do lag = 1, lmax
            ac(lag) = dot_product(centered(1:n - lag), centered(1 + lag:n)) / denom
        end do
        tau = 1.0_dp
        k = 0
        do while (k < lmax .and. real(k, dp) < 3.0_dp * tau)
            k = k + 1
            tau = tau + 2.0_dp * ac(k)
        end do
        if (clamp_negative_tau .and. tau < 0.0_dp) tau = 1.0_dp
        sample_variance = denom / real(n - 1, dp)
        value = sqrt(max(0.0_dp, sample_variance * tau / real(n, dp)))
    end subroutine sokal_sd

    pure subroutine mcmc_sd_vector(x, value)
        real(dp), intent(in) :: x(:) !! Univariate MCMC output series.
        real(dp), intent(out) :: value !! Sokal-window Monte Carlo standard deviation of its mean.

        call sokal_sd(x, value, .true.)
    end subroutine mcmc_sd_vector

    pure subroutine mcmc_sd_matrix(x, value)
        real(dp), intent(in) :: x(:, :) !! Multivariate MCMC output with draws in rows and variables in columns.
        real(dp), intent(out) :: value(size(x, 2)) !! Columnwise Monte Carlo standard deviations of sample means.

        integer :: j

        do j = 1, size(x, 2)
            call sokal_sd(x(:, j), value(j), .true.)
        end do
    end subroutine mcmc_sd_matrix

    pure subroutine mcmc_mean(x, mean_value, sd_value)
        real(dp), intent(in) :: x(:, :) !! MCMC draws with iterations in rows and variables in columns.
        real(dp), intent(out) :: mean_value(size(x, 2)) !! Columnwise ergodic sample means.
        real(dp), intent(out), optional :: sd_value(size(x, 2)) !! Optional Sokal-window Monte Carlo standard deviations.

        integer :: j

        if (size(x, 1) == 0) then
            mean_value = 0.0_dp
            if (present(sd_value)) sd_value = 0.0_dp
            return
        end if
        mean_value = sum(x, dim=1) / real(size(x, 1), dp)
        if (present(sd_value)) then
            do j = 1, size(x, 2)
                call sokal_sd(x(:, j), sd_value(j), .false.)
            end do
        end if
    end subroutine mcmc_mean

    pure subroutine erg_mean_vector(x, means, info, m)
        real(dp), intent(in) :: x(:) !! Univariate simulation output in chronological order.
        real(dp), allocatable, intent(out) :: means(:) !! Cumulative/ergodic means starting after the first `m` observations.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_argument` when the starting count is invalid.
        integer, intent(in), optional :: m !! Number of observations in the first reported mean; defaults to one.

        real(dp) :: running
        integer :: i, n, start_count

        n = size(x)
        start_count = 1
        if (present(m)) start_count = max(1, m)
        if (start_count > n) then
            allocate(means(0))
            info = dlm_invalid_argument
            return
        end if
        allocate(means(n - start_count + 1))
        running = sum(x(1:start_count))
        means(1) = running / real(start_count, dp)
        do i = start_count + 1, n
            running = running + x(i)
            means(i - start_count + 1) = running / real(i, dp)
        end do
        info = dlm_success
    end subroutine erg_mean_vector

    pure subroutine erg_mean_matrix(x, means, info, m)
        real(dp), intent(in) :: x(:, :) !! Multivariate simulation output with observations in rows.
        real(dp), allocatable, intent(out) :: means(:, :) !! Rowwise sequence of cumulative column means.
        integer, intent(out) :: info !! Zero on success or `dlm_invalid_argument` when `m` exceeds the number of rows.
        integer, intent(in), optional :: m !! Number of rows in the first reported mean; defaults to one.

        real(dp), allocatable :: running(:)
        integer :: i, n, start_count

        n = size(x, 1)
        start_count = 1
        if (present(m)) start_count = max(1, m)
        if (start_count > n) then
            allocate(means(0, size(x, 2)))
            info = dlm_invalid_argument
            return
        end if
        allocate(means(n - start_count + 1, size(x, 2)), running(size(x, 2)))
        running = sum(x(1:start_count, :), dim=1)
        means(1, :) = running / real(start_count, dp)
        do i = start_count + 1, n
            running = running + x(i, :)
            means(i - start_count + 1, :) = running / real(i, dp)
        end do
        info = dlm_success
    end subroutine erg_mean_matrix

end module dlm_mcmc
