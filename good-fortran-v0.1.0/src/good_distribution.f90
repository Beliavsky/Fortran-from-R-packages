! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

module good_distribution
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use good_kinds, only : dp
    use good_special, only : good_series_stats, quiet_nan
    implicit none
    private

    public :: dgood, pgood, qgood, rgood
    public :: goodmean, good_moments
    public :: good_logpmf

contains

    function good_logpmf(x, z, s, status) result(logp)
        integer, intent(in) :: x
        real(dp), intent(in) :: z, s
        integer, intent(out), optional :: status
        real(dp) :: logp

        real(dp) :: log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        integer :: istat

        if (x < 0 .or. z <= 0.0_dp .or. z >= 1.0_dp) then
            logp = -huge(1.0_dp)
            if (present(status)) status = -1
            return
        end if

        call good_series_stats(z, s, log_norm, mean_n, var_n, mean_log_n, &
                               cov_n_log_n, istat)
        logp = real(x + 1, dp) * log(z) - s * log(real(x + 1, dp)) - log_norm
        if (present(status)) status = istat
    end function good_logpmf

    function dgood(x, z, s, status) result(p)
        integer, intent(in) :: x
        real(dp), intent(in) :: z, s
        integer, intent(out), optional :: status
        real(dp) :: p

        real(dp) :: logp
        integer :: istat

        if (z <= 0.0_dp .or. z >= 1.0_dp) then
            p = quiet_nan()
            if (present(status)) status = -1
        else if (x < 0) then
            p = 0.0_dp
            if (present(status)) status = -2
        else
            logp = good_logpmf(x, z, s, istat)
            if (logp < log(tiny(1.0_dp))) then
                p = 0.0_dp
            else
                p = exp(logp)
            end if
            if (present(status)) status = istat
        end if
    end function dgood

    function pgood(q, z, s, lower_tail, status) result(p)
        real(dp), intent(in) :: q, z, s
        logical, intent(in), optional :: lower_tail
        integer, intent(out), optional :: status
        real(dp) :: p

        real(dp) :: log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        real(dp) :: term, cdf, logz
        integer :: k, nq, istat
        logical :: lower

        lower = .true.
        if (present(lower_tail)) lower = lower_tail

        if (z <= 0.0_dp .or. z >= 1.0_dp .or. q < 0.0_dp) then
            p = quiet_nan()
            if (present(status)) status = -1
            return
        end if

        nq = floor(q)
        call good_series_stats(z, s, log_norm, mean_n, var_n, mean_log_n, &
                               cov_n_log_n, istat)
        logz = log(z)
        cdf = 0.0_dp
        do k = 0, nq
            term = exp(real(k + 1, dp) * logz - s * log(real(k + 1, dp)) - log_norm)
            cdf = cdf + term
        end do
        cdf = min(1.0_dp, max(0.0_dp, cdf))
        if (lower) then
            p = cdf
        else
            p = max(0.0_dp, 1.0_dp - cdf)
        end if
        if (present(status)) status = istat
    end function pgood

    function qgood(p, z, s, lower_tail, status, max_q) result(q)
        real(dp), intent(in) :: p, z, s
        logical, intent(in), optional :: lower_tail
        integer, intent(out), optional :: status
        integer, intent(in), optional :: max_q
        real(dp) :: q

        real(dp) :: target, log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        real(dp) :: cdf, pmf, logz
        integer :: x, istat, qmax
        logical :: lower

        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        qmax = 10000000
        if (present(max_q)) qmax = max(1, max_q)

        if (z <= 0.0_dp .or. z >= 1.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            q = quiet_nan()
            if (present(status)) status = -1
            return
        end if

        if (lower) then
            target = p
        else
            target = 1.0_dp - p
        end if

        if (target <= 0.0_dp) then
            q = 0.0_dp
            if (present(status)) status = 0
            return
        else if (target >= 1.0_dp) then
            q = ieee_value(0.0_dp, ieee_positive_inf)
            if (present(status)) status = 0
            return
        end if

        call good_series_stats(z, s, log_norm, mean_n, var_n, mean_log_n, &
                               cov_n_log_n, istat)
        logz = log(z)
        cdf = 0.0_dp
        do x = 0, qmax
            pmf = exp(real(x + 1, dp) * logz - s * log(real(x + 1, dp)) - log_norm)
            cdf = cdf + pmf
            if (cdf >= target) then
                q = real(x, dp)
                if (present(status)) status = istat
                return
            end if
        end do

        q = real(qmax, dp)
        if (present(status)) status = 2
    end function qgood

    subroutine rgood(n, z, s, values, th, status)
        integer, intent(in) :: n
        real(dp), intent(in) :: z, s
        integer, intent(out) :: values(:)
        real(dp), intent(in), optional :: th
        integer, intent(out), optional :: status

        integer :: i, istat, qlo, qhi
        real(dp) :: threshold, u, qv

        if (size(values) < max(0, n)) then
            if (present(status)) status = -3
            return
        end if
        if (n <= 0 .or. z <= 0.0_dp .or. z >= 1.0_dp) then
            if (size(values) > 0) values = 0
            if (present(status)) status = -1
            return
        end if

        threshold = 1.0e-6_dp
        if (present(th)) threshold = th
        if (threshold < 0.0_dp) threshold = 1.0e-6_dp
        threshold = max(1.0e-12_dp, min(0.49_dp, threshold))

        qv = qgood(threshold, z, s, status=istat)
        qlo = max(0, int(qv))
        qv = qgood(1.0_dp - threshold, z, s, status=istat)
        qhi = max(qlo, int(qv))

        do i = 1, n
            call random_number(u)
            qv = qgood(u, z, s, status=istat, max_q=max(1000, qhi * 10 + 100))
            values(i) = max(0, int(qv))
        end do
        if (present(status)) status = istat
    end subroutine rgood

    function goodmean(z, s, status) result(mu)
        real(dp), intent(in) :: z, s
        integer, intent(out), optional :: status
        real(dp) :: mu

        real(dp) :: log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        integer :: istat

        call good_series_stats(z, s, log_norm, mean_n, var_n, mean_log_n, &
                               cov_n_log_n, istat)
        if (istat < 0) then
            mu = quiet_nan()
        else
            mu = mean_n - 1.0_dp
        end if
        if (present(status)) status = istat
    end function goodmean

    pure subroutine good_moments(z, s, mean, variance, status)
        real(dp), intent(in) :: z, s
        real(dp), intent(out) :: mean, variance
        integer, intent(out), optional :: status

        real(dp) :: log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        integer :: istat

        call good_series_stats(z, s, log_norm, mean_n, var_n, mean_log_n, &
                               cov_n_log_n, istat)
        if (istat < 0) then
            mean = quiet_nan()
            variance = quiet_nan()
        else
            mean = mean_n - 1.0_dp
            variance = var_n
        end if
        if (present(status)) status = istat
    end subroutine good_moments

end module good_distribution
