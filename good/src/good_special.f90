! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

module good_special
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use good_kinds, only : dp
    implicit none
    private

    public :: good_series_stats
    public :: normal_sf_two_sided
    public :: chi_square_sf
    public :: quiet_nan

contains

    pure function quiet_nan() result(x)
        real(dp) :: x

        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function quiet_nan

    pure subroutine good_series_stats(z, s, log_norm, mean_n, var_n, mean_log_n, &
                                      cov_n_log_n, status, tol, max_terms)
        real(dp), intent(in) :: z, s
        real(dp), intent(out) :: log_norm, mean_n, var_n, mean_log_n, cov_n_log_n
        integer, intent(out) :: status
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: max_terms

        integer :: k, nmax, stable_count
        real(dp) :: tolerance, logz, logterm, scale, w, rescale
        real(dp) :: sumw, sumk, sumk2, sumlogk, sumklogk
        real(dp) :: rk, rbound, tail_rel, logk, tiny_w
        logical :: converged

        tolerance = 1.0e-13_dp
        if (present(tol)) tolerance = max(tol, epsilon(1.0_dp))
        nmax = 200000
        if (present(max_terms)) nmax = max(10, max_terms)

        if (z <= 0.0_dp .or. z >= 1.0_dp) then
            log_norm = quiet_nan()
            mean_n = quiet_nan()
            var_n = quiet_nan()
            mean_log_n = quiet_nan()
            cov_n_log_n = quiet_nan()
            status = -1
            return
        end if

        logz = log(z)
        scale = -huge(1.0_dp)
        sumw = 0.0_dp
        sumk = 0.0_dp
        sumk2 = 0.0_dp
        sumlogk = 0.0_dp
        sumklogk = 0.0_dp
        stable_count = 0
        converged = .false.

        do k = 1, nmax
            logk = log(real(k, dp))
            logterm = real(k, dp) * logz - s * logk

            if (logterm > scale) then
                if (sumw > 0.0_dp) then
                    rescale = exp(scale - logterm)
                    sumw = sumw * rescale
                    sumk = sumk * rescale
                    sumk2 = sumk2 * rescale
                    sumlogk = sumlogk * rescale
                    sumklogk = sumklogk * rescale
                end if
                scale = logterm
                w = 1.0_dp
            else
                w = exp(logterm - scale)
            end if

            sumw = sumw + w
            sumk = sumk + real(k, dp) * w
            sumk2 = sumk2 + real(k, dp) * real(k, dp) * w
            sumlogk = sumlogk + logk * w
            sumklogk = sumklogk + real(k, dp) * logk * w

            rk = z * ((real(k + 1, dp) / real(k, dp)) ** (-s))
            if (rk < 1.0_dp) then
                if (s >= 0.0_dp) then
                    rbound = z
                else
                    rbound = rk
                end if
                tiny_w = exp(logterm - scale) / sumw
                if (rbound < 1.0_dp) then
                    tail_rel = tiny_w * rk / max(1.0_dp - rbound, tiny(1.0_dp))
                else
                    tail_rel = huge(1.0_dp)
                end if
                if (tail_rel < tolerance) then
                    stable_count = stable_count + 1
                else
                    stable_count = 0
                end if
                if (stable_count >= 3) then
                    converged = .true.
                    exit
                end if
            end if
        end do

        log_norm = scale + log(sumw)
        mean_n = sumk / sumw
        var_n = max(0.0_dp, sumk2 / sumw - mean_n * mean_n)
        mean_log_n = sumlogk / sumw
        cov_n_log_n = sumklogk / sumw - mean_n * mean_log_n

        if (converged) then
            status = 0
        else
            status = 1
        end if
    end subroutine good_series_stats

    pure function normal_sf_two_sided(z) result(p)
        real(dp), intent(in) :: z
        real(dp) :: p

        p = erfc(abs(z) / sqrt(2.0_dp))
    end function normal_sf_two_sided

    pure function chi_square_sf(x, df) result(p)
        real(dp), intent(in) :: x
        integer, intent(in) :: df
        real(dp) :: p

        if (x <= 0.0_dp) then
            p = 1.0_dp
        else if (df <= 0) then
            p = quiet_nan()
        else
            p = gamma_q(0.5_dp * real(df, dp), 0.5_dp * x)
        end if
    end function chi_square_sf

    pure function gamma_q(a, x) result(q)
        real(dp), intent(in) :: a, x
        real(dp) :: q

        if (a <= 0.0_dp .or. x < 0.0_dp) then
            q = quiet_nan()
        else if (x <= 0.0_dp) then
            q = 1.0_dp
        else if (x < a + 1.0_dp) then
            q = 1.0_dp - gamma_p_series(a, x)
        else
            q = gamma_q_cf(a, x)
        end if
        q = min(1.0_dp, max(0.0_dp, q))
    end function gamma_q

    pure function gamma_p_series(a, x) result(p)
        real(dp), intent(in) :: a, x
        real(dp) :: p
        integer, parameter :: itmax = 10000
        real(dp), parameter :: eps = 8.0_dp * epsilon(1.0_dp)
        integer :: n
        real(dp) :: ap, del, sum

        ap = a
        del = 1.0_dp / a
        sum = del
        do n = 1, itmax
            ap = ap + 1.0_dp
            del = del * x / ap
            sum = sum + del
            if (abs(del) <= abs(sum) * eps) exit
        end do
        p = sum * exp(-x + a * log(x) - log_gamma(a))
    end function gamma_p_series

    pure function gamma_q_cf(a, x) result(q)
        real(dp), intent(in) :: a, x
        real(dp) :: q
        integer, parameter :: itmax = 10000
        real(dp), parameter :: eps = 8.0_dp * epsilon(1.0_dp)
        real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
        integer :: i
        real(dp) :: an, b, c, d, del, h

        b = x + 1.0_dp - a
        c = 1.0_dp / fpmin
        d = 1.0_dp / max(abs(b), fpmin)
        if (b < 0.0_dp) d = -d
        h = d
        do i = 1, itmax
            an = -real(i, dp) * (real(i, dp) - a)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < fpmin) d = sign(fpmin, d)
            c = b + an / c
            if (abs(c) < fpmin) c = sign(fpmin, c)
            d = 1.0_dp / d
            del = d * c
            h = h * del
            if (abs(del - 1.0_dp) <= eps) exit
        end do
        q = exp(-x + a * log(x) - log_gamma(a)) * h
    end function gamma_q_cf

end module good_special
