! SPDX-License-Identifier: GPL-3.0-or-later
module nbconv_approximations
    use nbconv_kinds, only : dp
    use nbconv_math, only : nbconv_pi, negbin_pmf_mu
    implicit none
    private

    public :: nb_sum_moments
    public :: nb_sum_saddlepoint

contains

    function nb_sum_moments(mus, phis, counts) result(pmf)
        real(dp), intent(in) :: mus(:), phis(:)
        integer, intent(in) :: counts(:)
        real(dp), allocatable :: pmf(:)
        real(dp) :: mu_moment, phi_moment, denom
        integer :: i

        call validate_mu_inputs(mus, phis)
        if (any(counts < 0)) error stop "nb_sum_moments: counts must be nonnegative"
        allocate(pmf(size(counts)))
        mu_moment = sum(mus)
        if (mu_moment <= 0.0_dp) then
            pmf = 0.0_dp
            do i = 1, size(counts)
                if (counts(i) == 0) pmf(i) = 1.0_dp
            end do
            return
        end if
        denom = sum(mus * mus / phis)
        if (denom <= 0.0_dp) then
            phi_moment = huge(1.0_dp)
        else
            phi_moment = mu_moment * mu_moment / denom
        end if
        do i = 1, size(counts)
            pmf(i) = negbin_pmf_mu(counts(i), mu_moment, phi_moment)
        end do
    end function nb_sum_moments

    function nb_sum_saddlepoint(mus, phis, counts, normalize) result(pmf)
        real(dp), intent(in) :: mus(:), phis(:)
        integer, intent(in) :: counts(:)
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: pmf(:)
        real(dp) :: logp0, total, t
        logical :: do_normalize
        integer :: i, x

        call validate_mu_inputs(mus, phis)
        if (any(counts < 0)) error stop "nb_sum_saddlepoint: counts must be nonnegative"
        allocate(pmf(size(counts)))
        if (size(counts) == 0) return

        do_normalize = .true.
        if (present(normalize)) do_normalize = normalize

        if (sum(mus) <= 0.0_dp) then
            pmf = 0.0_dp
            do i = 1, size(counts)
                if (counts(i) == 0) pmf(i) = 1.0_dp
            end do
            return
        end if

        logp0 = 0.0_dp
        do i = 1, size(mus)
            if (mus(i) > 0.0_dp) logp0 = logp0 + phis(i) * log(phis(i) / (phis(i) + mus(i)))
        end do

        do i = 1, size(counts)
            x = counts(i)
            if (x == 0) then
                pmf(i) = exp(logp0)
            else
                t = saddlepoint_root(real(x, dp), mus, phis)
                pmf(i) = exp(-0.5_dp * log(2.0_dp * nbconv_pi * k2(t, mus, phis)) &
                     + cgf(t, mus, phis) - t * real(x, dp))
            end if
        end do

        if (do_normalize) then
            total = sum(pmf)
            if (total > 0.0_dp) pmf = pmf / total
        end if
    end function nb_sum_saddlepoint

    function saddlepoint_root(x, mus, phis) result(root)
        real(dp), intent(in) :: x, mus(:), phis(:)
        real(dp) :: root
        real(dp) :: lo, hi, mid, fmid, tmax, gap, scale
        integer :: i, iter

        tmax = huge(1.0_dp)
        do i = 1, size(mus)
            if (mus(i) > 0.0_dp) tmax = min(tmax, log(1.0_dp + phis(i) / mus(i)))
        end do
        scale = max(1.0_dp, abs(tmax))
        gap = 1.0e-6_dp * scale
        hi = tmax - gap
        do while (k1(hi, mus, phis) <= x .and. gap > 1024.0_dp * epsilon(1.0_dp) * scale)
            gap = 0.1_dp * gap
            hi = tmax - gap
        end do
        lo = -100.0_dp
        do while (k1(lo, mus, phis) > x)
            lo = 2.0_dp * lo
            if (lo < -10000.0_dp) exit
        end do
        if (k1(hi, mus, phis) <= x) then
            error stop "nb_sum_saddlepoint: failed to bracket saddlepoint"
        end if

        do iter = 1, 200
            mid = 0.5_dp * (lo + hi)
            fmid = k1(mid, mus, phis) - x
            if (fmid > 0.0_dp) then
                hi = mid
            else
                lo = mid
            end if
            if (abs(hi - lo) <= 128.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        root = 0.5_dp * (lo + hi)
    end function saddlepoint_root

    pure function cgf(t, mus, phis) result(v)
        real(dp), intent(in) :: t, mus(:), phis(:)
        real(dp) :: v, et, den
        integer :: i

        et = exp(t)
        v = 0.0_dp
        do i = 1, size(mus)
            if (mus(i) <= 0.0_dp) cycle
            den = phis(i) + mus(i) * (1.0_dp - et)
            v = v + phis(i) * (log(phis(i)) - log(den))
        end do
    end function cgf

    pure function k1(t, mus, phis) result(v)
        real(dp), intent(in) :: t, mus(:), phis(:)
        real(dp) :: v, et, den
        integer :: i

        et = exp(t)
        v = 0.0_dp
        do i = 1, size(mus)
            if (mus(i) <= 0.0_dp) cycle
            den = phis(i) + mus(i) - mus(i) * et
            v = v + phis(i) * mus(i) * et / den
        end do
    end function k1

    pure function k2(t, mus, phis) result(v)
        real(dp), intent(in) :: t, mus(:), phis(:)
        real(dp) :: v, et, den
        integer :: i

        et = exp(t)
        v = 0.0_dp
        do i = 1, size(mus)
            if (mus(i) <= 0.0_dp) cycle
            den = phis(i) + mus(i) - mus(i) * et
            v = v + phis(i) * mus(i) * (phis(i) + mus(i)) * et / (den * den)
        end do
    end function k2

    subroutine validate_mu_inputs(mus, phis)
        real(dp), intent(in) :: mus(:), phis(:)

        if (size(mus) /= size(phis)) error stop "nbconv: mus and phis must have equal length"
        if (size(mus) < 1) error stop "nbconv: parameter arrays must be nonempty"
        if (any(mus < 0.0_dp)) error stop "nbconv: mus must be nonnegative"
        if (any(phis <= 0.0_dp)) error stop "nbconv: phis must be positive"
    end subroutine validate_mu_inputs

end module nbconv_approximations
