! SPDX-License-Identifier: GPL-2.0-or-later
module adequacy_stats
    use adequacy_kinds, only: dp
    use adequacy_math, only: mean_real, median_real, sample_variance, sort_real
    implicit none
    private

    type, public :: descriptive_result
        real(dp) :: mean = 0.0_dp
        real(dp) :: median = 0.0_dp
        real(dp), allocatable :: mode(:)
        real(dp) :: variance = 0.0_dp
        real(dp) :: skewness = 0.0_dp
        real(dp) :: kurtosis = 0.0_dp
        real(dp) :: minimum = 0.0_dp
        real(dp) :: maximum = 0.0_dp
        integer :: n = 0
    end type descriptive_result

    public :: descriptive, ttt_curve

contains

    function descriptive(x) result(out)
        real(dp), intent(in) :: x(:)
        type(descriptive_result) :: out
        real(dp) :: mu, m2
        integer :: n

        n = size(x)
        out%n = n
        if (n == 0) then
            allocate(out%mode(0))
            return
        end if

        mu = mean_real(x)
        m2 = sum((x - mu)**2) / real(n, dp)
        out%mean = mu
        out%median = median_real(x)
        out%variance = sample_variance(x)
        out%minimum = minval(x)
        out%maximum = maxval(x)
        if (m2 > 0.0_dp) then
            out%skewness = (sum((x - mu)**3) / real(n, dp)) / m2**1.5_dp
            out%kurtosis = (sum((x - mu)**4) / real(n, dp)) / (m2*m2) - 3.0_dp
        end if
        call estimate_mode(x, out%mode)
    end function descriptive

    subroutine estimate_mode(x, modes)
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: modes(:)
        logical :: all_integer
        real(dp), allocatable :: work(:), uniq(:)
        integer, allocatable :: counts(:)
        real(dp) :: tol, lo, hi, width
        integer :: i, j, n, nu, max_count, nb, bin

        n = size(x)
        tol = sqrt(epsilon(1.0_dp))
        all_integer = all(abs(x - anint(x)) < tol)

        if (all_integer) then
            work = x
            call sort_real(work)
            allocate(uniq(n), counts(n))
            nu = 0
            do i = 1, n
                if (i == 1) then
                    nu = nu + 1
                    uniq(nu) = work(i)
                    counts(nu) = 1
                else if (abs(work(i) - work(i-1)) > tol) then
                    nu = nu + 1
                    uniq(nu) = work(i)
                    counts(nu) = 1
                else
                    counts(nu) = counts(nu) + 1
                end if
            end do
            max_count = maxval(counts(1:nu))
            if (max_count == 1) then
                allocate(modes(0))
                return
            end if
            allocate(modes(count(counts(1:nu) == max_count)))
            j = 0
            do i = 1, nu
                if (counts(i) == max_count) then
                    j = j + 1
                    modes(j) = uniq(i)
                end if
            end do
            return
        end if

        lo = minval(x)
        hi = maxval(x)
        if (hi <= lo) then
            allocate(modes(1))
            modes(1) = lo
            return
        end if
        nb = max(1, ceiling(log(real(n, dp))/log(2.0_dp) + 1.0_dp))
        width = (hi - lo) / real(nb, dp)
        allocate(counts(nb))
        counts = 0
        do i = 1, n
            bin = min(nb, int((x(i) - lo) / width) + 1)
            bin = max(1, bin)
            counts(bin) = counts(bin) + 1
        end do
        max_count = maxval(counts)
        if (all(counts == counts(1))) then
            allocate(modes(0))
            return
        end if
        allocate(modes(count(counts == max_count)))
        j = 0
        do i = 1, nb
            if (counts(i) == max_count) then
                j = j + 1
                modes(j) = lo + (real(i, dp) - 0.5_dp) * width
            end if
        end do
    end subroutine estimate_mode

    subroutine ttt_curve(x, r, tnorm)
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: r(:), tnorm(:)
        real(dp), allocatable :: sx(:), trn(:)
        integer :: i, n

        n = size(x)
        allocate(r(n), tnorm(n), sx(n), trn(n))
        if (n == 0) return
        sx = x
        call sort_real(sx)
        trn(1) = real(n, dp) * sx(1)
        r(1) = 1.0_dp / real(n, dp)
        do i = 2, n
            trn(i) = trn(i-1) + real(n-i+1, dp) * (sx(i) - sx(i-1))
            r(i) = real(i, dp) / real(n, dp)
        end do
        if (abs(trn(n)) > tiny(1.0_dp)) then
            tnorm = trn / trn(n)
        else
            tnorm = 0.0_dp
        end if
    end subroutine ttt_curve

end module adequacy_stats
