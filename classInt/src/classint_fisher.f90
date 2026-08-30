! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_fisher
    use classint_kinds, only: dp
    use classint_utils, only: sort_real
    implicit none
    private

    public :: fisher_exact, jenks_breaks

contains

    subroutine fisher_exact(x, k, stats, breaks)
        real(dp), intent(in) :: x(:) !! Numeric sample clustered into contiguous subsequences after sorting in ascending order.
        integer, intent(in) :: k !! Requested number of Fisher classes; must be between one and the sample size.
        real(dp), allocatable, intent(out) :: stats(:, :) !! K-by-4 upstream fish rows: min, max, mean, and population SD.
        real(dp), allocatable, intent(out) :: breaks(:) !! K+1 ascending midpoint boundaries between optimal classes.
        real(dp), allocatable :: srt(:)
        real(dp), allocatable :: work(:, :)
        integer, allocatable :: first(:, :)
        integer, allocatable :: lower(:)
        integer, allocatable :: upper(:)
        real(dp) :: ss
        real(dp) :: s
        real(dp) :: variance
        real(dp) :: candidate
        real(dp) :: mu
        real(dp) :: varpop
        integer :: m
        integer :: i
        integer :: ii
        integer :: iii
        integer :: j
        integer :: c
        integer :: top
        integer :: lo
        integer :: hi
        integer :: row

        m = size(x)
        if (m < 1) error stop "fisher_exact: empty sample"
        if (k < 1 .or. k > m) error stop "fisher_exact: invalid class count"
        srt = x
        call sort_real(srt)
        allocate(work(m, k), source=huge(1.0_dp))
        allocate(first(m, k), source=1)
        do j = 1, k
            work(1, j) = 0.0_dp
            first(1, j) = 1
        end do
        do i = 1, m
            ss = 0.0_dp
            s = 0.0_dp
            variance = 0.0_dp
            do ii = 1, i
                iii = i - ii + 1
                ss = ss + srt(iii) * srt(iii)
                s = s + srt(iii)
                variance = max(0.0_dp, ss - s * s / real(ii, dp))
                if (iii > 1) then
                    do j = 2, k
                        candidate = variance + work(iii - 1, j - 1)
                        if (work(i, j) >= candidate) then
                            first(i, j) = iii
                            work(i, j) = candidate
                        end if
                    end do
                end if
            end do
            work(i, 1) = variance
            first(i, 1) = 1
        end do

        allocate(lower(k), upper(k))
        top = m
        do c = k, 1, -1
            lo = first(top, c)
            lower(c) = lo
            upper(c) = top
            top = lo - 1
        end do
        if (top /= 0) error stop "fisher_exact: partition reconstruction failed"

        allocate(stats(k, 4))
        do c = 1, k
            lo = lower(c)
            hi = upper(c)
            mu = sum(srt(lo:hi)) / real(hi - lo + 1, dp)
            varpop = sum((srt(lo:hi) - mu)**2) / real(hi - lo + 1, dp)
            row = k - c + 1
            stats(row, 1) = srt(lo)
            stats(row, 2) = srt(hi)
            stats(row, 3) = mu
            stats(row, 4) = sqrt(max(0.0_dp, varpop))
        end do

        allocate(breaks(k + 1))
        breaks(1) = stats(k, 1)
        do c = k, 2, -1
            breaks(k - c + 2) = 0.5_dp * (stats(c, 2) + stats(c - 1, 1))
        end do
        breaks(k + 1) = stats(1, 2)
    end subroutine fisher_exact

    function jenks_breaks(x, k) result(breaks)
        real(dp), intent(in) :: x(:) !! Numeric sample sorted internally for the R Jenks dynamic-programming implementation.
        integer, intent(in) :: k !! Requested number of right-closed Jenks classes; must be between one and the sample size.
        real(dp), allocatable :: breaks(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: cost(:, :)
        integer, allocatable :: lower(:, :)
        integer, allocatable :: endpoint(:)
        real(dp) :: s1
        real(dp) :: s2
        real(dp) :: variance
        real(dp) :: value
        real(dp) :: candidate
        integer :: n
        integer :: l
        integer :: m
        integer :: i3
        integer :: i4
        integer :: j
        integer :: current

        n = size(x)
        if (n < 1) error stop "jenks_breaks: empty sample"
        if (k < 1 .or. k > n) error stop "jenks_breaks: invalid class count"
        d = x
        call sort_real(d)
        allocate(lower(n, k), source=1)
        allocate(cost(n, k), source=0.0_dp)
        if (n > 1) cost(2:n, :) = huge(1.0_dp)
        variance = 0.0_dp
        do l = 2, n
            s1 = 0.0_dp
            s2 = 0.0_dp
            do m = 1, l
                i3 = l - m + 1
                value = d(i3)
                s2 = s2 + value * value
                s1 = s1 + value
                variance = max(0.0_dp, s2 - s1 * s1 / real(m, dp))
                i4 = i3 - 1
                if (i4 /= 0) then
                    do j = 2, k
                        candidate = variance + cost(i4, j - 1)
                        if (cost(l, j) >= candidate) then
                            lower(l, j) = i3
                            cost(l, j) = candidate
                        end if
                    end do
                end if
            end do
            lower(l, 1) = 1
            cost(l, 1) = variance
        end do
        allocate(endpoint(k))
        endpoint(k) = n
        current = n
        do j = k, 2, -1
            endpoint(j - 1) = lower(current, j) - 1
            current = endpoint(j - 1)
        end do
        allocate(breaks(k + 1))
        breaks(1) = d(1)
        do j = 1, k
            breaks(j + 1) = d(endpoint(j))
        end do
    end function jenks_breaks
end module classint_fisher
