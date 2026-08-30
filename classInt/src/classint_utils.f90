! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_utils
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use classint_kinds, only: dp
    use e1071, only: rng_state, rng_integer
    implicit none
    private

    public :: finite_values, sort_real, unique_sorted, sample_without_replacement
    public :: mean_dp, sample_sd, quantile_r, sturges_classes, lower_text

contains

    pure function finite_values(x) result(y)
        !! Returns the finite elements of `x` in their original order.
        real(dp), intent(in) :: x(:) !! Observations from which non-finite values are omitted.
        real(dp), allocatable :: y(:) !! Packed vector containing only finite observations.
        integer :: i, n

        n = count(ieee_is_finite(x))
        allocate (y(n))
        n = 0
        do i = 1, size(x)
            if (ieee_is_finite(x(i))) then
                n = n + 1
                y(n) = x(i)
            end if
        end do
    end function finite_values

    pure subroutine sort_real(x)
        !! Sorts a real vector in nondecreasing order using insertion sort.
        real(dp), intent(inout) :: x(:) !! Real vector sorted in nondecreasing order in place.
        integer :: i, j
        real(dp) :: key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_real

    pure function unique_sorted(x) result(y)
        !! Returns the distinct elements of `x` in nondecreasing order.
        real(dp), intent(in) :: x(:) !! Finite numeric values; ordering is arbitrary and duplicates are removed exactly.
        real(dp), allocatable :: y(:) !! Sorted vector containing one copy of each distinct value.
        real(dp), allocatable :: tmp(:)
        integer :: i, n

        tmp = x
        call sort_real(tmp)
        if (size(tmp) == 0) then
            allocate (y(0))
            return
        end if
        n = 1
        do i = 2, size(tmp)
            if (tmp(i) < tmp(i - 1) .or. tmp(i) > tmp(i - 1)) n = n + 1
        end do
        allocate (y(n))
        y(1) = tmp(1)
        n = 1
        do i = 2, size(tmp)
            if (tmp(i) < tmp(i - 1) .or. tmp(i) > tmp(i - 1)) then
                n = n + 1
                y(n) = tmp(i)
            end if
        end do
    end function unique_sorted

    pure function mean_dp(x) result(mu)
        !! Computes the arithmetic mean of a nonempty finite sample.
        real(dp), intent(in) :: x(:) !! Finite sample whose arithmetic mean is requested; the vector must be nonempty.
        real(dp) :: mu !! Arithmetic mean of `x`.

        if (size(x) < 1) error stop "mean_dp: empty input"
        mu = sum(x) / real(size(x), dp)
    end function mean_dp

    pure function sample_sd(x) result(sd)
        !! Computes the sample standard deviation using denominator `n - 1`.
        real(dp), intent(in) :: x(:) !! Finite sample; returns the R-style standard deviation using denominator n-1.
        real(dp) :: sd !! Sample standard deviation, or zero when fewer than two values are supplied.
        real(dp) :: mu

        if (size(x) < 2) then
            sd = 0.0_dp
            return
        end if
        mu = mean_dp(x)
        sd = sqrt(sum((x - mu)**2) / real(size(x) - 1, dp))
    end function sample_sd

    pure function quantile_r(x, p, qtype) result(q)
        !! Computes one of R's nine sample-quantile estimators.
        real(dp), intent(in) :: x(:) !! Finite observations used to compute an R-compatible sample quantile.
        real(dp), intent(in) :: p !! Quantile probability in the closed interval [0,1].
        integer, intent(in) :: qtype !! R quantile algorithm number from 1 through 9.
        real(dp) :: q !! Requested sample quantile.
        real(dp), allocatable :: s(:)
        real(dp) :: g, gamma, h, m
        integer :: j, n

        if (size(x) < 1) error stop "quantile_r: empty input"
        if (p < 0.0_dp .or. p > 1.0_dp) error stop "quantile_r: p outside [0,1]"
        if (qtype < 1 .or. qtype > 9) error stop "quantile_r: type must be 1..9"
        s = x
        call sort_real(s)
        n = size(s)
        if (p <= 0.0_dp) then
            q = s(1)
            return
        end if
        if (p >= 1.0_dp) then
            q = s(n)
            return
        end if

        select case (qtype)
        case (1)
            h = real(n, dp) * p
            j = floor(h)
            g = h - real(j, dp)
            gamma = merge(0.0_dp, 1.0_dp, abs(g) <= epsilon(1.0_dp))
        case (2)
            h = real(n, dp) * p
            j = floor(h)
            g = h - real(j, dp)
            gamma = merge(0.5_dp, 1.0_dp, abs(g) <= epsilon(1.0_dp))
        case (3)
            h = real(n, dp) * p - 0.5_dp
            j = floor(h)
            g = h - real(j, dp)
            if (abs(g) <= epsilon(1.0_dp) .and. modulo(j, 2) == 0) then
                gamma = 0.0_dp
            else
                gamma = 1.0_dp
            end if
        case default
            select case (qtype)
            case (4)
                m = 0.0_dp
            case (5)
                m = 0.5_dp
            case (6)
                m = p
            case (7)
                m = 1.0_dp - p
            case (8)
                m = (p + 1.0_dp) / 3.0_dp
            case (9)
                m = p / 4.0_dp + 3.0_dp / 8.0_dp
            end select
            h = real(n, dp) * p + m
            j = floor(h)
            gamma = h - real(j, dp)
        end select

        if (j < 1) then
            q = s(1)
        else if (j >= n) then
            q = s(n)
        else
            q = (1.0_dp - gamma) * s(j) + gamma * s(j + 1)
        end if
    end function quantile_r

    pure elemental function sturges_classes(n) result(k)
        !! Computes the number of histogram classes prescribed by Sturges' rule.
        integer, intent(in) :: n !! Number of finite observations used by Sturges' histogram class-count rule.
        integer :: k

        if (n < 1) error stop "sturges_classes: n must be positive"
        k = ceiling(log(real(n, dp)) / log(2.0_dp) + 1.0_dp)
        k = max(1, k)
    end function sturges_classes

    subroutine sample_without_replacement(x, nout, rng, y)
        !! Draws a random sample without replacement using a caller-owned RNG state.
        real(dp), intent(in) :: x(:) !! Population vector from which observations are sampled without replacement.
        integer, intent(in) :: nout !! Requested sample size between zero and the population length.
        type(rng_state), intent(inout) :: rng !! Mutable deterministic RNG state advanced by the sampling operation.
        real(dp), allocatable, intent(out) :: y(:) !! Sampled values in draw order, with length nout.
        integer, allocatable :: idx(:)
        integer :: i, j, tmp

        if (nout < 0 .or. nout > size(x)) error stop "sample_without_replacement: invalid sample size"
        allocate (idx(size(x)))
        do i = 1, size(x)
            idx(i) = i
        end do
        do i = 1, nout
            j = i - 1 + rng_integer(rng, size(x) - i + 1)
            tmp = idx(i)
            idx(i) = idx(j)
            idx(j) = tmp
        end do
        allocate (y(nout))
        do i = 1, nout
            y(i) = x(idx(i))
        end do
    end subroutine sample_without_replacement

    pure elemental function lower_text(text) result(out)
        !! Converts ASCII uppercase letters in `text` to lowercase.
        character(len=*), intent(in) :: text !! Character value normalized to lower-case ASCII for option matching.
        character(len=len(text)) :: out
        integer :: i
        integer :: code

        out = text
        do i = 1, len(text)
            code = iachar(out(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
        end do
    end function lower_text
end module classint_utils
