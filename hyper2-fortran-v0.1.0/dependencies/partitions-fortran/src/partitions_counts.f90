! Translation of computational code from R package partitions 1.10-9.
! Upstream authors: Robin K. S. Hankin; contributor Paul Egeler.
! Upstream DESCRIPTION declares: License: GPL.
! This Fortran translation is distributed under the same GPL terms.

module partitions_counts
    use partitions_kinds, only : i8
    implicit none
    private

    public :: partition_count
    public :: partition_numbers
    public :: distinct_partition_count
    public :: distinct_partition_numbers
    public :: restricted_partition_count
    public :: block_partition_count
    public :: factorial_i8
    public :: multiset_permutation_count
    public :: set_partition_count

contains

    pure function checked_add(a, b) result(c)
        integer(i8), intent(in) :: a, b
        integer(i8) :: c

        if (b > 0_i8) then
            if (a > huge(a) - b) error stop "integer overflow in partitions count"
        else if (b < 0_i8) then
            if (a < -huge(a) - b) error stop "integer underflow in partitions count"
        end if
        c = a + b
    end function checked_add

    pure function checked_mul(a, b) result(c)
        integer(i8), intent(in) :: a, b
        integer(i8) :: c

        if (a == 0_i8 .or. b == 0_i8) then
            c = 0_i8
            return
        end if
        if (a > 0_i8 .and. b > 0_i8) then
            if (a > huge(a) / b) error stop "integer overflow in partitions count"
        else if (a < 0_i8 .and. b < 0_i8) then
            if (a < huge(a) / b) error stop "integer overflow in partitions count"
        else
            if (a == -1_i8) then
                if (b == -huge(b) - 1_i8) error stop "integer overflow in partitions count"
            else if (b == -1_i8) then
                if (a == -huge(a) - 1_i8) error stop "integer overflow in partitions count"
            else if (a < 0_i8) then
                if (a < (-huge(a) - 1_i8) / b) error stop "integer overflow in partitions count"
            else
                if (b < (-huge(b) - 1_i8) / a) error stop "integer overflow in partitions count"
            end if
        end if
        c = a * b
    end function checked_mul

    function partition_numbers(n) result(p)
        integer, intent(in) :: n
        integer(i8), allocatable :: p(:)
        integer :: i, m, g1, g2, sgn
        integer(i8) :: value

        if (n < 0) error stop "partition_numbers: n must be nonnegative"
        allocate(p(0:n))
        p = 0_i8
        p(0) = 1_i8

        do i = 1, n
            value = 0_i8
            m = 1
            do
                g1 = m * (3 * m - 1) / 2
                if (g1 > i) exit
                if (mod(m, 2) == 1) then
                    sgn = 1
                else
                    sgn = -1
                end if
                value = checked_add(value, int(sgn, i8) * p(i - g1))
                g2 = m * (3 * m + 1) / 2
                if (g2 <= i) value = checked_add(value, int(sgn, i8) * p(i - g2))
                m = m + 1
            end do
            p(i) = value
        end do
    end function partition_numbers

    function partition_count(n) result(value)
        integer, intent(in) :: n
        integer(i8) :: value
        integer(i8), allocatable :: p(:)

        if (n < 0) then
            value = 0_i8
            return
        end if
        allocate(p(0:n))
        p = partition_numbers(n)
        value = p(n)
    end function partition_count

    function distinct_partition_numbers(n) result(q)
        integer, intent(in) :: n
        integer(i8), allocatable :: q(:)
        integer :: i, r, f, sgn

        if (n < 0) error stop "distinct_partition_numbers: n must be nonnegative"
        allocate(q(0:n))
        q = 0_i8
        q(0) = 1_i8
        if (n >= 1) q(1) = 1_i8

        do i = 2, n
            sgn = 1
            f = 5
            r = 2
            do while (i >= r)
                q(i) = checked_add(q(i), int(sgn, i8) * q(i - r))
                if (i == 2 * r) q(i) = checked_add(q(i), -int(sgn, i8))
                r = r + f
                f = f + 3
                sgn = -sgn
            end do

            sgn = 1
            f = 4
            r = 1
            do while (i >= r)
                q(i) = checked_add(q(i), int(sgn, i8) * q(i - r))
                if (i == 2 * r) q(i) = checked_add(q(i), -int(sgn, i8))
                r = r + f
                f = f + 3
                sgn = -sgn
            end do
        end do
    end function distinct_partition_numbers

    function distinct_partition_count(n) result(value)
        integer, intent(in) :: n
        integer(i8) :: value
        integer(i8), allocatable :: q(:)

        if (n < 0) then
            value = 0_i8
            return
        end if
        allocate(q(0:n))
        q = distinct_partition_numbers(n)
        value = q(n)
    end function distinct_partition_count

    function restricted_partition_count(m, n, include_zero) result(value)
        integer, intent(in) :: m, n
        logical, intent(in), optional :: include_zero
        integer(i8) :: value
        integer(i8), allocatable :: dp_count(:,:)
        logical :: inc0
        integer :: part, s, k, kmax

        if (m < 0 .or. n < 0) error stop "restricted_partition_count: negative argument"
        inc0 = .false.
        if (present(include_zero)) inc0 = include_zero

        if (m == 0) then
            if (n == 0) then
                value = 1_i8
            else
                value = 0_i8
            end if
            return
        end if
        if (.not. inc0 .and. m > n) then
            value = 0_i8
            return
        end if

        kmax = min(m, n)
        allocate(dp_count(0:n, 0:kmax))
        dp_count = 0_i8
        dp_count(0, 0) = 1_i8

        do part = 1, n
            do s = part, n
                do k = 1, min(kmax, s)
                    dp_count(s, k) = checked_add(dp_count(s, k), dp_count(s - part, k - 1))
                end do
            end do
        end do

        if (inc0) then
            value = 0_i8
            if (n == 0) value = 1_i8
            do k = 1, kmax
                value = checked_add(value, dp_count(n, k))
            end do
        else
            if (m <= kmax) then
                value = dp_count(n, m)
            else
                value = 0_i8
            end if
        end if
    end function restricted_partition_count

    function block_partition_count(f, n, include_fewer) result(value)
        integer, intent(in) :: f(:)
        integer, intent(in) :: n
        logical, intent(in), optional :: include_fewer
        integer(i8) :: value
        integer(i8), allocatable :: a(:), b(:)
        logical :: fewer
        integer :: i, j, k, upper

        if (any(f < 0)) error stop "block_partition_count: capacities must be nonnegative"
        if (n < 0) error stop "block_partition_count: n must be nonnegative"
        fewer = .false.
        if (present(include_fewer)) fewer = include_fewer

        allocate(a(0:n), b(0:n))
        a = 0_i8
        a(0) = 1_i8
        do i = 1, size(f)
            b = 0_i8
            do j = 0, n
                if (a(j) == 0_i8) cycle
                upper = min(f(i), n - j)
                do k = 0, upper
                    b(j + k) = checked_add(b(j + k), a(j))
                end do
            end do
            a = b
        end do

        if (fewer) then
            value = 0_i8
            do j = 0, n
                value = checked_add(value, a(j))
            end do
        else
            value = a(n)
        end if
    end function block_partition_count

    function factorial_i8(n) result(value)
        integer, intent(in) :: n
        integer(i8) :: value
        integer :: i

        if (n < 0) error stop "factorial_i8: n must be nonnegative"
        value = 1_i8
        do i = 2, n
            value = checked_mul(value, int(i, i8))
        end do
    end function factorial_i8

    function multiset_permutation_count(v) result(value)
        integer, intent(in) :: v(:)
        integer(i8) :: value
        integer, allocatable :: w(:), denoms(:)
        integer :: i, j, n, run, nd

        n = size(v)
        if (n == 0) then
            value = 1_i8
            return
        end if
        allocate(w(n), denoms(n))
        w = v
        denoms = 1
        call insertion_sort(w)
        nd = 0
        i = 1
        do while (i <= n)
            run = 1
            j = i + 1
            do while (j <= n)
                if (w(j) /= w(i)) exit
                run = run + 1
                j = j + 1
            end do
            nd = nd + 1
            denoms(nd) = run
            i = j
        end do
        value = factorial_ratio_count(n, denoms(1:nd))
    end function multiset_permutation_count

    function set_partition_count(shape) result(value)
        integer, intent(in) :: shape(:)
        integer(i8) :: value
        integer, allocatable :: s(:), denoms(:)
        integer :: n, i, j, run, nd

        if (any(shape <= 0)) error stop "set_partition_count: block sizes must be positive"
        n = sum(shape)
        allocate(denoms(2 * size(shape)))
        nd = size(shape)
        denoms(1:nd) = shape

        allocate(s(size(shape)))
        s = shape
        call insertion_sort(s)
        i = 1
        do while (i <= size(s))
            run = 1
            j = i + 1
            do while (j <= size(s))
                if (s(j) /= s(i)) exit
                run = run + 1
                j = j + 1
            end do
            nd = nd + 1
            denoms(nd) = run
            i = j
        end do
        value = factorial_ratio_count(n, denoms(1:nd))
    end function set_partition_count

    function factorial_ratio_count(n, denoms) result(value)
        integer, intent(in) :: n
        integer, intent(in) :: denoms(:)
        integer(i8) :: value
        logical, allocatable :: is_prime(:)
        integer :: p, k, exponent, d

        if (n < 0 .or. any(denoms < 0)) error stop "factorial_ratio_count: negative argument"
        if (sum(denoms) < 0) error stop "factorial_ratio_count: invalid denominator"
        if (n < 2) then
            value = 1_i8
            return
        end if
        allocate(is_prime(0:n))
        is_prime = .true.
        is_prime(0:1) = .false.
        p = 2
        do while (p * p <= n)
            if (is_prime(p)) then
                k = p * p
                do while (k <= n)
                    is_prime(k) = .false.
                    k = k + p
                end do
            end if
            p = p + 1
        end do

        value = 1_i8
        do p = 2, n
            if (.not. is_prime(p)) cycle
            exponent = factorial_valuation(n, p)
            do d = 1, size(denoms)
                exponent = exponent - factorial_valuation(denoms(d), p)
            end do
            if (exponent < 0) error stop "factorial_ratio_count: nonintegral ratio"
            do k = 1, exponent
                value = checked_mul(value, int(p, i8))
            end do
        end do
    end function factorial_ratio_count

    pure integer function factorial_valuation(n, p) result(value)
        integer, intent(in) :: n, p
        integer :: q

        value = 0
        q = n
        do while (q > 0)
            q = q / p
            value = value + q
        end do
    end function factorial_valuation

    pure subroutine insertion_sort(x)
        integer, intent(inout) :: x(:)
        integer :: i, j, key

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
    end subroutine insertion_sort

end module partitions_counts
