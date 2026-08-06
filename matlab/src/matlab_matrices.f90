! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module matlab_matrices
    use matlab_kinds, only : dp
    use matlab_array, only : rot90
    implicit none
    private

    public :: hilb
    public :: vander
    public :: vander_complex
    public :: magic
    public :: pascal
    public :: rosser

contains

    function hilb(n) result(h)
        integer, intent(in) :: n
        real(dp), allocatable :: h(:, :)
        integer :: i, j

        allocate(h(max(n, 0), max(n, 0)))
        do j = 1, n
            do i = 1, n
                h(i, j) = 1.0_dp / real(i + j - 1, dp)
            end do
        end do
    end function hilb

    function vander(v) result(a)
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: a(:, :)
        integer :: i, j, n

        n = size(v)
        allocate(a(n, n))
        do i = 1, n
            do j = 1, n
                a(i, j) = v(i) ** (n - j)
            end do
        end do
    end function vander

    function vander_complex(v) result(a)
        complex(dp), intent(in) :: v(:)
        complex(dp), allocatable :: a(:, :)
        integer :: i, j, n

        n = size(v)
        allocate(a(n, n))
        do i = 1, n
            do j = 1, n
                a(i, j) = v(i) ** (n - j)
            end do
        end do
    end function vander_complex

    recursive function magic(n) result(m)
        integer, intent(in) :: n
        real(dp), allocatable :: m(:, :)
        real(dp), allocatable :: a(:, :)
        real(dp) :: tmp
        integer :: i, j, p, k, value

        if (n <= 0) then
            allocate(m(0, 0))
        else if (n == 1) then
            allocate(m(1, 1))
            m = 1.0_dp
        else if (modulo(n, 2) == 1) then
            allocate(m(n, n))
            do j = 1, n
                do i = 1, n
                    m(i, j) = real(n * modulo(i + j - (n + 3) / 2, n) + &
                                      modulo(i + 2 * j - 2, n) + 1, dp)
                end do
            end do
        else if (modulo(n, 4) == 0) then
            allocate(m(n, n))
            do i = 1, n
                do j = 1, n
                    value = (i - 1) * n + j
                    if (modulo(i, 4) / 2 == modulo(j, 4) / 2) then
                        value = n * n + 1 - value
                    end if
                    m(i, j) = real(value, dp)
                end do
            end do
        else
            p = n / 2
            a = magic(p)
            allocate(m(n, n))
            m(1:p, 1:p) = a
            m(1:p, p + 1:n) = a + real(2 * p * p, dp)
            m(p + 1:n, 1:p) = a + real(3 * p * p, dp)
            m(p + 1:n, p + 1:n) = a + real(p * p, dp)
            k = (n - 2) / 4
            do i = 1, p
                do j = 1, k
                    tmp = m(i, j)
                    m(i, j) = m(i + p, j)
                    m(i + p, j) = tmp
                end do
                do j = n - k + 2, n
                    tmp = m(i, j)
                    m(i, j) = m(i + p, j)
                    m(i + p, j) = tmp
                end do
            end do
            i = k + 1
            tmp = m(i, 1)
            m(i, 1) = m(i + p, 1)
            m(i + p, 1) = tmp
            tmp = m(i, k + 1)
            m(i, k + 1) = m(i + p, k + 1)
            m(i + p, k + 1) = tmp
        end if
    end function magic

    function pascal(n, k) result(p)
        integer, intent(in) :: n
        integer, intent(in), optional :: k
        real(dp), allocatable :: p(:, :)
        real(dp), allocatable :: l(:, :), r(:, :)
        integer :: i, j, kk

        kk = 0
        if (present(k)) kk = k
        allocate(l(max(n, 0), max(n, 0)))
        l = 0.0_dp
        do i = 1, n
            l(i, i) = merge(1.0_dp, -1.0_dp, modulo(i - 1, 2) == 0)
            l(i, 1) = 1.0_dp
        end do
        do j = 2, n - 1
            do i = j + 1, n
                l(i, j) = l(i - 1, j) - l(i - 1, j - 1)
            end do
        end do

        select case (kk)
        case (0)
            allocate(p(n, n))
            p = matmul(l, transpose(l))
        case (1)
            allocate(p(n, n))
            p = l
        case (2)
            r = rot90(l, 3)
            allocate(p(n, n))
            p = r
            if (modulo(n, 2) == 0) p = -p
        case default
            allocate(p(0, 0))
        end select
    end function pascal

    function rosser() result(a)
        integer, allocatable :: a(:, :)

        allocate(a(8, 8))
        a = reshape([ &
             611, 196,-192, 407,  -8, -52, -49,  29, &
             196, 899, 113,-192, -71, -43,  -8, -44, &
            -192, 113, 899, 196,  61,  49,   8,  52, &
             407,-192, 196, 611,   8,  44,  59, -23, &
              -8, -71,  61,   8, 411,-599, 208, 208, &
             -52, -43,  49,  44,-599, 411, 208, 208, &
             -49,  -8,   8,  59, 208, 208,  99,-911, &
              29, -44,  52, -23, 208, 208,-911,  99], [8, 8], order=[2, 1])
    end function rosser
end module matlab_matrices
