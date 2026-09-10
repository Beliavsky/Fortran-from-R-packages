! SPDX-License-Identifier: GPL-2.0-only
!
! Savitzky-Golay routines translated from CRAN signal.
module signal_sgolay
    use signal_kinds, only : dp
    use signal_types, only : sgolay_filter
    use signal_utils, only : factorial_real, solve_linear_system
    implicit none
    private

    public :: savitzky_golay
    public :: savitzky_golay_filter

contains

    pure function savitzky_golay(p, n, m, ts) result(filt)
        integer, intent(in) :: p !! Polynomial order, strictly smaller than n.
        integer, intent(in) :: n !! Odd filter length.
        integer, intent(in), optional :: m !! Derivative order; defaults to zero.
        real(dp), intent(in), optional :: ts !! Sample spacing used to scale derivatives; defaults to one.
        type(sgolay_filter) :: filt
        real(dp), allocatable :: c(:, :)
        real(dp), allocatable :: gram(:, :)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: lambda(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: spacing
        logical :: ok
        integer :: derivative
        integer :: k
        integer :: row
        integer :: i
        integer :: j

        derivative = 0
        if (present(m)) derivative = m
        spacing = 1.0_dp
        if (present(ts)) spacing = ts
        filt%p = p
        filt%n = n
        filt%m = derivative
        filt%ts = spacing
        if (n <= 0 .or. modulo(n, 2) == 0 .or. p >= n .or. derivative < 0 .or. derivative > p) then
            allocate(filt%coefficients(0, 0))
            return
        end if
        allocate(filt%coefficients(n, n))
        filt%coefficients = 0.0_dp
        k = n / 2
        do row = 1, k + 1
            allocate(c(n, p + 1), gram(p + 1, p + 1), rhs(p + 1))
            do i = 1, n
                do j = 1, p + 1
                    c(i, j) = real(i - row, dp) ** (j - 1)
                end do
            end do
            gram = matmul(transpose(c), c)
            rhs = 0.0_dp
            rhs(derivative + 1) = 1.0_dp
            call solve_linear_system(gram, rhs, lambda, ok)
            if (ok) then
                weights = matmul(c, lambda)
                filt%coefficients(row, :) = weights
            end if
            deallocate(c, gram, rhs)
            if (allocated(lambda)) deallocate(lambda)
            if (allocated(weights)) deallocate(weights)
        end do
        do row = k + 2, n
            filt%coefficients(row, :) = (-1.0_dp) ** derivative * &
                filt%coefficients(n - row + 1, n:1:-1)
        end do
        if (derivative > 0) then
            filt%coefficients = filt%coefficients * factorial_real(derivative) / spacing ** derivative
        end if
    end function savitzky_golay

    pure function savitzky_golay_filter(x, p, n, m, ts, supplied_filter) result(y)
        real(dp), intent(in) :: x(:) !! Input sequence to smooth or differentiate.
        integer, intent(in), optional :: p !! Polynomial order when supplied_filter is absent; defaults to three.
        integer, intent(in), optional :: n !! Odd filter length; defaults to p+3-p mod 2.
        integer, intent(in), optional :: m !! Derivative order; defaults to zero.
        real(dp), intent(in), optional :: ts !! Sample spacing; defaults to one.
        type(sgolay_filter), intent(in), optional :: supplied_filter !! Precomputed Savitzky-Golay coefficient matrix.
        real(dp), allocatable :: y(:)
        type(sgolay_filter) :: filt
        integer :: poly_order
        integer :: filter_length
        integer :: derivative
        integer :: k
        integer :: i
        real(dp) :: spacing

        if (present(supplied_filter)) then
            filt = supplied_filter
        else
            poly_order = 3
            if (present(p)) poly_order = p
            filter_length = poly_order + 3 - modulo(poly_order, 2)
            if (present(n)) filter_length = n
            derivative = 0
            if (present(m)) derivative = m
            spacing = 1.0_dp
            if (present(ts)) spacing = ts
            filt = savitzky_golay(poly_order, filter_length, derivative, spacing)
        end if
        filter_length = filt%n
        allocate(y(size(x)))
        y = 0.0_dp
        if (filter_length <= 0 .or. size(x) < filter_length) return
        k = filter_length / 2
        do i = 1, k
            y(i) = dot_product(filt%coefficients(i, :), x(1:filter_length))
        end do
        do i = k + 1, size(x) - k
            y(i) = dot_product(filt%coefficients(k + 1, :), x(i - k:i + k))
        end do
        do i = size(x) - k + 1, size(x)
            y(i) = dot_product(filt%coefficients(filter_length - (size(x) - i), :), &
                x(size(x) - filter_length + 1:size(x)))
        end do
    end function savitzky_golay_filter

end module signal_sgolay
