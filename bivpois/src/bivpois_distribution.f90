! Computational translation of the R package bivpois.
! License: GPL-2.0-or-later.
module bivpois_distribution
    use bivpois_kinds, only : dp
    use bivpois_math, only : rpois_one
    implicit none
    private

    type, public :: bp_grid
        integer, allocatable :: x(:)
        integer, allocatable :: y(:)
        real(dp), allocatable :: probability(:, :)
    end type bp_grid

    type, public :: bp_table
        integer :: xmin = 0
        integer :: xmax = -1
        integer :: ymin = 0
        integer :: ymax = -1
        integer, allocatable :: count(:, :)
    end type bp_table

    public :: dbp_scalar, dbp, rbp, bp_probability_grid, make_bp_table

contains

    pure real(dp) function xloglambda(n, lambda) result(v)
        integer, intent(in) :: n
        real(dp), intent(in) :: lambda
        if (n == 0) then
            v = 0.0_dp
        else if (lambda > 0.0_dp) then
            v = real(n, dp) * log(lambda)
        else
            v = -huge(1.0_dp)
        end if
    end function xloglambda

    pure real(dp) function logaddexp(a, b) result(c)
        real(dp), intent(in) :: a, b
        real(dp) :: m
        if (a <= -0.5_dp * huge(1.0_dp)) then
            c = b
        else if (b <= -0.5_dp * huge(1.0_dp)) then
            c = a
        else
            m = max(a, b)
            c = m + log(exp(a - m) + exp(b - m))
        end if
    end function logaddexp

    pure real(dp) function dbp_scalar(x, y, lambda, logged) result(den)
        integer, intent(in) :: x, y
        real(dp), intent(in) :: lambda(3)
        logical, intent(in), optional :: logged
        logical :: want_log
        integer :: k, kmax
        real(dp) :: lt, ls

        want_log = .true.
        if (present(logged)) want_log = logged
        if (any(lambda < 0.0_dp) .or. x < 0 .or. y < 0) then
            if (want_log) then
                den = -huge(1.0_dp)
            else
                den = 0.0_dp
            end if
            return
        end if

        kmax = min(x, y)
        ls = -huge(1.0_dp)
        do k = 0, kmax
            lt = -(lambda(1) + lambda(2) + lambda(3)) &
                 + xloglambda(x - k, lambda(1)) - log_gamma(real(x - k + 1, dp)) &
                 + xloglambda(y - k, lambda(2)) - log_gamma(real(y - k + 1, dp)) &
                 + xloglambda(k, lambda(3)) - log_gamma(real(k + 1, dp))
            ls = logaddexp(ls, lt)
        end do
        if (want_log) then
            den = ls
        else
            den = exp(ls)
        end if
    end function dbp_scalar

    subroutine dbp(x1, x2, lambda, den, logged)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: lambda(3)
        real(dp), intent(out) :: den(:)
        logical, intent(in), optional :: logged
        integer :: i
        logical :: want_log

        if (size(x2) /= size(x1) .or. size(den) /= size(x1)) error stop "dbp: size mismatch"
        want_log = .true.
        if (present(logged)) want_log = logged
        do i = 1, size(x1)
            den(i) = dbp_scalar(x1(i), x2(i), lambda, want_log)
        end do
    end subroutine dbp

    subroutine rbp(n, lambda, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: lambda(3)
        integer, intent(out) :: x(n, 2)
        integer :: i, z3

        if (n < 0) error stop "rbp: n must be nonnegative"
        if (any(lambda < 0.0_dp)) error stop "rbp: lambda must be nonnegative"
        do i = 1, n
            z3 = rpois_one(lambda(3))
            x(i, 1) = rpois_one(lambda(1)) + z3
            x(i, 2) = rpois_one(lambda(2)) + z3
        end do
    end subroutine rbp

    function bp_probability_grid(x1, x2, lambda, padding) result(grid)
        integer, intent(in) :: x1(:), x2(:)
        real(dp), intent(in) :: lambda(3)
        integer, intent(in), optional :: padding
        type(bp_grid) :: grid
        integer :: pad, lo1, hi1, lo2, hi2, i, j

        if (size(x1) == 0 .or. size(x2) /= size(x1)) error stop "bp_probability_grid: invalid input"
        pad = 3
        if (present(padding)) pad = max(0, padding)
        lo1 = max(minval(x1) - pad, 0)
        hi1 = maxval(x1) + pad
        lo2 = max(minval(x2) - pad, 0)
        hi2 = maxval(x2) + pad
        allocate(grid%x(lo1:hi1), grid%y(lo2:hi2))
        allocate(grid%probability(lo1:hi1, lo2:hi2))
        do i = lo1, hi1
            grid%x(i) = i
        end do
        do j = lo2, hi2
            grid%y(j) = j
        end do
        do j = lo2, hi2
            do i = lo1, hi1
                grid%probability(i, j) = dbp_scalar(i, j, lambda, .false.)
            end do
        end do
    end function bp_probability_grid

    function make_bp_table(x1, x2) result(tab)
        integer, intent(in) :: x1(:), x2(:)
        type(bp_table) :: tab
        integer :: i

        if (size(x1) == 0 .or. size(x2) /= size(x1)) error stop "make_bp_table: invalid input"
        tab%xmin = minval(x1)
        tab%xmax = maxval(x1)
        tab%ymin = minval(x2)
        tab%ymax = maxval(x2)
        allocate(tab%count(tab%xmin:tab%xmax, tab%ymin:tab%ymax), source=0)
        do i = 1, size(x1)
            tab%count(x1(i), x2(i)) = tab%count(x1(i), x2(i)) + 1
        end do
    end function make_bp_table

end module bivpois_distribution
