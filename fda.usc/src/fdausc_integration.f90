module fdausc_integration
    use r_kinds, only : dp
    implicit none
    private

    public :: integrate_curve, integrate_curves, inprod_curves, norm_curves, cosine_proximity

contains

    pure real(dp) function integrate_curve(x, y, method) result(value)
        real(dp), intent(in) :: x(:) !! Strictly increasing grid values for the numerical integral.
        real(dp), intent(in) :: y(:) !! Function values on x; size must equal size(x).
        integer, intent(in), optional :: method !! Code: 1 trapezoid, 2 composite Simpson, 3 extended Simpson; default 2.
        integer :: n
        integer :: meth
        integer :: i
        real(dp) :: h
        real(dp), allocatable :: xx(:)
        real(dp), allocatable :: yy(:)
        logical :: equi

        n = size(x)
        if (size(y) /= n .or. n < 2) then
            value = 0.0_dp
            return
        end if
        meth = 2
        if (present(method)) meth = method
        if (n == 2) meth = 1
        equi = grid_is_equidistant(x)

        select case (meth)
        case (1)
            value = 0.0_dp
            do i = 2, n
                value = value + 0.5_dp*(x(i) - x(i - 1))*(y(i) + y(i - 1))
            end do
        case (3)
            if (.not. equi) then
                call interpolate_equal_grid(x, y, 2*n - 1, xx, yy)
                value = extended_simpson_equal(xx, yy)
            else
                value = extended_simpson_equal(x, y)
            end if
        case default
            if (.not. equi) then
                call interpolate_equal_grid(x, y, 2*n - 1, xx, yy)
                value = composite_simpson_equal(xx, yy)
            else if (mod(n, 2) == 1) then
                value = composite_simpson_equal(x, y)
            else
                ! R's implementation assumes an odd number for CSR. Preserve all
                ! observations by applying Simpson to the first n-1 points and a
                ! trapezoid to the final interval.
                value = composite_simpson_equal(x(:n - 1), y(:n - 1))
                h = x(n) - x(n - 1)
                value = value + 0.5_dp*h*(y(n) + y(n - 1))
            end if
        end select
    end function integrate_curve

    pure subroutine integrate_curves(x, curves, values, method)
        real(dp), intent(in) :: x(:) !! Common increasing grid for all rows of curves.
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row and one grid point per column.
        real(dp), intent(out) :: values(:) !! Integrals of the curve rows; size must equal size(curves,1).
        integer, intent(in), optional :: method !! Integration selector passed to integrate_curve.
        integer :: i
        integer :: meth

        meth = 2
        if (present(method)) meth = method
        do i = 1, min(size(values), size(curves, 1))
            values(i) = integrate_curve(x, curves(i, :), meth)
        end do
    end subroutine integrate_curves

    pure subroutine inprod_curves(x, curves1, curves2, products, weights, method)
        real(dp), intent(in) :: x(:) !! Common increasing argument grid for both curve matrices.
        real(dp), intent(in) :: curves1(:, :) !! First functional sample, curves in rows.
        real(dp), intent(in) :: curves2(:, :) !! Second functional sample, curves in rows.
        real(dp), intent(out) :: products(:, :) !! Pairwise integrated products with shape [size(curves1,1),size(curves2,1)].
        real(dp), intent(in), optional :: weights(:) !! Optional pointwise weights; a constant vector represents scalar weighting.
        integer, intent(in), optional :: method !! Integration selector passed to integrate_curve.
        integer :: i
        integer :: j
        integer :: meth
        real(dp), allocatable :: integrand(:)

        meth = 2
        if (present(method)) meth = method
        allocate(integrand(size(x)))
        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                integrand = curves1(i, :)*curves2(j, :)
                if (present(weights)) integrand = integrand*weights
                products(i, j) = integrate_curve(x, integrand, meth)
            end do
        end do
    end subroutine inprod_curves

    pure subroutine norm_curves(x, curves, norms, p, weights, method)
        real(dp), intent(in) :: x(:) !! Common increasing argument grid for curve rows.
        real(dp), intent(in) :: curves(:, :) !! Functional observations, one curve per row.
        real(dp), intent(out) :: norms(:) !! Lp norms, one per input curve.
        real(dp), intent(in), optional :: p !! Norm order; p=0 requests the supremum norm, default 2.
        real(dp), intent(in), optional :: weights(:) !! Optional pointwise weights applied before integration.
        integer, intent(in), optional :: method !! Integration selector passed to integrate_curve.
        integer :: i
        integer :: meth
        real(dp) :: pp
        real(dp), allocatable :: integrand(:)

        pp = 2.0_dp
        if (present(p)) pp = p
        meth = 2
        if (present(method)) meth = method
        allocate(integrand(size(x)))
        do i = 1, size(curves, 1)
            if (abs(pp) <= tiny(1.0_dp)) then
                norms(i) = maxval(abs(curves(i, :)))
            else
                integrand = abs(curves(i, :))**pp
                if (present(weights)) integrand = integrand*weights
                norms(i) = max(0.0_dp, integrate_curve(x, integrand, meth))**(1.0_dp/pp)
            end if
        end do
    end subroutine norm_curves

    pure subroutine cosine_proximity(x, curves1, curves2, rho, as_distance, method)
        real(dp), intent(in) :: x(:) !! Common increasing argument grid for both curve matrices.
        real(dp), intent(in) :: curves1(:, :) !! First functional sample, curves in rows.
        real(dp), intent(in) :: curves2(:, :) !! Second functional sample, curves in rows.
        real(dp), intent(out) :: rho(:, :) !! Pairwise cosine correlations or one-minus-absolute-correlation distances.
        logical, intent(in), optional :: as_distance !! If true, return 1-abs(correlation); otherwise return correlation.
        integer, intent(in), optional :: method !! Integration selector passed to functional inner products.
        real(dp), allocatable :: inner(:, :)
        real(dp), allocatable :: n1(:)
        real(dp), allocatable :: n2(:)
        real(dp) :: denom
        logical :: distance
        integer :: i
        integer :: j
        integer :: meth

        distance = .false.
        if (present(as_distance)) distance = as_distance
        meth = 2
        if (present(method)) meth = method
        allocate(inner(size(curves1, 1), size(curves2, 1)))
        allocate(n1(size(curves1, 1)), n2(size(curves2, 1)))
        call inprod_curves(x, curves1, curves2, inner, method=meth)
        call norm_curves(x, curves1, n1, method=meth)
        call norm_curves(x, curves2, n2, method=meth)
        do i = 1, size(curves1, 1)
            do j = 1, size(curves2, 1)
                denom = n1(i)*n2(j)
                if (denom <= tiny(1.0_dp)) then
                    rho(i, j) = 0.0_dp
                else
                    rho(i, j) = inner(i, j)/denom
                end if
                if (distance) rho(i, j) = 1.0_dp - abs(rho(i, j))
            end do
        end do
    end subroutine cosine_proximity

    pure logical function grid_is_equidistant(x) result(equi)
        real(dp), intent(in) :: x(:) !! Grid values tested for constant spacing within roundoff-scaled tolerance.
        integer :: i
        real(dp) :: h
        real(dp) :: tol

        if (size(x) <= 2) then
            equi = .true.
            return
        end if
        h = (x(size(x)) - x(1))/real(size(x) - 1, dp)
        tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(h), maxval(abs(x)))
        equi = .true.
        do i = 2, size(x)
            if (abs((x(i) - x(i - 1)) - h) > tol) then
                equi = .false.
                exit
            end if
        end do
    end function grid_is_equidistant

    pure real(dp) function composite_simpson_equal(x, y) result(value)
        real(dp), intent(in) :: x(:) !! Equally spaced grid with an odd number of points.
        real(dp), intent(in) :: y(:) !! Function values on x.
        integer :: i
        real(dp) :: h

        if (size(x) < 2) then
            value = 0.0_dp
            return
        end if
        if (size(x) == 2) then
            value = 0.5_dp*(x(2) - x(1))*(y(1) + y(2))
            return
        end if
        h = (x(size(x)) - x(1))/real(size(x) - 1, dp)
        value = y(1) + y(size(y))
        do i = 2, size(y) - 1
            if (mod(i, 2) == 0) then
                value = value + 4.0_dp*y(i)
            else
                value = value + 2.0_dp*y(i)
            end if
        end do
        value = h*value/3.0_dp
    end function composite_simpson_equal

    pure real(dp) function extended_simpson_equal(x, y) result(value)
        real(dp), intent(in) :: x(:) !! Equally spaced grid with at least five points.
        real(dp), intent(in) :: y(:) !! Function values on x.
        integer :: n
        real(dp) :: h

        n = size(x)
        if (n <= 4) then
            value = composite_simpson_equal(x, y)
            return
        end if
        h = (x(n) - x(1))/real(n - 1, dp)
        value = 17.0_dp*(y(1) + y(n)) + 59.0_dp*(y(2) + y(n - 1)) &
            + 43.0_dp*(y(3) + y(n - 2)) + 49.0_dp*(y(4) + y(n - 3))
        if (n > 8) value = value + 48.0_dp*sum(y(5:n - 4))
        value = h*value/48.0_dp
    end function extended_simpson_equal

    pure subroutine interpolate_equal_grid(x, y, nout, xx, yy)
        real(dp), intent(in) :: x(:) !! Strictly increasing source grid for linear interpolation.
        real(dp), intent(in) :: y(:) !! Source values paired with x.
        integer, intent(in) :: nout !! Requested number of equally spaced output points, at least two.
        real(dp), allocatable, intent(out) :: xx(:) !! Equally spaced output grid spanning x(1) to x(end).
        real(dp), allocatable, intent(out) :: yy(:) !! Linearly interpolated values on xx.
        integer :: i
        integer :: j
        real(dp) :: t

        allocate(xx(nout), yy(nout))
        do i = 1, nout
            xx(i) = x(1) + real(i - 1, dp)*(x(size(x)) - x(1))/real(nout - 1, dp)
        end do
        j = 1
        do i = 1, nout
            do while (j < size(x) - 1 .and. xx(i) > x(j + 1))
                j = j + 1
            end do
            if (xx(i) <= x(1)) then
                yy(i) = y(1)
            else if (xx(i) >= x(size(x))) then
                yy(i) = y(size(y))
            else
                t = (xx(i) - x(j))/(x(j + 1) - x(j))
                yy(i) = (1.0_dp - t)*y(j) + t*y(j + 1)
            end if
        end do
    end subroutine interpolate_equal_grid

end module fdausc_integration
