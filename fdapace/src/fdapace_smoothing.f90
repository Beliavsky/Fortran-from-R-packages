module fdapace_smoothing
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use fdapace_kinds, only : dp
    use fdapace_math, only : bilinear_interp, factorial_real, linear_interp, pi_dp, &
                             weighted_poly_fit_1d, weighted_poly_fit_2d
    implicit none
    private

    public :: convert_support
    public :: convert_support_covariance
    public :: convert_support_functions
    public :: cumtrapz_rcpp
    public :: lwls1d
    public :: lwls2d
    public :: lwls2d_deriv
    public :: norm_curv_to_area
    public :: trapz_rcpp

contains

    pure real(dp) function trapz_rcpp(x, y) result(value)
        real(dp), intent(in) :: x(:) !! Sorted integration coordinates.
        real(dp), intent(in) :: y(:) !! Function values corresponding one-to-one with x.
        integer :: i

        if (size(x) /= size(y)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        value = 0.0_dp
        do i = 1, size(x) - 1
            value = value + 0.5_dp * (x(i + 1) - x(i)) * (y(i) + y(i + 1))
        end do
    end function trapz_rcpp

    pure function cumtrapz_rcpp(x, y) result(values)
        real(dp), intent(in) :: x(:) !! Sorted integration coordinates.
        real(dp), intent(in) :: y(:) !! Function values corresponding one-to-one with x.
        real(dp), allocatable :: values(:)
        integer :: i

        allocate(values(size(x)))
        if (size(x) /= size(y)) then
            values = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (size(values) == 0) return
        values(1) = 0.0_dp
        do i = 2, size(x)
            values(i) = values(i - 1) + 0.5_dp * (x(i) - x(i - 1)) * (y(i) + y(i - 1))
        end do
    end function cumtrapz_rcpp

    pure function norm_curv_to_area(y, x, area) result(values)
        real(dp), intent(in) :: y(:) !! Function values to rescale.
        real(dp), intent(in) :: x(:) !! Integration coordinates corresponding one-to-one with y.
        real(dp), intent(in), optional :: area !! Desired trapezoidal area; defaults to one.
        real(dp), allocatable :: values(:)
        real(dp) :: target
        real(dp) :: current

        target = 1.0_dp
        if (present(area)) target = area
        allocate(values(size(y)))
        current = trapz_rcpp(x, y)
        if (current == 0.0_dp .or. ieee_is_nan(current)) then
            values = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            values = target * y / current
        end if
    end function norm_curv_to_area

    pure function convert_support(from_grid, to_grid, mu) result(values)
        real(dp), intent(in) :: from_grid(:) !! Sorted source support coordinates.
        real(dp), intent(in) :: to_grid(:) !! Sorted target support coordinates contained within from_grid.
        real(dp), intent(in) :: mu(:) !! Source function values corresponding one-to-one with from_grid.
        real(dp), allocatable :: values(:)
        integer :: i

        call check_support(from_grid, to_grid)
        if (size(mu) /= size(from_grid)) error stop "convert_support: incompatible vector length"
        allocate(values(size(to_grid)))
        do i = 1, size(to_grid)
            values(i) = linear_interp(from_grid, mu, to_grid(i))
        end do
    end function convert_support

    pure function convert_support_functions(from_grid, to_grid, phi) result(values)
        real(dp), intent(in) :: from_grid(:) !! Sorted source support coordinates.
        real(dp), intent(in) :: to_grid(:) !! Sorted target support coordinates contained within from_grid.
        real(dp), intent(in) :: phi(:,:) !! Source functions in columns, with rows corresponding to from_grid.
        real(dp), allocatable :: values(:,:)
        integer :: i
        integer :: j

        call check_support(from_grid, to_grid)
        if (size(phi, 1) /= size(from_grid)) error stop "convert_support_functions: incompatible matrix shape"
        allocate(values(size(to_grid), size(phi, 2)))
        do j = 1, size(phi, 2)
            do i = 1, size(to_grid)
                values(i, j) = linear_interp(from_grid, phi(:, j), to_grid(i))
            end do
        end do
    end function convert_support_functions

    pure function convert_support_covariance(from_grid, to_grid, cov, is_cross_cov) result(values)
        real(dp), intent(in) :: from_grid(:) !! Sorted source support coordinates for both covariance dimensions.
        real(dp), intent(in) :: to_grid(:) !! Sorted target support coordinates contained within from_grid.
        real(dp), intent(in) :: cov(:,:) !! Source covariance or cross-covariance surface on from_grid by from_grid.
        logical, intent(in), optional :: is_cross_cov !! If true, retain asymmetry instead of symmetrizing the result.
        real(dp), allocatable :: values(:,:)
        logical :: cross
        integer :: i
        integer :: j

        call check_support(from_grid, to_grid)
        if (size(cov, 1) /= size(from_grid) .or. size(cov, 2) /= size(from_grid)) then
            error stop "convert_support_covariance: incompatible covariance shape"
        end if
        cross = .false.
        if (present(is_cross_cov)) cross = is_cross_cov
        allocate(values(size(to_grid), size(to_grid)))
        do j = 1, size(to_grid)
            do i = 1, size(to_grid)
                values(i, j) = bilinear_interp(from_grid, from_grid, cov, to_grid(i), to_grid(j))
            end do
        end do
        if (.not. cross) values = 0.5_dp * (values + transpose(values))
    end function convert_support_covariance

    pure subroutine check_support(from_grid, to_grid)
        real(dp), intent(in) :: from_grid(:) !! Source support used for interpolation-range validation.
        real(dp), intent(in) :: to_grid(:) !! Target support used for interpolation-range validation.
        real(dp) :: buff

        if (size(from_grid) < 2 .or. size(to_grid) < 1) error stop "convert_support: empty support"
        buff = epsilon(1.0_dp) * max(1.0_dp, maxval(abs(from_grid))) * 3.0_dp
        if (minval(to_grid) < from_grid(1) - buff .or. maxval(to_grid) > from_grid(size(from_grid)) + buff) then
            error stop "convert_support: insufficient source support"
        end if
    end subroutine check_support

    function lwls1d(bw, kernel_type, xin, yin, xout, win, npoly, nder) result(values)
        real(dp), intent(in) :: bw !! Positive local-smoothing bandwidth.
        character(len=*), intent(in) :: kernel_type !! Kernel name: epan, rect, gauss, gausvar, or quar.
        real(dp), intent(in) :: xin(:) !! Increasing measurement coordinates.
        real(dp), intent(in) :: yin(:) !! Measurement values corresponding to xin; NaNs are omitted.
        real(dp), intent(in) :: xout(:) !! Increasing coordinates at which to evaluate the smoother.
        real(dp), intent(in), optional :: win(:) !! Optional observation weights; defaults to one and NaNs are omitted.
        integer, intent(in), optional :: npoly !! Local polynomial degree; defaults to one.
        integer, intent(in), optional :: nder !! Derivative order to estimate; defaults to zero.
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: xuse(:)
        real(dp), allocatable :: yuse(:)
        real(dp), allocatable :: wuse(:)
        real(dp) :: u
        integer :: degree
        integer :: deriv
        integer :: i
        integer :: info
        integer :: j
        integer :: nvalid
        integer :: pos

        if (bw <= 0.0_dp) error stop "lwls1d: bandwidth must be positive"
        if (size(xin) /= size(yin)) error stop "lwls1d: xin and yin lengths differ"
        if (present(win)) then
            if (size(win) /= size(xin)) error stop "lwls1d: win length differs"
        end if
        do i = 2, size(xin)
            if (xin(i) < xin(i - 1)) error stop "lwls1d: xin must be sorted"
        end do
        do i = 2, size(xout)
            if (xout(i) < xout(i - 1)) error stop "lwls1d: xout must be sorted"
        end do
        degree = 1
        if (present(npoly)) degree = npoly
        deriv = 0
        if (present(nder)) deriv = nder
        if (degree < 0 .or. deriv < 0 .or. deriv > degree) error stop "lwls1d: invalid polynomial or derivative degree"

        nvalid = 0
        do i = 1, size(xin)
            if (ieee_is_nan(xin(i)) .or. ieee_is_nan(yin(i))) cycle
            if (present(win)) then
                if (ieee_is_nan(win(i))) cycle
            end if
            nvalid = nvalid + 1
        end do
        allocate(xuse(nvalid), yuse(nvalid), wuse(nvalid))
        pos = 0
        do i = 1, size(xin)
            if (ieee_is_nan(xin(i)) .or. ieee_is_nan(yin(i))) cycle
            if (present(win)) then
                if (ieee_is_nan(win(i))) cycle
            end if
            pos = pos + 1
            xuse(pos) = xin(i)
            yuse(pos) = yin(i)
            if (present(win)) then
                wuse(pos) = win(i)
            else
                wuse(pos) = 1.0_dp
            end if
        end do
        allocate(values(size(xout)), beta(degree + 1), weights(nvalid))
        values = ieee_value(0.0_dp, ieee_quiet_nan)
        do i = 1, size(xout)
            do j = 1, nvalid
                u = (xuse(j) - xout(i)) / bw
                weights(j) = wuse(j) * kernel_1d(u, kernel_type)
                if (trim(kernel_type) /= "gauss" .and. trim(kernel_type) /= "gausvar") then
                    if (xuse(j) < xout(i) - bw .or. xuse(j) >= xout(i) + bw) weights(j) = 0.0_dp
                end if
            end do
            if (count(weights /= 0.0_dp) <= deriv) cycle
            call weighted_poly_fit_1d(xuse, yuse, weights, xout(i), degree, beta, info)
            if (info == 0) values(i) = beta(deriv + 1) * factorial_real(deriv)
        end do
    end function lwls1d

    function lwls2d(bw, kernel_type, xin, yin, xout1, xout2, win, crosscov) result(values)
        real(dp), intent(in) :: bw(:) !! One or two positive bandwidths; a scalar-equivalent size-one array is recycled.
        character(len=*), intent(in) :: kernel_type !! Kernel name: epan, rect, gauss, gausvar, or quar.
        real(dp), intent(in) :: xin(:,:) !! Measurement coordinates with shape (n,2).
        real(dp), intent(in) :: yin(:) !! Measurement values corresponding to xin rows.
        real(dp), intent(in) :: xout1(:) !! First output-coordinate grid.
        real(dp), intent(in) :: xout2(:) !! Second output-coordinate grid.
        real(dp), intent(in), optional :: win(:) !! Optional observation weights; defaults to one.
        logical, intent(in), optional :: crosscov !! True for a general/cross-covariance surface; false requests symmetric autocovariance output.
        real(dp), allocatable :: values(:,:)
        logical :: is_cross

        is_cross = .false.
        if (present(crosscov)) is_cross = crosscov
        values = lwls2d_deriv(bw, kernel_type, xin, yin, xout1, xout2, win, 1, 0, 0, .not. is_cross)
    end function lwls2d

    function lwls2d_deriv(bw, kernel_type, xin, yin, xout1, xout2, win, npoly, nder1, nder2, auto_cov) result(values)
        real(dp), intent(in) :: bw(:) !! One or two positive bandwidths; a size-one array is recycled across dimensions.
        character(len=*), intent(in) :: kernel_type !! Kernel name: epan, rect, gauss, gausvar, or quar.
        real(dp), intent(in) :: xin(:,:) !! Measurement coordinates with shape (n,2).
        real(dp), intent(in) :: yin(:) !! Measurement values corresponding to xin rows.
        real(dp), intent(in) :: xout1(:) !! First output-coordinate grid.
        real(dp), intent(in) :: xout2(:) !! Second output-coordinate grid.
        real(dp), intent(in), optional :: win(:) !! Optional observation weights; defaults to one.
        integer, intent(in), optional :: npoly !! Total local polynomial degree; defaults to one.
        integer, intent(in), optional :: nder1 !! Derivative order in the first direction; defaults to zero.
        integer, intent(in), optional :: nder2 !! Derivative order in the second direction; defaults to zero.
        logical, intent(in), optional :: auto_cov !! If true and derivative orders match, enforce a symmetric output surface.
        real(dp), allocatable :: values(:,:)
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: base_weights(:)
        real(dp) :: bw1
        real(dp) :: bw2
        real(dp) :: u
        real(dp) :: v
        integer :: degree
        integer :: d1
        integer :: d2
        integer :: i
        integer :: idx
        integer :: info
        integer :: j
        integer :: ncoef
        integer :: p
        logical :: symmetric

        if (size(xin, 2) /= 2 .or. size(xin, 1) /= size(yin)) error stop "lwls2d_deriv: incompatible input shapes"
        if (size(bw) /= 1 .and. size(bw) /= 2) error stop "lwls2d_deriv: bw must have length one or two"
        bw1 = bw(1)
        if (size(bw) == 1) then
            bw2 = bw(1)
        else
            bw2 = bw(2)
        end if
        if (bw1 <= 0.0_dp .or. bw2 <= 0.0_dp) error stop "lwls2d_deriv: bandwidths must be positive"
        if (present(win)) then
            if (size(win) /= size(yin)) error stop "lwls2d_deriv: win length differs"
        end if
        degree = 1
        if (present(npoly)) degree = npoly
        d1 = 0
        if (present(nder1)) d1 = nder1
        d2 = 0
        if (present(nder2)) d2 = nder2
        if (degree < d1 + d2 .or. d1 < 0 .or. d2 < 0) error stop "lwls2d_deriv: invalid derivative degree"
        symmetric = .false.
        if (present(auto_cov)) symmetric = auto_cov .and. d1 == d2
        if (symmetric .and. size(xout1) /= size(xout2)) error stop "lwls2d_deriv: symmetric grids differ in size"
        if (symmetric .and. any(abs(xout1 - xout2) > 10.0_dp * epsilon(1.0_dp))) then
            error stop "lwls2d_deriv: symmetric grids differ"
        end if

        ncoef = (degree + 1) * (degree + 2) / 2
        allocate(values(size(xout1), size(xout2)), beta(ncoef), weights(size(yin)), base_weights(size(yin)))
        if (present(win)) then
            base_weights = win
        else
            base_weights = 1.0_dp
        end if
        values = ieee_value(0.0_dp, ieee_quiet_nan)
        idx = (d1 + d2) * (d1 + d2 + 1) / 2 + d2 + 1
        do j = 1, size(xout2)
            do i = 1, size(xout1)
                if (symmetric .and. j < i) cycle
                do p = 1, size(yin)
                    u = (xin(p, 1) - xout1(i)) / bw1
                    v = (xin(p, 2) - xout2(j)) / bw2
                    weights(p) = base_weights(p) * kernel_2d(u, v, kernel_type)
                    if (trim(kernel_type) /= "gauss") then
                        if (abs(xin(p, 1) - xout1(i)) > bw1 + 1.0e-6_dp &
                            .or. abs(xin(p, 2) - xout2(j)) > bw2 + 1.0e-6_dp) weights(p) = 0.0_dp
                    end if
                end do
                if (count(weights /= 0.0_dp) < ncoef) cycle
                call weighted_poly_fit_2d(xin(:, 1), xin(:, 2), yin, weights, xout1(i), xout2(j), degree, beta, info)
                if (info == 0) values(i, j) = beta(idx) * factorial_real(d1) * factorial_real(d2)
            end do
        end do
        if (symmetric) then
            do j = 1, size(xout2)
                do i = j + 1, size(xout1)
                    values(i, j) = values(j, i)
                end do
            end do
        end if
    end function lwls2d_deriv

    pure real(dp) function kernel_1d(u, kernel_type) result(value)
        real(dp), intent(in) :: u !! Bandwidth-scaled coordinate difference.
        character(len=*), intent(in) :: kernel_type !! Kernel family name.
        real(dp) :: inv_sqrt_2pi

        inv_sqrt_2pi = 1.0_dp / sqrt(2.0_dp * pi_dp)
        select case (trim(kernel_type))
        case ("epan")
            value = 0.75_dp * (1.0_dp - u * u)
        case ("rect")
            value = 1.0_dp
        case ("gauss")
            value = exp(-0.5_dp * u * u) * inv_sqrt_2pi
        case ("gausvar")
            value = exp(-0.5_dp * u * u) * inv_sqrt_2pi * (1.25_dp - 0.25_dp * u * u)
        case ("quar")
            value = (15.0_dp / 16.0_dp) * (1.0_dp - u * u)**2
        case default
            value = 0.75_dp * (1.0_dp - u * u)
        end select
    end function kernel_1d

    pure real(dp) function kernel_2d(u, v, kernel_type) result(value)
        real(dp), intent(in) :: u !! Bandwidth-scaled first coordinate difference.
        real(dp), intent(in) :: v !! Bandwidth-scaled second coordinate difference.
        character(len=*), intent(in) :: kernel_type !! Kernel family name.
        real(dp) :: inv_sqrt_2pi

        inv_sqrt_2pi = 1.0_dp / sqrt(2.0_dp * pi_dp)
        select case (trim(kernel_type))
        case ("epan")
            value = (9.0_dp / 16.0_dp) * (1.0_dp - u * u) * (1.0_dp - v * v)
        case ("rect")
            value = 0.25_dp
        case ("gauss")
            value = exp(-0.5_dp * u * u) * inv_sqrt_2pi * exp(-0.5_dp * v * v) * inv_sqrt_2pi
        case ("gausvar")
            value = exp(-0.5_dp * u * u) * inv_sqrt_2pi * exp(-0.5_dp * v * v) * inv_sqrt_2pi &
                    * (1.25_dp - 0.25_dp * u * u) * (1.50_dp - 0.50_dp * v * v)
        case ("quar")
            value = (225.0_dp / 256.0_dp) * (1.0_dp - u * u)**2 * (1.0_dp - v * v)**2
        case default
            value = (9.0_dp / 16.0_dp) * (1.0_dp - u * u) * (1.0_dp - v * v)
        end select
    end function kernel_2d

end module fdapace_smoothing
