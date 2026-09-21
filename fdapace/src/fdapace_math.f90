module fdapace_math
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use fdapace_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: bilinear_interp
    public :: factorial_real
    public :: fill_missing_linear
    public :: fmm_spline_eval
    public :: gaussian_solve
    public :: jacobi_eigen_sym
    public :: linear_interp
    public :: mean_value
    public :: normal_random
    public :: pi_dp
    public :: sample_sd
    public :: seed_rng
    public :: sort_eigen_descending
    public :: sort_real
    public :: thin_svd
    public :: weighted_poly_fit_1d
    public :: weighted_poly_fit_2d

contains

    pure logical function increasing(x)
        real(dp), intent(in) :: x(:) !! Grid to test for strict monotone increase.
        integer :: i

        increasing = .true.
        do i = 2, size(x)
            if (x(i) <= x(i - 1)) then
                increasing = .false.
                return
            end if
        end do
    end function increasing

    pure real(dp) function factorial_real(n) result(value)
        integer, intent(in) :: n !! Nonnegative integer whose factorial is required.
        integer :: i

        value = 1.0_dp
        do i = 2, n
            value = value * real(i, dp)
        end do
    end function factorial_real

    pure real(dp) function mean_value(x) result(value)
        real(dp), intent(in) :: x(:) !! Finite sample values to average.

        if (size(x) == 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            value = sum(x) / real(size(x), dp)
        end if
    end function mean_value

    pure real(dp) function sample_sd(x) result(value)
        real(dp), intent(in) :: x(:) !! Finite sample values for the n-1 standard deviation.
        real(dp) :: mu

        if (size(x) < 2) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        mu = mean_value(x)
        value = sqrt(sum((x - mu)**2) / real(size(x) - 1, dp))
    end function sample_sd

    pure subroutine sort_real(x)
        real(dp), intent(inout) :: x(:) !! Values sorted into ascending order in place.
        integer :: i
        integer :: j
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

    pure real(dp) function linear_interp(x, y, xout, extrapolate_constant) result(value)
        real(dp), intent(in) :: x(:) !! Strictly increasing source grid.
        real(dp), intent(in) :: y(:) !! Source values corresponding one-to-one with x.
        real(dp), intent(in) :: xout !! Target coordinate for interpolation.
        logical, intent(in), optional :: extrapolate_constant !! If true, clamp outside values to the nearest endpoint.
        logical :: clamp
        integer :: lo
        integer :: hi
        integer :: mid
        real(dp) :: w

        value = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(x) /= size(y) .or. size(x) == 0) return
        if (size(x) == 1) then
            value = y(1)
            return
        end if
        if (.not. increasing(x)) return

        clamp = .false.
        if (present(extrapolate_constant)) clamp = extrapolate_constant
        if (xout <= x(1)) then
            if (xout == x(1) .or. clamp) value = y(1)
            return
        end if
        if (xout >= x(size(x))) then
            if (xout == x(size(x)) .or. clamp) value = y(size(y))
            return
        end if

        lo = 1
        hi = size(x)
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x(mid) <= xout) then
                lo = mid
            else
                hi = mid
            end if
        end do
        w = (xout - x(lo)) / (x(hi) - x(lo))
        value = (1.0_dp - w) * y(lo) + w * y(hi)
    end function linear_interp

    pure real(dp) function bilinear_interp(x, y, z, xout, yout) result(value)
        real(dp), intent(in) :: x(:) !! Strictly increasing first source grid.
        real(dp), intent(in) :: y(:) !! Strictly increasing second source grid.
        real(dp), intent(in) :: z(:,:) !! Surface values with shape (size(x), size(y)).
        real(dp), intent(in) :: xout !! Target coordinate along the first grid.
        real(dp), intent(in) :: yout !! Target coordinate along the second grid.
        integer :: ix
        integer :: iy
        real(dp) :: tx
        real(dp) :: ty

        value = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(z, 1) /= size(x) .or. size(z, 2) /= size(y)) return
        if (size(x) < 2 .or. size(y) < 2) return
        if (.not. increasing(x) .or. .not. increasing(y)) return
        if (xout < x(1) .or. xout > x(size(x))) return
        if (yout < y(1) .or. yout > y(size(y))) return

        ix = 1
        do while (ix < size(x) - 1 .and. x(ix + 1) < xout)
            ix = ix + 1
        end do
        iy = 1
        do while (iy < size(y) - 1 .and. y(iy + 1) < yout)
            iy = iy + 1
        end do
        tx = (xout - x(ix)) / (x(ix + 1) - x(ix))
        ty = (yout - y(iy)) / (y(iy + 1) - y(iy))
        value = (1.0_dp - tx) * (1.0_dp - ty) * z(ix, iy) &
                + tx * (1.0_dp - ty) * z(ix + 1, iy) &
                + (1.0_dp - tx) * ty * z(ix, iy + 1) &
                + tx * ty * z(ix + 1, iy + 1)
    end function bilinear_interp

    subroutine gaussian_solve(a, b, x, info)
        real(dp), intent(in) :: a(:,:) !! Square coefficient matrix.
        real(dp), intent(in) :: b(:) !! Right-hand-side vector with size equal to the matrix order.
        real(dp), intent(out) :: x(:) !! Solution vector; NaNs are returned when the system is singular.
        integer, intent(out) :: info !! Zero on success; nonzero when dimensions are invalid or the matrix is singular.
        real(dp), allocatable :: aa(:,:)
        real(dp), allocatable :: bb(:)
        real(dp) :: factor
        real(dp) :: pivot_abs
        real(dp) :: temp
        real(dp) :: tol
        integer :: i
        integer :: j
        integer :: k
        integer :: n
        integer :: p

        n = size(b)
        info = 0
        x = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(a, 1) /= n .or. size(a, 2) /= n .or. size(x) /= n) then
            info = -1
            return
        end if
        if (n == 0) return

        allocate(aa(n, n), bb(n))
        aa = a
        bb = b
        tol = epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aa))) * real(max(1, n), dp)

        do k = 1, n - 1
            p = k
            pivot_abs = abs(aa(k, k))
            do i = k + 1, n
                if (abs(aa(i, k)) > pivot_abs) then
                    p = i
                    pivot_abs = abs(aa(i, k))
                end if
            end do
            if (pivot_abs <= tol) then
                info = k
                return
            end if
            if (p /= k) then
                do j = k, n
                    temp = aa(k, j)
                    aa(k, j) = aa(p, j)
                    aa(p, j) = temp
                end do
                temp = bb(k)
                bb(k) = bb(p)
                bb(p) = temp
            end if
            do i = k + 1, n
                factor = aa(i, k) / aa(k, k)
                aa(i, k) = 0.0_dp
                do j = k + 1, n
                    aa(i, j) = aa(i, j) - factor * aa(k, j)
                end do
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do
        if (abs(aa(n, n)) <= tol) then
            info = n
            return
        end if

        x(n) = bb(n) / aa(n, n)
        do i = n - 1, 1, -1
            x(i) = (bb(i) - dot_product(aa(i, i + 1:n), x(i + 1:n))) / aa(i, i)
        end do
    end subroutine gaussian_solve

    subroutine jacobi_eigen_sym(a, values, vectors, info)
        real(dp), intent(in) :: a(:,:) !! Real symmetric matrix whose eigensystem is required.
        real(dp), intent(out) :: values(:) !! Eigenvalues, initially in unsorted order.
        real(dp), intent(out) :: vectors(:,:) !! Orthogonal eigenvectors stored as columns.
        integer, intent(out) :: info !! Zero on convergence; positive when the iteration limit is reached; negative on shape errors.
        real(dp), allocatable :: work(:,:)
        real(dp) :: app
        real(dp) :: apq
        real(dp) :: aqq
        real(dp) :: c
        real(dp) :: phi
        real(dp) :: s
        real(dp) :: t1
        real(dp) :: t2
        real(dp) :: tol
        integer :: i
        integer :: iter
        integer :: max_iter
        integer :: n
        integer :: p
        integer :: q

        n = size(a, 1)
        info = 0
        if (size(a, 2) /= n .or. size(values) /= n .or. size(vectors, 1) /= n .or. size(vectors, 2) /= n) then
            info = -1
            return
        end if
        allocate(work(n, n))
        work = 0.5_dp * (a + transpose(a))
        vectors = 0.0_dp
        do i = 1, n
            vectors(i, i) = 1.0_dp
        end do
        if (n <= 1) then
            if (n == 1) values(1) = work(1, 1)
            return
        end if

        tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(work)))
        max_iter = max(100, 50 * n * n)
        do iter = 1, max_iter
            p = 1
            q = 2
            apq = abs(work(p, q))
            do i = 1, n
                if (i >= n) cycle
                do q = i + 1, n
                    if (abs(work(i, q)) > apq) then
                        p = i
                        apq = abs(work(i, q))
                    end if
                end do
            end do
            q = p + 1
            apq = abs(work(p, q))
            do i = p + 1, n
                if (abs(work(p, i)) > apq) then
                    q = i
                    apq = abs(work(p, i))
                end if
            end do
            do i = 1, n - 1
                if (i == p) cycle
                if (maxval(abs(work(i, i + 1:n))) > apq) then
                    q = maxloc(abs(work(i, i + 1:n)), dim=1) + i
                    p = i
                    apq = abs(work(p, q))
                end if
            end do
            if (apq <= tol) exit

            app = work(p, p)
            aqq = work(q, q)
            phi = 0.5_dp * atan2(2.0_dp * work(p, q), aqq - app)
            c = cos(phi)
            s = sin(phi)

            do i = 1, n
                if (i /= p .and. i /= q) then
                    t1 = work(i, p)
                    t2 = work(i, q)
                    work(i, p) = c * t1 - s * t2
                    work(p, i) = work(i, p)
                    work(i, q) = s * t1 + c * t2
                    work(q, i) = work(i, q)
                end if
            end do
            work(p, p) = c * c * app - 2.0_dp * s * c * work(p, q) + s * s * aqq
            work(q, q) = s * s * app + 2.0_dp * s * c * work(p, q) + c * c * aqq
            work(p, q) = 0.0_dp
            work(q, p) = 0.0_dp

            do i = 1, n
                t1 = vectors(i, p)
                t2 = vectors(i, q)
                vectors(i, p) = c * t1 - s * t2
                vectors(i, q) = s * t1 + c * t2
            end do
        end do
        if (iter > max_iter) info = 1
        do i = 1, n
            values(i) = work(i, i)
        end do
    end subroutine jacobi_eigen_sym

    subroutine sort_eigen_descending(values, vectors)
        real(dp), intent(inout) :: values(:) !! Eigenvalues sorted into descending order.
        real(dp), intent(inout) :: vectors(:,:) !! Matching eigenvectors reordered by columns.
        integer :: i
        integer :: j
        integer :: best
        real(dp) :: tmp
        real(dp), allocatable :: col(:)

        allocate(col(size(vectors, 1)))
        do i = 1, size(values) - 1
            best = i
            do j = i + 1, size(values)
                if (values(j) > values(best)) best = j
            end do
            if (best /= i) then
                tmp = values(i)
                values(i) = values(best)
                values(best) = tmp
                col = vectors(:, i)
                vectors(:, i) = vectors(:, best)
                vectors(:, best) = col
            end if
        end do
    end subroutine sort_eigen_descending

    subroutine thin_svd(a, u, s, v, info)
        real(dp), intent(in) :: a(:,:) !! Input matrix to decompose.
        real(dp), allocatable, intent(out) :: u(:,:) !! Left singular vectors corresponding to returned singular values.
        real(dp), allocatable, intent(out) :: s(:) !! Singular values in descending order.
        real(dp), allocatable, intent(out) :: v(:,:) !! Right singular vectors corresponding to returned singular values.
        integer, intent(out) :: info !! Zero on success; nonzero when the symmetric eigensolver fails.
        real(dp), allocatable :: ata(:,:)
        real(dp), allocatable :: eval(:)
        real(dp), allocatable :: evec(:,:)
        real(dp) :: normu
        real(dp) :: tol
        integer :: i
        integer :: m
        integer :: n
        integer :: r

        m = size(a, 1)
        n = size(a, 2)
        r = min(m, n)
        allocate(ata(n, n), eval(n), evec(n, n))
        ata = matmul(transpose(a), a)
        call jacobi_eigen_sym(ata, eval, evec, info)
        if (info /= 0) then
            allocate(u(m, 0), s(0), v(n, 0))
            return
        end if
        call sort_eigen_descending(eval, evec)
        allocate(u(m, r), s(r), v(n, r))
        s = sqrt(max(eval(1:r), 0.0_dp))
        v = evec(:, 1:r)
        tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(s))
        do i = 1, r
            if (s(i) > tol) then
                u(:, i) = matmul(a, v(:, i)) / s(i)
                normu = sqrt(sum(u(:, i)**2))
                if (normu > 0.0_dp) u(:, i) = u(:, i) / normu
            else
                u(:, i) = 0.0_dp
            end if
        end do
    end subroutine thin_svd

    subroutine weighted_poly_fit_1d(x, y, w, x0, degree, beta, info)
        real(dp), intent(in) :: x(:) !! Abscissae used by the local polynomial fit.
        real(dp), intent(in) :: y(:) !! Responses corresponding to x.
        real(dp), intent(in) :: w(:) !! Nonnegative observation weights corresponding to x.
        real(dp), intent(in) :: x0 !! Local polynomial expansion point.
        integer, intent(in) :: degree !! Polynomial degree, zero or greater.
        real(dp), intent(out) :: beta(:) !! Coefficients for powers of (x-x0), from degree zero upward.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or a singular normal system.
        real(dp), allocatable :: normal(:,:)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: row(:)
        integer :: i
        integer :: j
        integer :: k
        integer :: p

        p = degree + 1
        info = 0
        if (size(x) /= size(y) .or. size(x) /= size(w) .or. size(beta) /= p .or. degree < 0) then
            info = -1
            return
        end if
        allocate(normal(p, p), rhs(p), row(p))
        normal = 0.0_dp
        rhs = 0.0_dp
        do i = 1, size(x)
            if (w(i) <= 0.0_dp) cycle
            row(1) = 1.0_dp
            do j = 2, p
                row(j) = row(j - 1) * (x(i) - x0)
            end do
            do j = 1, p
                rhs(j) = rhs(j) + w(i) * row(j) * y(i)
                do k = 1, p
                    normal(j, k) = normal(j, k) + w(i) * row(j) * row(k)
                end do
            end do
        end do
        call gaussian_solve(normal, rhs, beta, info)
    end subroutine weighted_poly_fit_1d

    subroutine weighted_poly_fit_2d(x, y, z, w, x0, y0, degree, beta, info)
        real(dp), intent(in) :: x(:) !! First coordinates used by the local polynomial fit.
        real(dp), intent(in) :: y(:) !! Second coordinates used by the local polynomial fit.
        real(dp), intent(in) :: z(:) !! Responses corresponding to each coordinate pair.
        real(dp), intent(in) :: w(:) !! Nonnegative observation weights corresponding to coordinate pairs.
        real(dp), intent(in) :: x0 !! Local polynomial expansion point for the first coordinate.
        real(dp), intent(in) :: y0 !! Local polynomial expansion point for the second coordinate.
        integer, intent(in) :: degree !! Total polynomial degree, zero or greater.
        real(dp), intent(out) :: beta(:) !! Coefficients ordered by total degree and then second-coordinate degree.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or a singular normal system.
        real(dp), allocatable :: normal(:,:)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: row(:)
        integer :: col
        integer :: d2
        integer :: i
        integer :: j
        integer :: k
        integer :: p
        integer :: total

        p = (degree + 1) * (degree + 2) / 2
        info = 0
        if (size(x) /= size(y) .or. size(x) /= size(z) .or. size(x) /= size(w) .or. size(beta) /= p .or. degree < 0) then
            info = -1
            return
        end if
        allocate(normal(p, p), rhs(p), row(p))
        normal = 0.0_dp
        rhs = 0.0_dp
        do i = 1, size(x)
            if (w(i) <= 0.0_dp) cycle
            col = 0
            do total = 0, degree
                do d2 = 0, total
                    col = col + 1
                    row(col) = (x(i) - x0)**(total - d2) * (y(i) - y0)**d2
                end do
            end do
            do j = 1, p
                rhs(j) = rhs(j) + w(i) * row(j) * z(i)
                do k = 1, p
                    normal(j, k) = normal(j, k) + w(i) * row(j) * row(k)
                end do
            end do
        end do
        call gaussian_solve(normal, rhs, beta, info)
    end subroutine weighted_poly_fit_2d

    subroutine fill_missing_linear(t, y, filled, info)
        real(dp), intent(in) :: t(:) !! Strictly increasing observation grid.
        real(dp), intent(in) :: y(:) !! Observations; NaNs are treated as missing values.
        real(dp), intent(out) :: filled(:) !! Linearly interpolated values with constant endpoint extrapolation.
        integer, intent(out) :: info !! Zero on success; nonzero if dimensions are invalid or all observations are missing.
        real(dp), allocatable :: tx(:)
        real(dp), allocatable :: yy(:)
        integer :: i
        integer :: nvalid
        integer :: pos

        info = 0
        if (size(t) /= size(y) .or. size(filled) /= size(y) .or. .not. increasing(t)) then
            info = -1
            return
        end if
        nvalid = count(.not. ieee_is_nan(y))
        if (nvalid == 0) then
            filled = ieee_value(0.0_dp, ieee_quiet_nan)
            info = 1
            return
        end if
        allocate(tx(nvalid), yy(nvalid))
        pos = 0
        do i = 1, size(y)
            if (.not. ieee_is_nan(y(i))) then
                pos = pos + 1
                tx(pos) = t(i)
                yy(pos) = y(i)
            end if
        end do
        do i = 1, size(y)
            filled(i) = linear_interp(tx, yy, t(i), .true.)
        end do
    end subroutine fill_missing_linear

    subroutine fmm_spline_eval(x, y, xout, yout, info)
        real(dp), intent(in) :: x(:) !! Strictly increasing spline knot coordinates.
        real(dp), intent(in) :: y(:) !! Spline values at the knots.
        real(dp), intent(in) :: xout(:) !! Coordinates at which to evaluate the FMM/not-a-knot cubic spline.
        real(dp), intent(out) :: yout(:) !! Interpolated spline values corresponding to xout.
        integer, intent(out) :: info !! Zero on success; nonzero on invalid input or failed spline system solution.
        real(dp), allocatable :: mat(:,:)
        real(dp), allocatable :: rhs(:)
        real(dp), allocatable :: m2(:)
        real(dp) :: a
        real(dp) :: b
        real(dp) :: h
        integer :: i
        integer :: j
        integer :: n

        n = size(x)
        info = 0
        if (n /= size(y) .or. size(xout) /= size(yout) .or. n < 2 .or. .not. increasing(x)) then
            info = -1
            return
        end if
        if (minval(xout) < x(1) .or. maxval(xout) > x(n)) then
            info = -2
            return
        end if
        if (n == 2) then
            do i = 1, size(xout)
                yout(i) = linear_interp(x, y, xout(i))
            end do
            return
        end if
        if (n == 3) then
            do i = 1, size(xout)
                yout(i) = lagrange3(x, y, xout(i))
            end do
            return
        end if

        allocate(mat(n, n), rhs(n), m2(n))
        mat = 0.0_dp
        rhs = 0.0_dp
        mat(1, 1) = -(x(3) - x(2))
        mat(1, 2) = x(3) - x(1)
        mat(1, 3) = -(x(2) - x(1))
        do i = 2, n - 1
            mat(i, i - 1) = x(i) - x(i - 1)
            mat(i, i) = 2.0_dp * (x(i + 1) - x(i - 1))
            mat(i, i + 1) = x(i + 1) - x(i)
            rhs(i) = 6.0_dp * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) &
                     - (y(i) - y(i - 1)) / (x(i) - x(i - 1)))
        end do
        mat(n, n - 2) = -(x(n) - x(n - 1))
        mat(n, n - 1) = x(n) - x(n - 2)
        mat(n, n) = -(x(n - 1) - x(n - 2))
        call gaussian_solve(mat, rhs, m2, info)
        if (info /= 0) return

        do i = 1, size(xout)
            if (xout(i) == x(n)) then
                j = n - 1
            else
                j = 1
                do while (j < n - 1 .and. x(j + 1) <= xout(i))
                    j = j + 1
                end do
            end if
            h = x(j + 1) - x(j)
            a = (x(j + 1) - xout(i)) / h
            b = (xout(i) - x(j)) / h
            yout(i) = a * y(j) + b * y(j + 1) &
                      + ((a**3 - a) * m2(j) + (b**3 - b) * m2(j + 1)) * h * h / 6.0_dp
        end do
    end subroutine fmm_spline_eval

    pure real(dp) function lagrange3(x, y, xout) result(value)
        real(dp), intent(in) :: x(3) !! Three distinct interpolation coordinates.
        real(dp), intent(in) :: y(3) !! Function values at the three interpolation coordinates.
        real(dp), intent(in) :: xout !! Coordinate at which the quadratic interpolant is evaluated.

        value = y(1) * (xout - x(2)) * (xout - x(3)) / ((x(1) - x(2)) * (x(1) - x(3))) &
                + y(2) * (xout - x(1)) * (xout - x(3)) / ((x(2) - x(1)) * (x(2) - x(3))) &
                + y(3) * (xout - x(1)) * (xout - x(2)) / ((x(3) - x(1)) * (x(3) - x(2)))
    end function lagrange3

    subroutine seed_rng(seed)
        integer, intent(in) :: seed !! Deterministic scalar seed used to initialize the processor's Fortran RNG state.
        integer, allocatable :: put(:)
        integer :: i
        integer :: n
        integer :: modulus

        call random_seed(size=n)
        allocate(put(n))
        modulus = huge(1) - 1
        do i = 1, n
            put(i) = modulo(abs(seed) + 104729 * i + 37 * i * i, modulus)
            if (put(i) == 0) put(i) = i
        end do
        call random_seed(put=put)
    end subroutine seed_rng

    real(dp) function normal_random() result(value)
        real(dp) :: u1
        real(dp) :: u2

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi_dp * u2)
    end function normal_random

end module fdapace_math
